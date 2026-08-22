-- Fit each table to the page it is printed on, and keep every table and figure
-- on the same page as the title above it.
--
-- WHAT WENT WRONG BEFORE. The first version of this filter rotated a table onto
-- a landscape page whenever it had five columns or one cell longer than 120
-- characters. Almost every table in this manual is two or three columns with one
-- prose column, and a prose column is routinely 200 to 500 characters, so the
-- rule caught thirty of them. A two-column table wraps perfectly well in
-- portrait; rotating it costs the reader a page turn and a head tilt for
-- nothing, and the \clearpage that pins the rotated page left the heading and
-- the sentence introducing the table stranded at the bottom of the page before.
--
-- WHAT ACTUALLY OVERLAPS. Not width in the abstract. pandoc sizes the columns of
-- a pipe table in proportion to how wide they are IN THE SOURCE, so a table
-- whose description column holds 500 characters and whose flag column holds 25
-- gives the flag column 4% of the text block, about 7mm. `--max-events-per-file`
-- does not fit in 7mm, and a \raggedright p-column sets what does not fit into
-- the column beside it. That is the overlap, and it happens at any page
-- orientation, so rotating was treating a symptom.
--
-- WHAT THIS DOES INSTEAD. Every column is measured: the longest cell, and the
-- longest unbreakable word in it. Each column is first given the space its
-- longest word needs, so nothing can protrude, and whatever is left over is
-- shared out among the columns holding prose, which are the ones that can use
-- it. The result is written back as explicit column widths. A table is rotated
-- only when even the longest-word minimum will not fit across a portrait page,
-- which in this manual is no table at all. The branch stays because the option
-- reference grows.
--
-- CHARACTERS, NOT MILLIMETRES. Widths are computed in characters of the body
-- face and converted to fractions of the text block at the end. One character of
-- 10pt Latin Modern is close enough to half its point size that 86 characters is
-- a good estimate of a 160mm line, and what a cell holds is something the cell
-- gives directly. Characters are not all one width, and the ones that are wider
-- than average -- monospace, capitals, digits -- are exactly what the narrow
-- columns of these tables hold, so each is weighted rather than counted.

-- Characters across the text block, per typesetting mode. Portrait is A4 less
-- 25mm margins at 10pt; landscape is the same page turned, at 8pt.
local CHARS = { plain = 86, small = 96, foot = 108, land = 172 }
-- Both gutters of a column, in characters of that mode's face: 2\tabcolsep is
-- 12pt, and pandoc's width formula subtracts it before applying the fractions.
local PAD   = { plain = 2.4, small = 2.7, foot = 3.0, land = 3.0 }
local SIZE  = { plain = nil, small = '\\small', foot = '\\footnotesize',
                land = nil }
-- Tried in order, first fit wins, so a table is shrunk before it is rotated and
-- rotated only when shrinking has run out.
local ORDER = { 'plain', 'small', 'foot', 'land' }

-- No column narrower than this, whatever its content: a 3-character column of
-- one-word cells still needs somewhere to put a word that turns out to be long.
local MIN_COL = 10

-- WHAT A CHARACTER IS WORTH. One unit is an average character of the body face.
-- A flat count is not good enough at column scale: the first version of this
-- gave a column holding cyCONDOR eight units, and eight capitals are half again
-- as wide as eight average characters, so the word was set into the column
-- beside it. Ratios are from Latin Modern at the same size.
local MONO_W  = 1.2   -- monospace, every character the same width
local UPPER_W = 1.5   -- a capital
local DIGIT_W = 1.1   -- a digit, which is set on the figure width
-- A paragraph longer than this is prose the table happens to follow, not the
-- line introducing it, so it is left where it is rather than dragged along.
local TITLE_CHARS = 320
-- Ceiling on what a title may reserve. Without it a mis-detected title asks for
-- more than a page and every table starts on a fresh one.
local MAX_RESERVE = 14

local function raw(s) return pandoc.RawBlock('latex', s) end

-- A Cell is a wrapper around a list of blocks, not a list of blocks. Passing the
-- Cell itself to stringify raises "table expected, got pandoc Cell", which
-- aborts pandoc, and with continue-on-error on the build step that failure
-- reported as success, so the only symptom was a 404 on the download link.
local function prose_width(s)
  local upper = select(2, s:gsub('%u', ''))
  local digit = select(2, s:gsub('%d', ''))
  return (#s - upper - digit) + upper * UPPER_W + digit * DIGIT_W
end

local function measure(blocks)
  local div  = pandoc.Div(blocks)
  local text = pandoc.utils.stringify(div)
  local code = {}
  pandoc.walk_block(div, { Code = function(c) code[#code + 1] = c.text end })
  -- stringify returns the code along with everything else, measured as prose.
  -- Its share is taken back out and re-measured as monospace rather than
  -- counted twice.
  local codetext = table.concat(code)
  local width = prose_width(text) - prose_width(codetext) + #codetext * MONO_W
  local token = 0
  for word in text:gmatch('%S+') do
    local w = prose_width(word)
    if w > token then token = w end
  end
  for _, c in ipairs(code) do
    for word in c:gmatch('%S+') do
      local w = #word * MONO_W
      if w > token then token = w end
    end
  end
  return width, token
end

local function col_stats(tbl)
  local n = #tbl.colspecs
  local nat, tok = {}, {}
  for i = 1, n do nat[i], tok[i] = 0, 0 end
  local function scan(rows)
    for _, row in ipairs(rows or {}) do
      local i = 1
      for _, cell in ipairs(row.cells or {}) do
        if i <= n then
          local w, t = measure(cell.contents)
          -- A spanning cell is charged to no single column: it is as wide as
          -- the columns under it together, and charging it to the first one
          -- would inflate that column by the width of all of them.
          if (cell.col_span or 1) == 1 then
            if w > nat[i] then nat[i] = w end
            if t > tok[i] then tok[i] = t end
          end
        end
        i = i + (cell.col_span or 1)
      end
    end
  end
  if tbl.head then scan(tbl.head.rows) end
  for _, body in ipairs(tbl.bodies or {}) do
    scan(body.head)
    scan(body.body)
  end
  if tbl.foot then scan(tbl.foot.rows) end
  return nat, tok
end

-- Longest word first, prose second. Every column is given the width its longest
-- unbreakable word needs, because a word wider than its column is printed into
-- the next one. What remains is shared among the columns that were cut, in
-- proportion to how much each lost, which sends it to the prose columns without
-- naming them. Returns nil when the minimum alone does not fit.
local function allocate(nat, tok, avail)
  local n, w, floor_sum = #nat, {}, 0
  for i = 1, n do
    w[i] = math.min(nat[i], math.max(tok[i] + 1, MIN_COL))
    floor_sum = floor_sum + w[i]
  end
  if floor_sum > avail then return nil end
  local short, total = {}, 0
  for i = 1, n do
    short[i] = nat[i] - w[i]
    total = total + short[i]
  end
  if total > 0 then
    local slack = avail - floor_sum
    for i = 1, n do
      w[i] = w[i] + math.min(short[i], slack * short[i] / total)
    end
  end
  return w
end

-- Sets the column widths and reports which mode the table has to be set in.
local function fit(tbl)
  local n = #tbl.colspecs
  if n == 0 then return 'plain' end
  local nat, tok = col_stats(tbl)
  for _, mode in ipairs(ORDER) do
    local avail = CHARS[mode] - PAD[mode] * n
    if avail > 0 then
      local w = allocate(nat, tok, avail)
      if w then
        for i = 1, n do tbl.colspecs[i][2] = w[i] / avail end
        return mode
      end
    end
  end
  -- Wider than a rotated page even at its narrowest. Nothing can prevent this
  -- from overflowing, so the overflow is put where it does least harm: each
  -- column keeps its longest word and the whole is scaled to the page.
  local total = 0
  for i = 1, n do total = total + math.max(tok[i] + 1, MIN_COL) end
  for i = 1, n do
    tbl.colspecs[i][2] = math.max(tok[i] + 1, MIN_COL) / total
  end
  return 'land'
end

-- KEEPING THE TITLE WITH WHAT IT TITLES. Every table and figure in this manual
-- is written as a heading, one line saying what it holds, then the thing itself.
-- LaTeX will break a page between any two of those three. The heading and the
-- line are taken back off the output, the space they and the first rows need is
-- reserved with \cyneed, and they are put back immediately before the table,
-- so a page break happens above the heading or not at all.
local function grab_title(out)
  local title = {}
  local last = out[#out]
  if last and last.t == 'Para' and
     #pandoc.utils.stringify(last) <= TITLE_CHARS then
    table.insert(title, 1, table.remove(out))
  end
  last = out[#out]
  -- Level 1 is a chapter, which starts its own page anyway.
  if last and last.t == 'Header' and last.level >= 2 then
    table.insert(title, 1, table.remove(out))
  end
  return title
end

-- What the title above a table or figure will occupy, in lines of body text.
-- A heading is charged four lines beyond its own text: it is set larger than the
-- body, and the space a \section leaves above and below itself is another two.
-- Undercounting here is not a rounding error, it is the heading on the page
-- before, so the estimate is deliberately generous.
local function reserve(title, base)
  local lines = base
  for _, b in ipairs(title) do
    local n = #pandoc.utils.stringify(b)
    lines = lines + math.max(1, math.ceil(n / CHARS.plain))
    lines = lines + (b.t == 'Header' and 4 or 1)
  end
  return math.min(lines, MAX_RESERVE)
end

local function image_src(block)
  local src
  pandoc.walk_block(block, {
    Image = function(im) if not src then src = im.src end end })
  return src
end

local function is_figure(block)
  if block.t == 'Figure' then return image_src(block) ~= nil end
  if block.t ~= 'Para' then return false end
  for _, inl in ipairs(block.content) do
    if inl.t ~= 'Image' and inl.t ~= 'Space' and inl.t ~= 'SoftBreak' then
      return false
    end
  end
  return image_src(block) ~= nil
end

local function process(blocks)
  local out = {}
  for _, block in ipairs(blocks) do
    -- Recurse first, so a table inside a div is fitted too.
    if block.t == 'Div' or block.t == 'BlockQuote' then
      block.content = process(block.content)
    end

    if block.t == 'Table' then
      local mode  = fit(block)
      local title = grab_title(out)
      if mode == 'land' then
        -- \clearpage on both sides, or LaTeX floats the rotated page to
        -- wherever it next fits, which is routinely several pages after the
        -- paragraph introducing the table.
        out[#out + 1] = raw('\\clearpage\\begin{landscape}\\begin{footnotesize}')
        for _, b in ipairs(title) do out[#out + 1] = b end
        out[#out + 1] = block
        out[#out + 1] = raw('\\end{footnotesize}\\end{landscape}\\clearpage')
      else
        -- 4 lines beyond the title: the rule, the header row and the first row
        -- of the table, so a table cannot start with nothing but its own header
        -- at the foot of a page.
        out[#out + 1] = raw(string.format('\\cyneed{%d\\baselineskip}',
                                          reserve(title, 4)))
        for _, b in ipairs(title) do out[#out + 1] = b end
        if SIZE[mode] then out[#out + 1] = raw('\\begingroup' .. SIZE[mode]) end
        out[#out + 1] = block
        if SIZE[mode] then out[#out + 1] = raw('\\endgroup') end
      end
    elseif is_figure(block) then
      -- The figure's own height is not known here, only in LaTeX, so the
      -- reservation is made there: \cyfigneed measures the image and asks for
      -- that much plus the lines the title needs.
      -- 5 lines beyond the title: two for the caption under the image, two for
      -- the space a float leaves above and below itself, one for rounding. The
      -- caption is the part that is easy to forget, and forgetting it puts the
      -- image on one page and its own caption on the next.
      local title = grab_title(out)
      out[#out + 1] = raw(string.format('\\cyfigneed{%s}{%d}',
                                        image_src(block), reserve(title, 5)))
      for _, b in ipairs(title) do out[#out + 1] = b end
      out[#out + 1] = block
    else
      out[#out + 1] = block
    end
  end
  return out
end

function Pandoc(doc)
  -- Only for PDF. Guarding keeps the filter safe to reuse against the site.
  if not FORMAT:match('latex') then return nil end
  doc.blocks = process(doc.blocks)
  return doc
end
