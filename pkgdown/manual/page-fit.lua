-- Lay the manual out in two columns: fit each table to a column or to the full
-- page width, decide which figures span both columns, and keep every table and
-- figure with the title above it.
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
-- TWO COLUMNS. The text is set in two columns with the multicol package. A
-- table that fits one column at a readable size, and is short enough not to
-- need a page break, stays in the column. Anything wider or longer is set
-- across the full page: the columns are closed above it, balanced, and opened
-- again below. The reason for the split is longtable. It is the only table
-- pandoc writes that can break across pages, and it cannot be used inside
-- multicols, because both take over LaTeX's output routine. So a table in a
-- column is rewritten as a plain tabular, which cannot break, and only a table
-- set full width may run on to the next page.
--
-- CHARACTERS, NOT MILLIMETRES. Widths are computed in characters of the body
-- face and converted to fractions of the line at the end. What a cell holds is
-- something the cell gives directly. Characters are not all one width, and the
-- ones that are wider than average -- monospace, capitals, digits -- are exactly
-- what the narrow columns of these tables hold, so each is weighted rather than
-- counted.

-- WHERE A TABLE CAN GO, tried in order, first fit wins. A table is shrunk before
-- it leaves its column, and rotated only when the full width has run out.
--
-- The page is A4 less 14mm margins, a 182mm line, split into two 88mm columns
-- 6mm apart. The class is 8pt with Bricolage Grotesque set at 96%, so body text
-- is about 7.7pt and \small about 6.7pt. The scale was chosen so that one
-- character averages the same width as the Inter it replaced: about a fifth of
-- Inter's point size in millimetres, measured on the earlier A5 build (68
-- characters on a 120mm line at 8.5pt). A column holds 56 characters of body
-- text. Landscape is the same page turned, a 263mm line.
-- \footnotesize is not offered: at 5.6pt it is smaller than anyone should be
-- asked to read.
--
-- `pad` is both gutters of a column in characters of that mode's face: the
-- preamble sets \tabcolsep to 4pt, so 8pt per column.
local MODES = {
  { name = 'col',        where = 'col',  chars = 56,  pad = 1.8, size = nil },
  { name = 'col-small',  where = 'col',  chars = 64,  pad = 2.1, size = '\\small' },
  { name = 'full',       where = 'full', chars = 116, pad = 1.8, size = nil },
  { name = 'full-small', where = 'full', chars = 132, pad = 2.1, size = '\\small' },
  { name = 'land',       where = 'land', chars = 160, pad = 1.8, size = nil },
}
local COL_CHARS, FULL_CHARS = 56, 116

-- A table in a column cannot break, so it must fit in one. A column is about 75
-- lines of body text; a table estimated at more than this goes full width, where
-- longtable can carry it on to the next page. Kept well under a column so the
-- heading and the line introducing it fit above the table too.
local COL_MAX_LINES = 40

-- A FIGURE SPANS BOTH COLUMNS WHEN IT IS A WIDE STRIP. At this width-to-height
-- ratio a figure set in one 88mm column is under 44mm tall, and its axis labels
-- shrink with it to below what reads. The concept illustrations (16:9) and the
-- panel figures (up to 1.7) stay in the column, where they are 50mm or more tall
-- and their lettering is still larger than the body text.
local WIDE_FIGURE = 1.9

-- No column narrower than this, whatever its content: a 3-character column of
-- one-word cells still needs somewhere to put a word that turns out to be long.
local MIN_COL = 8

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

local OPEN_COLS  = '\\begin{multicols}{2}'
local CLOSE_COLS = '\\end{multicols}'

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

-- Every row of the table, header and body alike, as a list of rows.
local function all_rows(tbl)
  local rows = {}
  local function add(list) for _, r in ipairs(list or {}) do rows[#rows + 1] = r end end
  if tbl.head then add(tbl.head.rows) end
  for _, body in ipairs(tbl.bodies or {}) do
    add(body.head)
    add(body.body)
  end
  if tbl.foot then add(tbl.foot.rows) end
  return rows
end

-- Per column: the widest cell and the longest word. Per row: each cell's width,
-- kept so the height of the table can be estimated once the widths are known.
local function col_stats(tbl)
  local n = #tbl.colspecs
  local nat, tok, widths = {}, {}, {}
  for i = 1, n do nat[i], tok[i] = 0, 0 end
  for r, row in ipairs(all_rows(tbl)) do
    widths[r] = {}
    local i = 1
    for _, cell in ipairs(row.cells or {}) do
      if i <= n then
        local w, t = measure(cell.contents)
        -- A spanning cell is charged to no single column: it is as wide as the
        -- columns under it together, and charging it to the first one would
        -- inflate that column by the width of all of them.
        if (cell.col_span or 1) == 1 then
          if w > nat[i] then nat[i] = w end
          if t > tok[i] then tok[i] = t end
          widths[r][i] = w
        end
      end
      i = i + (cell.col_span or 1)
    end
  end
  return nat, tok, widths
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

-- Lines the table will take at these column widths: each row is as tall as its
-- most wrapped cell, plus a line for the rules.
local function table_lines(widths, w)
  local lines = 1
  for _, row in ipairs(widths) do
    local tallest = 1
    for i, cw in pairs(row) do
      if w[i] and w[i] > 0 then
        tallest = math.max(tallest, math.ceil(cw / w[i]))
      end
    end
    lines = lines + tallest
  end
  return lines
end

-- Sets the column widths and returns the mode the table has to be set in. In a
-- nested block (a block quote) the columns cannot be closed around a table, so
-- only the column modes are tried, and the narrowest is forced if none fits.
local function fit(tbl, nested)
  local n = #tbl.colspecs
  if n == 0 then return MODES[1] end
  local nat, tok, widths = col_stats(tbl)
  for _, mode in ipairs(MODES) do
    if not (nested and mode.where ~= 'col') then
      local avail = mode.chars - mode.pad * n
      if avail > 0 then
        local w = allocate(nat, tok, avail)
        if w and (mode.where ~= 'col' or nested or
                  table_lines(widths, w) <= COL_MAX_LINES) then
          for i = 1, n do tbl.colspecs[i][2] = w[i] / avail end
          return mode
        end
      end
    end
  end
  -- Wider than a rotated page even at its narrowest. Nothing can prevent this
  -- from overflowing, so the overflow is put where it does least harm: each
  -- column keeps its longest word and the whole is scaled to the line.
  local total = 0
  for i = 1, n do total = total + math.max(tok[i] + 1, MIN_COL) end
  for i = 1, n do
    tbl.colspecs[i][2] = math.max(tok[i] + 1, MIN_COL) / total
  end
  return nested and MODES[2] or MODES[#MODES]
end

local function replace_plain(s, old, new)
  local a, b = s:find(old, 1, true)
  if not a then return s end
  return s:sub(1, a - 1) .. new .. s:sub(b + 1)
end

-- LONGTABLE TO TABULAR, for a table set inside a column. pandoc writes every
-- table as a longtable: the header once for the first page (\endfirsthead),
-- again for continuation pages (\endhead), and the closing rule as the last
-- foot (\endlastfoot), all before the first body row. A tabular has none of
-- that, so the repeated header is dropped and the closing rule is moved to the
-- end, where a tabular expects it. Done on pandoc's own output rather than by
-- writing the table from scratch, so cell contents, code and links are set
-- exactly as they are everywhere else.
local function as_tabular(tbl)
  local s = pandoc.write(pandoc.Pandoc({ tbl }), 'latex')
  s = replace_plain(s, '\\begin{longtable}[]', '\\begin{tabular}')
  local fh = s:find('\\endfirsthead', 1, true)
  local eh, eh_end = s:find('\\endhead', 1, true)
  if fh and eh then
    s = s:sub(1, fh - 1) .. s:sub(eh_end + 1)
  elseif eh then
    s = s:sub(1, eh - 1) .. s:sub(eh_end + 1)
  end
  s = replace_plain(s, '\\bottomrule\\noalign{}\n\\endlastfoot', '')
  s = replace_plain(s, '\\end{longtable}', '\\bottomrule\\noalign{}\n\\end{tabular}')
  -- \nopagebreak first: the space above the table would otherwise be a legal
  -- break between the table and the line introducing it.
  return '\\par\\nopagebreak\\smallskip\\noindent\n' .. s .. '\n\\par\\smallskip'
end

-- KEEPING THE TITLE WITH WHAT IT TITLES. Every table and figure in this manual
-- is written as a heading, one line saying what it holds, then the thing itself.
-- LaTeX will break a page between any two of those three. The heading and the
-- line are taken back off the output and put back immediately before the table,
-- with what keeps them together around all three, so a break happens above the
-- heading or not at all.
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

-- What the title above a full-width table or figure will occupy, in lines of
-- body text. A heading is charged four lines beyond its own text: it is set
-- larger than the body, and the space a \section leaves above and below itself
-- is another two. Undercounting here is not a rounding error, it is the heading
-- on the page before, so the estimate is deliberately generous.
local function reserve(title, base)
  local lines = base
  for _, b in ipairs(title) do
    local n = #pandoc.utils.stringify(b)
    lines = lines + math.max(1, math.ceil(n / FULL_CHARS))
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

-- Width over height, read from the PNG header: bytes 17 to 24 are the width and
-- height as big-endian integers. Anything that is not a readable PNG returns
-- nil and is set in the column, which is never wrong, only sometimes small.
local function png_aspect(src)
  local f = io.open(src, 'rb')
  if not f then return nil end
  local head = f:read(24)
  f:close()
  if not head or #head < 24 or head:sub(2, 4) ~= 'PNG' then return nil end
  local function u32(i)
    local a, b, c, d = head:byte(i, i + 3)
    return ((a * 256 + b) * 256 + c) * 256 + d
  end
  local w, h = u32(17), u32(21)
  if h == 0 then return nil end
  return w / h
end

local function process(blocks, nested)
  local out = {}
  -- multicols is opened lazily, just before the first block that goes in it, so
  -- a full-width table straight after a chapter heading does not leave an empty
  -- column environment behind it.
  local want_cols, cols_open = not nested, false
  local function open_cols()
    if want_cols and not cols_open then
      out[#out + 1] = raw(OPEN_COLS)
      cols_open = true
    end
  end
  local function close_cols()
    if not cols_open then return end
    -- The title just taken off the output may have been all the columns held.
    if out[#out] and out[#out].t == 'RawBlock' and out[#out].text == OPEN_COLS then
      table.remove(out)
    else
      out[#out + 1] = raw(CLOSE_COLS)
    end
    cols_open = false
  end
  local function emit(b) out[#out + 1] = b end

  for _, block in ipairs(blocks) do
    -- Recurse first, so a table inside a block quote is fitted too.
    if block.t == 'Div' or block.t == 'BlockQuote' then
      block.content = process(block.content, true)
    end

    if not nested and block.t == 'Header' and block.level == 1 then
      -- A chapter title spans both columns and starts the columns afresh.
      close_cols()
      emit(block)
    elseif block.t == 'Table' then
      local mode  = fit(block, nested)
      local title = grab_title(out)
      if mode.where == 'col' then
        open_cols()
        emit(raw('\\cykeepbegin'))
        for _, b in ipairs(title) do emit(b) end
        if mode.size then emit(raw('\\begingroup' .. mode.size)) end
        emit(raw(as_tabular(block)))
        if mode.size then emit(raw('\\endgroup')) end
        emit(raw('\\cykeepend'))
      elseif mode.where == 'land' then
        close_cols()
        -- \clearpage on both sides, or LaTeX floats the rotated page to
        -- wherever it next fits, which is routinely several pages after the
        -- paragraph introducing the table.
        emit(raw('\\clearpage\\begin{landscape}'))
        for _, b in ipairs(title) do emit(b) end
        emit(block)
        emit(raw('\\end{landscape}\\clearpage'))
      else
        close_cols()
        -- 4 lines beyond the title: the rule, the header row and the first row
        -- of the table, so a table cannot start with nothing but its own header
        -- at the foot of a page.
        emit(raw(string.format('\\cyneed{%d\\baselineskip}', reserve(title, 4))))
        for _, b in ipairs(title) do emit(b) end
        if mode.size then emit(raw('\\begingroup' .. mode.size)) end
        emit(block)
        if mode.size then emit(raw('\\endgroup')) end
      end
    elseif is_figure(block) then
      local src   = image_src(block)
      local title = grab_title(out)
      local ratio = png_aspect(src)
      if not nested and ratio and ratio >= WIDE_FIGURE then
        close_cols()
        -- The figure's own height is not known here, only in LaTeX, so the
        -- reservation is made there: \cyfigneed measures the image and asks
        -- for that much plus the lines the title needs.
        -- 5 lines beyond the title: two for the caption under the image, two
        -- for the space a float leaves above and below itself, one for
        -- rounding. The caption is the part that is easy to forget, and
        -- forgetting it puts the image on one page and its own caption on the
        -- next.
        emit(raw(string.format('\\cyfigneed{%s}{%d}', src, reserve(title, 5))))
        for _, b in ipairs(title) do emit(b) end
        emit(block)
      else
        open_cols()
        emit(raw('\\cykeepbegin'))
        for _, b in ipairs(title) do emit(b) end
        -- The space above the figure is otherwise a legal break.
        if #title > 0 then emit(raw('\\nopagebreak')) end
        emit(block)
        emit(raw('\\cykeepend'))
      end
    else
      open_cols()
      emit(block)
    end
  end
  close_cols()
  return out
end

function Pandoc(doc)
  -- Only for PDF. Guarding keeps the filter safe to reuse against the site.
  if not FORMAT:match('latex') then return nil end
  doc.blocks = process(doc.blocks, false)
  return doc
end
