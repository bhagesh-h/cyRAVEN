-- Put every table on its own landscape page.
--
-- WHY A FILTER AND NOT A LATEX OPTION. The tables in this documentation are wide
-- -- clinical_association.csv has fifteen columns, the option reference has a
-- description column that runs to a paragraph -- and on a portrait page LaTeX has
-- two ways to deal with that, both bad. It either overruns the margin, so the
-- right-hand columns are printed off the paper, or it compresses the columns
-- until adjacent cells touch and the table becomes unreadable. That is the
-- overlap this exists to prevent.
--
-- Rotating the page gives roughly 40% more usable width, which is the difference
-- between those tables fitting and not. It cannot be set globally in the preamble
-- because `landscape` is an environment that has to WRAP each table, and pandoc
-- emits the tables itself.
--
-- WHY NOT EVERY TABLE. A three-column table is perfectly readable in portrait,
-- and rotating it costs the reader a page turn and a head tilt for nothing. The
-- threshold below is the point where a table starts needing the width.
--
-- WHY \clearpage ON BOTH SIDES. Without it LaTeX floats the rotated page to
-- wherever it next fits, which is routinely several pages after the paragraph
-- that introduces the table. Clearing before and after pins it in place: the
-- table appears where it was written, at the cost of some white space above it.

-- Columns from which a table is wide enough to earn a rotated page. Four fits
-- portrait comfortably; five is where the description columns in the option
-- reference start to crowd.
local MIN_COLS = 5

-- A table with few columns can still be wide if one of them holds prose. This is
-- the character count in the widest single cell past which the table is rotated
-- regardless of how many columns it has.
local MIN_CELL_CHARS = 120

-- A Cell is a wrapper, not a list of blocks: stringify takes its .contents.
-- Passing the Cell itself raises "table expected, got pandoc Cell", which aborts
-- pandoc entirely -- and with continue-on-error on the build step, that failure
-- reported as success and the only symptom was a 404 on the download link.
local function cell_text(cell)
  return pandoc.utils.stringify(cell.contents or cell)
end

local function widest_cell(tbl)
  local widest = 0
  local function scan(rows)
    for _, row in ipairs(rows or {}) do
      for _, cell in ipairs(row.cells or {}) do
        local n = #cell_text(cell)
        if n > widest then widest = n end
      end
    end
  end
  for _, head in ipairs(tbl.head and tbl.head.rows or {}) do
    for _, cell in ipairs(head.cells or {}) do
      local n = #cell_text(cell)
      if n > widest then widest = n end
    end
  end
  for _, body in ipairs(tbl.bodies or {}) do
    scan(body.body)
  end
  return widest
end

local function n_cols(tbl)
  return #(tbl.colspecs or {})
end

function Table(tbl)
  -- Only for PDF. The same document is not rendered to HTML by this path, but
  -- guarding keeps the filter safe to reuse.
  if not FORMAT:match('latex') then return nil end

  local cols = n_cols(tbl)
  local widest = widest_cell(tbl)
  if cols < MIN_COLS and widest < MIN_CELL_CHARS then return nil end

  -- \small as well as the rotation: a fifteen-column table needs both, and on a
  -- table that only just crossed the threshold the smaller type is not noticeable
  -- against the body text of a table.
  return {
    pandoc.RawBlock('latex', '\\clearpage\\begin{landscape}\\begin{footnotesize}'),
    tbl,
    pandoc.RawBlock('latex', '\\end{footnotesize}\\end{landscape}\\clearpage')
  }
end
