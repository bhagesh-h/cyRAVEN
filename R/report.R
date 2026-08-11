# SECTION 13 -- RUN REPORT
# =============================================================================
#
# WHY THIS FILE EXISTS. A run writes several dozen tables and figures into a
# directory. The reading ORDER is the product: the diagnostics article sets out
# which checks can invalidate which results, and the statistics article sets out
# that a difference is read as adjusted p, then effect size, then against the
# gate's own uncertainty, then against the counting uncertainty. A directory
# listing enforces none of that. It presents a failed staining QC and a headline
# p-value as two files of equal standing, sorted alphabetically.
#
# The field says the same thing about adoption. The blockers reported for
# automated analysis are accessibility and user confidence rather than accuracy
# (Popp et al. 2025, Cytometry A 107:189).
#
# WHY IT IS WRITTEN BY HAND RATHER THAN THROUGH rmarkdown. Rendering Markdown to
# HTML needs pandoc, which is a system binary rather than an R package, and the
# runtime image deliberately does not carry one. A report that only works outside
# the container would be a report the documented execution path cannot produce.
# The HTML here is assembled directly, so it needs nothing that is not already
# present.
#
# WHY EVERYTHING IS EMBEDDED. The report is one file that carries every figure
# and every table inside it, as data URIs and as JSON. It references nothing.
# That is what makes it the thing you can attach to an email, put in a
# supplement, or archive on its own and still have the whole result: a report
# whose images live beside it becomes a page of broken icons the moment it is
# moved, and a result that cannot survive being moved is not a record.
#
# The cost is size, and it is a real one. Figures are embedded once, at full
# resolution, and displayed scaled down by CSS; the same bytes serve the
# on-screen figure and the full-resolution download, so nothing is stored twice.
# The final size is logged, because a 40 MB HTML file should not be a surprise.
#
# WHY THE TABLES ARE JSON RATHER THAN MARKUP. Every table is searchable, can be
# paged at 10/50/100/all rows, and can be exported to CSV. Doing that over
# pre-rendered <tr> elements means the export re-parses the DOM and the search
# hides rather than filters. Carrying the data and rendering from it makes all
# three operations read the same array, so what you export is what you filtered.

#' Size above which a table is named rather than embedded
#'
#' Set with `options(cyRAVEN.report_table_max_mb = )`. The default of 8 MB sits
#' above every summary table a run writes and below the per-cell exports, whose
#' row count is the number of cells rather than the number of samples.
#' @keywords internal
report_table_max_bytes <- function() {
  mb <- getOption("cyRAVEN.report_table_max_mb", 8)
  if (!is.numeric(mb) || length(mb) != 1L || is.na(mb) || mb <= 0) mb <- 8
  mb * 1024^2
}

#' Escape text for inclusion in HTML
#' @param x character vector
#' @keywords internal
html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

#' Escape a string for embedding inside a <script> block as JSON
#'
#' `</script>` anywhere inside the data would end the block early, so the slash
#' of any closing tag is escaped. JSON treats `\/` as `/`, so the data is
#' unchanged by it.
#' @param x character vector.
#' @keywords internal
json_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\r", " ", x, fixed = TRUE)
  x <- gsub("\n", " ", x, fixed = TRUE)
  x <- gsub("\t", " ", x, fixed = TRUE)
  # Control characters JSON forbids unescaped.
  x <- gsub("[\001-\037]", " ", x)
  gsub("</", "<\\/", x, fixed = TRUE)
}

#' Render a data.frame as a JSON array-of-arrays with a header
#' @param d a data.frame
#' @keywords internal
json_table <- function(d) {
  cols <- paste0("[", paste0("\"", json_escape(names(d)), "\"", collapse = ","), "]")
  if (!nrow(d)) return(paste0("{\"cols\":", cols, ",\"rows\":[]}"))
  cells <- lapply(d, function(col) paste0("\"", json_escape(col), "\""))
  rows <- do.call(paste, c(cells, sep = ","))
  paste0("{\"cols\":", cols, ",\"rows\":[",
         paste0("[", rows, "]", collapse = ","), "]}")
}

#' Legacy HTML table renderer
#'
#' Retained because the interactive tables are rendered from JSON in the
#' browser; this is the static fallback used when a table is too small for
#' searching to be worth the controls.
#' @param d a data.frame
#' @param max_rows rows beyond which the table is truncated, with a note
#' @keywords internal
html_table <- function(d, max_rows = 20L) {
  if (is.null(d) || !nrow(d)) return("<p class='none'>No rows.</p>")
  note <- ""
  if (nrow(d) > max_rows) {
    note <- sprintf("<p class='none'>Showing %d of %d rows; the file has the rest.</p>",
                    max_rows, nrow(d))
    d <- utils::head(d, max_rows)
  }
  hdr <- paste0("<th>", html_escape(names(d)), "</th>", collapse = "")
  body <- apply(d, 1, function(r)
    paste0("<tr>", paste0("<td>", html_escape(r), "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead><tr>", hdr, "</tr></thead><tbody>",
         paste(body, collapse = ""), "</tbody></table>", note)
}

#' Build one report section, embedding its figures and tables
#'
#' @param outdir the results directory.
#' @param id anchor id, unique per section.
#' @param title section heading, stated as what the section reports.
#' @param description what the figures and tables below show and how to read them.
#' @param figures figure filenames, embedded when present.
#' @param tables table filenames, embedded when present.
#' @param body optional raw HTML inserted before the figures.
#' @param open whether the section starts expanded.
#' @return list(html, nav, n_fig, n_tab, bytes, used)
#' @keywords internal
report_section <- function(outdir, id, title, description, figures = character(0),
                           tables = character(0), body = NULL, open = FALSE) {
  fig_present <- figures[file.exists(file.path(outdir, figures))]
  tab_present <- tables[file.exists(file.path(outdir, tables))]

  # A per-cell export carries one row per CELL, so on a real cohort it is
  # hundreds of thousands of rows. It is read by software rather than by a
  # person scrolling a report, and embedding it would multiply the file size for
  # something nobody reads there. Such a table is NAMED with its size instead of
  # being dropped, so its absence is a stated fact rather than a silent gap.
  oversize <- character(0)
  if (length(tab_present)) {
    big <- vapply(tab_present, function(t)
      file.size(file.path(outdir, t)) > report_table_max_bytes(), logical(1))
    oversize <- tab_present[big]
    tab_present <- tab_present[!big]
  }

  # A section with nothing to show is omitted rather than rendered empty: an
  # empty heading reads as "this check found nothing", which is a different
  # claim from "this check did not run".
  if (!length(fig_present) && !length(tab_present) && !length(oversize) &&
      is.null(body))
    return(list(html = "", nav = "", n_fig = 0L, n_tab = 0L, bytes = 0,
                used = character(0)))

  h <- c(sprintf("<section class='sec' id='%s'>", id),
         sprintf("<details%s><summary><span class='chev' aria-hidden='true'></span>",
                 if (open) " open" else ""),
         sprintf("<span class='sec-title'>%s</span>", html_escape(title)),
         sprintf("<span class='sec-count'>%s</span></summary>",
                 html_escape(paste0(length(fig_present), " fig / ",
                                    length(tab_present), " tab"))),
         "<div class='sec-body'>",
         sprintf("<p class='q'>%s</p>", description))
  if (!is.null(body)) h <- c(h, body)

  nav <- character(0)
  bytes <- 0
  for (f in fig_present) {
    p <- file.path(outdir, f)
    uri <- file_data_uri(p, "image/png")
    if (is.na(uri)) next
    bytes <- bytes + file.size(p)
    fid <- paste0("fig-", gsub("[^A-Za-z0-9]+", "-", sub("[.]png$", "", f)))
    h <- c(h, sprintf(paste0(
      "<figure class='fig' id='%s'>",
      "<div class='figbox'><img src='%s' alt='%s' loading='lazy' ",
      "onclick='cyZoom(this)'/>",
      "<button class='zoom' title='Zoom' onclick='cyZoom(this.parentNode.querySelector(\"img\"))'>&#9974;</button>",
      "</div>",
      "<figcaption><span class='fn'>%s</span>",
      "<a class='dl' download='%s' href='%s'>Full resolution PNG</a>",
      "</figcaption></figure>"),
      fid, uri, html_escape(f), html_escape(f), html_escape(f), uri))
    nav <- c(nav, sprintf("<a class='nav-fig' href='#%s'>%s</a>", fid,
                          html_escape(sub("[.]png$", "", f))))
  }

  for (t in tab_present) {
    p <- file.path(outdir, t)
    d <- tryCatch(utils::read.csv(p, stringsAsFactors = FALSE,
                                  check.names = FALSE),
                  error = function(e) NULL)
    if (is.null(d)) next
    bytes <- bytes + file.size(p)
    tid <- paste0("tab-", gsub("[^A-Za-z0-9]+", "-", sub("[.]csv$", "", t)))
    h <- c(h, sprintf(paste0(
      "<div class='tab' id='%s'>",
      "<div class='tabhead'><h3>%s</h3>",
      "<span class='dim'>%d rows &times; %d cols</span>",
      "<input class='search' type='search' placeholder='Search this table' ",
      "oninput='cyFilter(\"%s\")' aria-label='Search %s'/>",
      "<select class='pagesel' onchange='cyFilter(\"%s\")' aria-label='Rows to show'>",
      "<option value='10'>10 rows</option><option value='50' selected>50 rows</option>",
      "<option value='100'>100 rows</option><option value='0'>All rows</option></select>",
      "<button class='dl' onclick='cyCsv(\"%s\")'>Export CSV</button>",
      "</div><div class='tabwrap'><table></table></div>",
      "<p class='none shown'></p></div>"),
      tid, html_escape(t), nrow(d), ncol(d), tid, html_escape(t), tid, tid))
    h <- c(h, sprintf("<script type='application/json' id='%s-data'>%s</script>",
                      tid, json_table(d)))
    nav <- c(nav, sprintf("<a class='nav-tab' href='#%s'>%s</a>", tid,
                          html_escape(t)))
  }

  for (t in oversize) {
    p <- file.path(outdir, t)
    nr <- tryCatch(length(readLines(p, warn = FALSE)) - 1L,
                   error = function(e) NA_integer_)
    h <- c(h, sprintf(paste0(
      "<p class='none'><b>%s</b> is in the results directory but not embedded ",
      "here: %s row(s), %s MB. It is a per-cell export, one row per cell, read ",
      "by software rather than read in a report; carrying it would multiply the ",
      "size of this file. Every summary derived from it is embedded above.</p>"),
      html_escape(t), format(nr, big.mark = ","),
      format(round(file.size(p) / 1024^2, 1), nsmall = 1)))
  }

  h <- c(h, "</div></details></section>")
  list(html = paste(h, collapse = "\n"),
       nav = paste0(sprintf("<div class='nav-sec'><a class='nav-h' href='#%s'>%s</a>",
                            id, html_escape(title)),
                    paste(nav, collapse = ""), "</div>"),
       n_fig = length(fig_present), n_tab = length(tab_present), bytes = bytes,
       used = c(fig_present, tab_present))
}

#' Stylesheet for the run report
#' @keywords internal
report_css <- function() {
  paste0(
  ":root{--fg:#1a1a1a;--mut:#5b6470;--line:#e2e5ea;--bg:#fff;--panel:#f7f8fa;",
  "--accent:#0a7d4a;--warn:#8a6100;--stop:#a4231c;",
  "--ui:ui-sans-serif,system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;",
  "--mono:ui-monospace,SFMono-Regular,'SF Mono',Menlo,Consolas,monospace}",
  # One font stack for every element, set once on the root and inherited. Tables,
  # buttons, inputs and selects do NOT inherit font by default in any browser,
  # so they are named explicitly rather than left to the user-agent stylesheet.
  "*{box-sizing:border-box}",
  "html{font-family:var(--ui);font-size:15px;color:var(--fg);background:var(--bg)}",
  "body{margin:0;line-height:1.55;font-family:var(--ui)}",
  "button,input,select,table,th,td{font-family:var(--ui);font-size:inherit;color:inherit}",
  "code,kbd,.mono,.fn,.dim{font-family:var(--mono)}",
  # Layout: fixed sidebar, scrolling main.
  ".wrap{display:flex;align-items:flex-start;min-height:100vh}",
  "nav.side{position:sticky;top:0;flex:0 0 268px;height:100vh;overflow-y:auto;",
  "border-right:1px solid var(--line);background:var(--panel);padding:1rem .75rem}",
  "nav.side h2{font-size:.78rem;text-transform:uppercase;letter-spacing:.08em;",
  "color:var(--mut);margin:.2rem 0 .6rem .4rem}",
  ".nav-sec{margin-bottom:.35rem}",
  "nav.side a{display:block;text-decoration:none;color:var(--fg);border-radius:5px;",
  "padding:.22rem .45rem;font-size:.84rem;overflow:hidden;text-overflow:ellipsis;",
  "white-space:nowrap}",
  "nav.side a:hover{background:#e8ebf0}",
  "a.nav-h{font-weight:600;font-size:.87rem}",
  "a.nav-fig,a.nav-tab{padding-left:1.35rem;color:var(--mut);font-family:var(--mono);",
  "font-size:.76rem}",
  "a.nav-fig:before{content:'\\25A6\\00a0';color:var(--accent)}",
  "a.nav-tab:before{content:'\\2263\\00a0';color:var(--mut)}",
  "main{flex:1 1 auto;min-width:0;max-width:1180px;padding:1.5rem 2rem 5rem}",
  "h1{font-size:1.6rem;margin:0 0 .15rem}",
  # Collapsible sections.
  "section.sec{border:1px solid var(--line);border-radius:8px;margin:.7rem 0;",
  "background:var(--bg);overflow:hidden}",
  "summary{cursor:pointer;padding:.7rem .9rem;display:flex;align-items:center;",
  "gap:.6rem;background:var(--panel);user-select:none;list-style:none}",
  "summary::-webkit-details-marker{display:none}",
  "summary:hover{background:#eef1f5}",
  ".chev{width:0;height:0;border-left:6px solid var(--mut);",
  "border-top:4.5px solid transparent;border-bottom:4.5px solid transparent;",
  "transition:transform .15s;flex:0 0 auto}",
  "details[open] .chev{transform:rotate(90deg)}",
  ".sec-title{font-weight:600;flex:1 1 auto}",
  ".sec-count{color:var(--mut);font-size:.76rem;font-family:var(--mono)}",
  ".sec-body{padding:.4rem 1rem 1.1rem}",
  "p.q{color:var(--mut);margin:.4rem 0 1rem}",
  "p.none{color:var(--mut);font-size:.82rem}",
  # Figures: every figure occupies the same box whatever its native aspect
  # ratio, so the page does not jump between a wide strip and a tall grid.
  "figure.fig{margin:1rem 0 1.4rem}",
  ".figbox{position:relative;height:420px;border:1px solid var(--line);",
  "border-radius:6px;background:var(--panel);display:flex;align-items:center;",
  "justify-content:center;overflow:hidden}",
  ".figbox img{max-width:100%;max-height:100%;width:auto;height:auto;",
  "object-fit:contain;cursor:zoom-in;display:block}",
  "button.zoom{position:absolute;top:.4rem;right:.4rem;border:1px solid var(--line);",
  "background:rgba(255,255,255,.92);border-radius:5px;cursor:pointer;",
  "padding:.15rem .4rem;font-size:.95rem;line-height:1}",
  "button.zoom:hover{background:#fff}",
  "figcaption{display:flex;justify-content:space-between;align-items:center;",
  "gap:1rem;margin-top:.4rem;font-size:.78rem;color:var(--mut)}",
  "a.dl,button.dl{font-size:.76rem;color:var(--accent);text-decoration:none;",
  "border:1px solid var(--accent);border-radius:5px;padding:.16rem .5rem;",
  "background:none;cursor:pointer;white-space:nowrap}",
  "a.dl:hover,button.dl:hover{background:var(--accent);color:#fff}",
  # Lightbox.
  "#lb{position:fixed;inset:0;background:rgba(12,14,18,.93);display:none;",
  "z-index:99;overflow:auto;cursor:zoom-out}",
  "#lb.on{display:block}",
  "#lb img{display:block;margin:0 auto;max-width:none}",
  "#lbbar{position:fixed;top:0;left:0;right:0;padding:.5rem .9rem;display:flex;",
  "gap:.5rem;align-items:center;background:rgba(12,14,18,.86);color:#eef1f5;",
  "font-size:.8rem;z-index:100}",
  "#lbbar button{background:none;color:#eef1f5;border:1px solid #556;",
  "border-radius:5px;padding:.15rem .55rem;cursor:pointer}",
  "#lbbar button:hover{background:#2a2f38}",
  "#lbname{font-family:var(--mono);flex:1 1 auto;overflow:hidden;",
  "text-overflow:ellipsis;white-space:nowrap}",
  # Tables.
  ".tab{margin:1.2rem 0}",
  ".tabhead{display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;",
  "margin-bottom:.35rem}",
  ".tabhead h3{margin:0;font-size:.86rem;font-family:var(--mono);font-weight:600;",
  "flex:0 1 auto}",
  ".dim{color:var(--mut);font-size:.74rem;flex:1 1 auto}",
  "input.search,select.pagesel{border:1px solid var(--line);border-radius:5px;",
  "padding:.2rem .45rem;font-size:.78rem;background:#fff}",
  "input.search{width:12rem}",
  ".tabwrap{max-height:26rem;overflow:auto;border:1px solid var(--line);",
  "border-radius:6px}",
  ".tabwrap table{border-collapse:collapse;width:100%;font-size:.78rem}",
  ".tabwrap th,.tabwrap td{border-bottom:1px solid var(--line);padding:.26rem .5rem;",
  "text-align:left;white-space:nowrap}",
  ".tabwrap th{background:var(--panel);position:sticky;top:0;font-weight:600;",
  "cursor:pointer}",
  ".tabwrap th:hover{background:#e8ebf0}",
  ".tabwrap tr:nth-child(even) td{background:#fbfcfd}",
  # Banners.
  # Diagnosis block on a failed run.
  "h3{font-size:.9rem;margin:1.1rem 0 .3rem;text-transform:uppercase;",
  "letter-spacing:.05em;color:var(--mut)}",
  ".sec-body h3:first-of-type{margin-top:.2rem}",
  "pre{font-family:var(--mono);font-size:.79rem;line-height:1.45;overflow-x:auto;",
  "border:1px solid var(--line);border-radius:6px;padding:.6rem .75rem;",
  "background:var(--panel);margin:.2rem 0;white-space:pre-wrap;word-break:break-word}",
  "pre.err{background:#fdecea;border-color:#f5b5ae;color:#7d1a15;white-space:pre-wrap}",
  "pre.log{max-height:22rem;overflow-y:auto;white-space:pre}",
  "pre.cmd{background:#14181d;border-color:#14181d;color:#e6edf3}",
  ".banner{padding:.7rem .9rem;border-radius:6px;margin:.8rem 0;font-size:.9rem}",
  ".ok{background:#eefaf0;border:1px solid #b6e3c2}",
  ".warn{background:#fff6e5;border:1px solid #f0d9a8}",
  ".stop{background:#fdecea;border:1px solid #f5b5ae}",
  "@media print{nav.side{display:none}details{open:true}}",
  # The report is read on laptops and projected in meetings; below 900px the
  # sidebar becomes a normal block above the content rather than disappearing.
  "@media(max-width:900px){.wrap{display:block}nav.side{position:static;height:auto;",
  "width:auto;border-right:none;border-bottom:1px solid var(--line)}",
  "main{padding:1rem}.figbox{height:300px}}")
}

#' Behaviour for the run report
#' @keywords internal
report_js <- function() {
  paste0(
  "var CY={};\n",
  # Tables are parsed once on load from their JSON blocks and kept in memory,
  # so search, paging and export all read the same array.
  "function cyInit(){document.querySelectorAll(\"script[type='application/json']\")",
  ".forEach(function(s){CY[s.id.replace(/-data$/,'')]=JSON.parse(s.textContent);});",
  "Object.keys(CY).forEach(function(id){CY[id].sort=-1;cyFilter(id);});}\n",
  "function cyRows(id){var d=CY[id];var box=document.getElementById(id);",
  "var q=box.querySelector('input.search').value.trim().toLowerCase();",
  "var rows=d.rows;",
  "if(q){var terms=q.split(/\\s+/);rows=rows.filter(function(r){",
  "var hay=r.join(' ').toLowerCase();",
  "return terms.every(function(t){return hay.indexOf(t)>-1;});});}",
  # Sorting is numeric when every value in the column parses as a number, and
  # lexical otherwise. Deciding per column rather than per cell stops a column
  # of p-values sorting as text because one cell reads "NA".
  "if(d.sort>=0){var c=d.sort,dir=d.dir||1;",
  "var num=rows.every(function(r){return r[c]===''||!isNaN(parseFloat(r[c]));});",
  "rows=rows.slice().sort(function(a,b){var x=a[c],y=b[c];",
  "if(num){x=parseFloat(x);y=parseFloat(y);",
  "if(isNaN(x))return 1;if(isNaN(y))return -1;return (x-y)*dir;}",
  "return x.localeCompare(y)*dir;});}",
  "return rows;}\n",
  "function cyFilter(id){var d=CY[id];var box=document.getElementById(id);",
  "var n=parseInt(box.querySelector('select.pagesel').value,10);",
  "var rows=cyRows(id);var total=rows.length;",
  "var show=(n===0)?rows:rows.slice(0,n);",
  "var h='<thead><tr>'+d.cols.map(function(c,i){",
  "var mark=(d.sort===i)?(d.dir===1?' \\u25B2':' \\u25BC'):'';",
  "return \"<th onclick='cySort(\\\"\"+id+\"\\\",\"+i+\")'>\"+cyEsc(c)+mark+'</th>';",
  "}).join('')+'</tr></thead><tbody>';",
  "h+=show.map(function(r){return '<tr>'+r.map(function(v){",
  "return '<td>'+cyEsc(v)+'</td>';}).join('')+'</tr>';}).join('');",
  "box.querySelector('table').innerHTML=h+'</tbody>';",
  "box.querySelector('p.shown').textContent='Showing '+show.length+' of '+total+",
  "(total===d.rows.length?'':' matching')+' row'+(total===1?'':'s')+",
  "'. Export CSV writes what is shown after search and sorting.';}\n",
  "function cySort(id,i){var d=CY[id];",
  "if(d.sort===i){d.dir=(d.dir===1)?-1:1;}else{d.sort=i;d.dir=1;}cyFilter(id);}\n",
  "function cyEsc(s){return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;')",
  ".replace(/>/g,'&gt;');}\n",
  # Export writes exactly what is on screen after search and sorting, which is
  # why it reads cyRows() rather than the untouched array.
  "function cyCsv(id){var d=CY[id];var box=document.getElementById(id);",
  "var n=parseInt(box.querySelector('select.pagesel').value,10);",
  "var rows=cyRows(id);if(n>0)rows=rows.slice(0,n);",
  "var q=function(v){v=String(v);return /[\",\\n]/.test(v)?'\"'+v.replace(/\"/g,'\"\"')+'\"':v;};",
  "var csv=[d.cols.map(q).join(',')].concat(rows.map(function(r){",
  "return r.map(q).join(',');})).join('\\n');",
  "var b=new Blob([csv],{type:'text/csv;charset=utf-8'});",
  "var a=document.createElement('a');a.href=URL.createObjectURL(b);",
  "a.download=box.querySelector('h3').textContent;document.body.appendChild(a);",
  "a.click();document.body.removeChild(a);URL.revokeObjectURL(a.href);}\n",
  # Zoom opens the embedded full-resolution image; there is no second copy and
  # no request, so it works from a file:// path and with no network.
  "var LBZ=1;\n",
  "function cyZoom(img){var lb=document.getElementById('lb');",
  "var i=document.getElementById('lbimg');i.src=img.src;LBZ=1;cyZset(1);",
  "document.getElementById('lbname').textContent=img.getAttribute('alt')||'';",
  "document.getElementById('lbdl').href=img.src;",
  "document.getElementById('lbdl').download=img.getAttribute('alt')||'figure.png';",
  "lb.classList.add('on');document.body.style.overflow='hidden';}\n",
  "function cyZset(z){LBZ=Math.min(8,Math.max(0.1,z));var i=document.getElementById('lbimg');",
  "i.style.width=(LBZ*100)+'%';document.getElementById('lbpct').textContent=",
  "Math.round(LBZ*100)+'%';}\n",
  "function cyClose(){document.getElementById('lb').classList.remove('on');",
  "document.body.style.overflow='';}\n",
  "document.addEventListener('keydown',function(e){",
  "if(e.key==='Escape')cyClose();",
  "if(document.getElementById('lb').classList.contains('on')){",
  "if(e.key==='+'||e.key==='=')cyZset(LBZ*1.25);",
  "if(e.key==='-')cyZset(LBZ/1.25);}});\n",
  "function cyAll(open){document.querySelectorAll('details').forEach(function(d){",
  "d.open=open;});}\n",
  "document.addEventListener('DOMContentLoaded',cyInit);\n")
}

#' Write a single self-contained HTML report of everything a run produced
#'
#' Presents the outputs in the order the documentation says they must be read,
#' because each stage can invalidate the ones after it. Sections whose files are
#' absent are omitted rather than shown empty.
#'
#' The file references nothing: every figure is embedded at full resolution as a
#' data URI and every table as JSON, so the report can be moved, emailed or
#' archived on its own. Figures can be zoomed and downloaded at full resolution;
#' tables can be searched, sorted, paged at 10/50/100/all rows and exported to
#' CSV.
#'
#' When `failure` is supplied the same report is written for a run that stopped
#' early, with a diagnosis section first and every output produced up to the
#' failure embedded below it. The partial output is usually where the evidence
#' is, so it is kept rather than discarded.
#'
#' @param outdir the results directory, which is also where the report is written
#' @param opt parsed options, used for the header
#' @param verdicts per-sample staining verdicts, used for the summary banner
#' @param failure a condition or message; when given, the report is written as
#'   the record of a failed run rather than a completed one.
#' @return the path, invisibly, or NULL when there is nothing to report
#' @export
write_run_report <- function(outdir, opt = NULL, verdicts = NULL,
                             failure = NULL) {
  if (!dir.exists(outdir)) return(invisible(NULL))
  path <- file.path(outdir, "report.html")
  failed <- !is.null(failure)

  # The banner. A run whose staining QC excluded samples says so before anything
  # else, because every frequency below is computed without them.
  banner <- ""
  if (failed) banner <- failure_block(failure, outdir, opt)
  else if (!is.null(verdicts) && length(verdicts)) {
    st <- vapply(verdicts, function(v) v$qc_status %||% "pass", character(1))
    ctl <- vapply(verdicts, function(v) isTRUE(v$is_control), logical(1))
    nfail <- sum(st == "failed" & !ctl)
    banner <- if (nfail)
      sprintf("<div class='banner stop'><b>%d of %d sample(s) failed staining QC and are excluded from every test below.</b> See staining_qc.csv for the reason each was excluded.</div>",
              nfail, length(st))
    else
      sprintf("<div class='banner ok'>All %d declared sample(s) passed staining QC.</div>",
              sum(!ctl))
  }

  secs <- list(
    report_section(outdir, "s1", "1. Gate placement",
      paste("Per sample, the scatter gate boundary and every marker threshold,",
            "each drawn on the density it was derived from. A cut sitting in the",
            "trough between two modes is determined by the data. A cut on the",
            "flank of a single mode is a quantile fallback: not invalid, but",
            "carrying no evidence that the two populations separate. That",
            "distinction is visible here and in no downstream table, which is why",
            "this section comes first."),
      figures = c("recon_diagnostics.png", "gating_qc.png"), open = TRUE),

    report_section(outdir, "s2", "2. Acquisition stability",
      paste("The Time channel binned into equal-width intervals, tracking the",
            "event rate and each channel's median. A sustained trough is a partial",
            "clog, a spike is usually a bubble, a step is a settings change; any",
            "of them makes one file two instruments over its run, so a single",
            "threshold suits neither half. acquisition_qc.csv gives one verdict",
            "per file. acquisition_qc_impact.csv states how far each population",
            "would move if the flagged intervals were dropped, and that is the",
            "number the decision rests on: compare it against the same",
            "population's gate uncertainty in section 6. Nothing is removed",
            "unless --drop-unstable-events was given."),
      figures = "acquisition_qc.png",
      tables = c("acquisition_qc.csv", "acquisition_qc_impact.csv",
                 "acquisition_qc_bins.csv")),

    report_section(outdir, "s3", "3. Staining quality control",
      paste("One verdict per sample. A sample with no resolvable CD45+ mode has",
            "no usable parent gate, so its percentages are fractions of an",
            "arbitrary scatter region rather than of leukocytes. Such samples are",
            "excluded from every test below and the exclusion is recorded here,",
            "unless --include-qc-failed was given, in which case the verdict",
            "column still records it while qc_status reads pass."),
      tables = "staining_qc.csv"),

    report_section(outdir, "s4", "4. Phenotype concordance",
      paste("Measured marker intensity across the declared populations, z-scored",
            "over the run. Identity here is declared before the data are",
            "examined, so this figure serves the inverse function of its",
            "equivalent in clustering-first analysis: a population that does not",
            "express the markers its own definition requires falsifies the gate",
            "that produced it. The cohort heatmap shows each group's share of",
            "every population after normalising the groups to a common cell",
            "count, so an uneven split is not an artefact of unequal group",
            "sizes."),
      figures = c("population_marker_heatmap.png", "cohort_composition_heatmap.png")),

    report_section(outdir, "s5", "5. Threshold provenance and spillover spreading",
      paste("thresholds_used.csv records how every cut was obtained: a density",
            "minimum, a quantile fallback, an unstained or fluorescence-minus-one",
            "control, or a manual override. spreading_receivers.csv pairs each",
            "marker's fallback rate with how much wider its negative population",
            "becomes when another channel is bright. Reading the two columns",
            "together separates a cut that failed for an optical reason, which no",
            "gating strategy recovers and which a different fluorochrome",
            "assignment would fix, from one that failed for want of positive",
            "events or a titration problem."),
      tables = c("thresholds_used.csv", "threshold_scale_qc.csv",
                 "fmo_agreement.csv", "spreading_receivers.csv",
                 "spreading_pairs.csv")),

    report_section(outdir, "s6", "6. Gate uncertainty and detection limits",
      paste("Each frequency carries two separate uncertainties. The gate",
            "uncertainty is propagated from the thresholds behind the population",
            "and says how far the number moves when the cut moves;",
            "uncertainty_budget.csv attributes it to the individual markers, so",
            "the gate worth fixing is named. The counting uncertainty comes from",
            "the number of events observed. They are different guarantees: a cut",
            "through a wide empty gap is well determined however few events lie",
            "beyond it. detection_limits.png classifies every population against",
            "the limits of detection and quantification set by its own parent",
            "gate, and a population below them cannot be recovered by re-gating."),
      figures = c("frequency_uncertainty.png", "uncertainty_budget.png",
                  "detection_limits.png"),
      tables = c("uncertainty_budget.csv", "threshold_uncertainty.csv")),

    report_section(outdir, "s7", "7. Population abundance and the shared embedding",
      paste("pct_of_cd45_pos is the reportable quantity. count is an event count,",
            "set by how long the tube was run, and is not a cell number. The UMAP",
            "is computed once across all samples, so positions are comparable",
            "between them; per-sample embeddings would not be. The unsupervised",
            "clustering is computed without reference to the specification, which",
            "is what lets it contradict it: a cluster dominated by no declared",
            "label is a population the specification does not describe, and a",
            "declared label spread thinly across many clusters covers several",
            "distinct phenotypes."),
      figures = c("population_frequencies.png", "umap_overview.png",
                  "umap_overview_by_group.png", "umap_markers.png",
                  "umap_density.png", "umap_density_by_group.png",
                  "umap_multigraph_overlay.png", "unsupervised_clusters.png"),
      tables = c("population_frequencies.csv", "population_marker_mfi.csv",
                 "gate_counts.csv", "cluster_gate_agreement_populations.csv",
                 "cluster_gate_agreement_clusters.csv",
                 "subcluster_marker_shifts.csv", "unsupervised_clusters.csv",
                 "cells_umap.csv")),

    report_section(outdir, "s8", "8. Between-group differences",
      paste("Tests are on per-sample values with donors as the replicates, not on",
            "pooled cells. Read each row in four steps: the adjusted p-value,",
            "then Cliff's delta as the effect size, then difference_over_gate_u,",
            "which expresses the difference in units of the gate's own",
            "uncertainty, then difference_over_total_u, which adds the counting",
            "uncertainty. A difference many times the gate uncertainty is not an",
            "artefact of threshold placement, but it is still not a finding until",
            "it survives correction across populations. Compositional results",
            "test the same abundances as centred log-ratios, because percentages",
            "of one parent cannot all move independently."),
      figures = c("group_comparison.png", "marker_state.png",
                  "functional_markers.png", "population_ratios.png",
                  "absolute_counts.png", "absolute_counts_qc.png"),
      tables = c("group_comparison_stats.csv", "compositional_concordance.csv",
                 "compositional_clr_stats.csv", "marker_state_stats.csv",
                 "functional_markers.csv", "functional_markers_stats.csv",
                 "population_ratios.csv", "population_ratios_stats.csv",
                 "absolute_counts.csv", "absolute_counts_stats.csv",
                 "absolute_counts_raw.csv")),

    report_section(outdir, "s9", "9. Confounding and batch structure",
      paste("A variable confounds only when it both differs between the groups",
            "and associates with the outcome, and the two conditions are reported",
            "separately rather than combined into one verdict. threshold_drift",
            "asks whether the cuts themselves track the study group, which would",
            "make a difference definitional. batch_group_confounding.csv reports",
            "Cramer's V between batch and group: correction is refused above the",
            "configured threshold, because at that level of association removing",
            "the batch effect and removing the finding are the same operation."),
      figures = c("threshold_drift.png", "batch_diagnostic.png"),
      tables = c("threshold_drift_stats.csv", "confounding_diagnostics.csv",
                 "batch_group_confounding.csv", "batch_mixing_stats.csv",
                 "marker_batch_drift.csv")),

    report_section(outdir, "s10", "10. Agreement with an accepted baseline",
      paste("This run's per-sample values against a previously accepted run,",
            "written when --baseline was given. The within-run peer check finds",
            "one deviant tube against its peers. It cannot find a cohort that",
            "moved as a whole, because the peer median moves with it; only a",
            "baseline from a different run can."),
      tables = c("specification_conformance.csv", "specification_changes.csv")))

  # ---- completeness sweep ---------------------------------------------------
  # Every section above names its files, which is what puts them in reading
  # order. A file this run wrote but no section names would be invisible here
  # while appearing to be covered, so whatever is left over is collected rather
  # than dropped. If this section is ever non-empty for a standard run, a named
  # section is missing an entry.
  used <- unlist(lapply(secs, `[[`, "used"))
  rest_fig <- setdiff(list.files(outdir, "[.]png$"), used)
  rest_tab <- setdiff(list.files(outdir, "[.]csv$"), used)
  secs <- c(secs, list(report_section(outdir, "s11", "11. Further outputs",
    paste("Everything else this run wrote, in no particular order. These files",
          "are produced by optional flags or by stages the sections above do not",
          "cover; each is documented in the Output article."),
    figures = rest_fig, tables = rest_tab)))

  # Provenance describes the other sections, so it is only worth writing when
  # there are some. A directory holding nothing produces no report rather than
  # a page whose sole content is a note about how to read the content.
  if (length(Filter(function(s) nzchar(s$html), secs)) || failed)
  secs <- c(secs, list(report_section(outdir, "s12", "12. Provenance",
    paste("What produced this folder, and what this file is."),
    tables = "patient_metadata_english.csv",
    body = paste0(
      "<p>This report is self-contained. Every figure above is embedded at full",
      " resolution and every table is embedded in full, so the file can be",
      " moved, attached to an email or archived on its own with nothing lost;",
      " it references no other file and needs no network. Click any figure to",
      " zoom, or use its download link for the original PNG. Any table can be",
      " searched, sorted by clicking a column heading, shown 10, 50, 100 or all",
      " rows at a time, and exported to CSV exactly as filtered and sorted.</p>",
      "<p>The complete run record stays in the results directory as ",
      "<span class='mono'>run_manifest.txt</span>",
      if (file.exists(file.path(outdir, "miflowcyt.md")))
        " and <span class='mono'>miflowcyt.md</span>" else "",
      ", which record the package versions, the invocation and every option in",
      " force.</p>"))))

  secs <- Filter(function(s) nzchar(s$html), secs)
  # A failed run that produced nothing still gets a report: the diagnosis is
  # the point of it, and an empty directory is exactly the case where the user
  # has least else to go on.
  if (!length(secs) && !failed) return(invisible(NULL))

  nfig <- sum(vapply(secs, `[[`, integer(1), "n_fig"))
  ntab <- sum(vapply(secs, `[[`, integer(1), "n_tab"))

  head <- c(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>cyRAVEN run report</title>",
    sprintf("<style>%s</style>", report_css()), "</head><body>",
    "<div id='lb' onclick='if(event.target.id===\"lb\")cyClose()'>",
    "<div id='lbbar'><span id='lbname'></span>",
    "<button onclick='cyZset(LBZ/1.25)'>&minus;</button>",
    "<span id='lbpct'>100%</span>",
    "<button onclick='cyZset(LBZ*1.25)'>+</button>",
    "<button onclick='cyZset(1)'>Fit</button>",
    "<a class='dl' id='lbdl' download href='#'>Download</a>",
    "<button onclick='cyClose()'>Close</button></div>",
    "<img id='lbimg' alt='' onclick='event.stopPropagation()'/></div>",
    "<div class='wrap'>",
    "<nav class='side'><h2>Contents</h2>",
    paste(vapply(secs, `[[`, character(1), "nav"), collapse = "\n"),
    "</nav><main>",
    sprintf("<h1>cyRAVEN run report%s</h1>", if (failed) " &mdash; FAILED" else ""),
    sprintf("<p class='q'>%s &middot; cyRAVEN %s &middot; %d figures, %d tables, all embedded</p>",
            html_escape(format(Sys.time(), tz = "UTC", usetz = TRUE)),
            html_escape(tryCatch(as.character(utils::packageVersion("cyRAVEN")),
                                 error = function(e) "unknown")),
            nfig, ntab),
    banner,
    sprintf(paste0("<div class='banner warn'>%s",
            "<span style='float:right'>",
            "<button class='dl' onclick='cyAll(true)'>Expand all</button> ",
            "<button class='dl' onclick='cyAll(false)'>Collapse all</button>",
            "</span></div>"),
      if (failed)
        paste("What follows is everything the run wrote before it stopped, in the",
              "order a completed run is read in. Stages after the failure are",
              "absent, so a section missing here did not run rather than finding",
              "nothing.")
      else
        paste("Read the sections in the order given. Each one can invalidate the",
              "sections after it, so a result taken from the bottom without the",
              "top is not supported by this run.")))

  writeLines(c(head, vapply(secs, `[[`, character(1), "html"),
               "</main></div>", sprintf("<script>%s</script>", report_js()),
               "</body></html>"), path)
  sz <- file.size(path)
  log_msg("wrote report.html (", if (failed) "FAILED run, " else "",
          length(secs), " section(s), ", nfig,
          " figure(s) and ", ntab, " table(s) embedded, ",
          format(round(sz / 1024^2, 1), nsmall = 1), " MB, ",
          "self-contained: it references no other file)")
  invisible(path)
}
