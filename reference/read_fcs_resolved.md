# Read one FCS file and resolve marker symbols

WHY marker symbols, not channel names: detector/channel names are
instrument-specific and `alter.names = TRUE` mangles them (e.g.
`[[RB613]]-A` becomes `X..RB613...A`). The biological identity lives in
\$PnS, which on this instrument reads "CD45 : SparkUV-387 - Area" – the
symbol is the text before " : ". Building an embedding on channel names
(as the original template did) means no marker identity ever reaches the
analysis.

## Usage

``` r
read_fcs_resolved(path, sample_id = NULL, max_events = 0L)
```

## Arguments

- path:

  File path.

- sample_id:

  The sample id.

- max_events:

  The max events. Default `0L`.

## Value

list(exprs, marker_cols, all_cols, keywords, n_events, file, sample_id)

## Details

Note: FCS files larger than 99,999,999 bytes cannot express their data
offset in the 8-character HEADER field and write 0 there; the true
offsets live in the \$BEGINDATA/\$ENDDATA keywords. flowCore handles
this fallback internally.
