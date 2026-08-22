# Read a file and return it as a data URI

Read a file and return it as a data URI

## Usage

``` r
file_data_uri(path, mime = "image/png")
```

## Arguments

- path:

  File path.

- mime:

  MIME type, defaulting to PNG.

## Value

`data:<mime>;base64,...`, or NA when the file cannot be read.
