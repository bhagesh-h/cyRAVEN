# run_cyraven

Runs the pipeline. A run that stops on an error still writes
`report.html`, carrying the diagnosis and every output produced before
the failure, and then re-raises the original error unchanged so the exit
status and the message a caller sees are what they would have been.

## Usage

``` r
run_cyraven(opt)
```

## Arguments

- opt:

  Named list of options; anything omitted takes its command-line
  default. See build_option_list() for the full set.
