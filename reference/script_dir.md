# Directory this package is installed in

Recorded in the run manifest so a results folder can be traced back to
the code that produced it. As a script this resolved the `--file=`
argument; in a package the installed location is both simpler and more
truthful.

## Usage

``` r
script_dir()
```
