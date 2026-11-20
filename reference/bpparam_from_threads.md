# bpparam_from_threads

Build the BiocParallel back-end used to spread the motifs of one chunk
over workers.

## Usage

``` r
bpparam_from_threads(threads)
```

## Arguments

- threads:

  Thread count for parallel processing.

## Value

A `BiocParallelParam` object.

## Details

A forking back-end is used where the platform supports it and a socket
back-end on Windows, so that `threads` has the same meaning on every
platform. `threads = 1` runs serially in the current process.
