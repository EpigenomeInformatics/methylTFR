# create_sink

Create a temp sink for storing methylTFR results

## Usage

``` r
create_sink(
  files_list,
  motifs,
  temp_dir = "methylTFR_tmp",
  pattern = "methylTFR",
  fileext = ".h5",
  verbose = TRUE
)
```

## Arguments

- files_list:

  A character vector of file names

- motifs:

  A character vector of motifs

- temp_dir:

  A character vector specifying the temp directory

- pattern:

  A character vector specifying the pattern for the temp file

- fileext:

  A character vector specifying the file extension for the temp file

- verbose:

  A logical indicating whether to print messages

## Value

A methylTFR sink
