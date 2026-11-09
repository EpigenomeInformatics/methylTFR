# create_sink

Create a temporary HDF5 sink for storing methylTFR results.

## Usage

``` r
create_sink(
  files_list,
  motifs,
  tempDir = file.path(tempdir(), "methylTFR"),
  pattern = "methylTFR",
  fileext = ".h5",
  verbose = TRUE
)
```

## Arguments

- files_list:

  A character vector of sample names.

- motifs:

  A character vector of motifs.

- tempDir:

  Directory the HDF5 file is written to. Created if it does not exist.

- pattern:

  Prefix of the temporary file name.

- fileext:

  File extension of the temporary file.

- verbose:

  A logical indicating whether to print messages.

## Value

An `HDF5RealizationSink`. Its file is in `sink@filepath`.
