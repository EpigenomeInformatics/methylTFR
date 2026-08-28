# check_file_run_inputs

Validate the file-specific arguments of `run_methyltfr`.

## Usage

``` r
check_file_run_inputs(filetype, sampleColName, full_path)
```

## Arguments

- filetype:

  File type of the methylation call files.

- sampleColName:

  Column name holding the file names.

- full_path:

  if TRUE, the annotation file holds full paths.

## Value

Invisible `NULL`. Called for the errors it raises.
