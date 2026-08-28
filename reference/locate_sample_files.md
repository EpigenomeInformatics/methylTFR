# locate_sample_files

Build and check the list of per-sample methylation files.

## Usage

``` r
locate_sample_files(samples, sample_dir, sampleColName, full_path)
```

## Arguments

- samples:

  The sample annotation as a `data.frame`.

- sample_dir:

  Directory holding the methylation call files.

- sampleColName:

  Column name holding the file names.

- full_path:

  if TRUE, the annotation file holds full paths.

## Value

A character vector of existing file paths.
