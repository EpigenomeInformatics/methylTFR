# resolve_annotation_file

Work out the path of the sample annotation file.

## Usage

``` r
resolve_annotation_file(annfile, sample_ann, sample_dir)
```

## Arguments

- annfile:

  Explicit path to the annotation file, or NULL.

- sample_ann:

  Name of the annotation file inside `sample_dir`.

- sample_dir:

  Directory holding the methylation call files.

## Value

The path of the annotation file as a character string.
