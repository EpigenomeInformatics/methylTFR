# methylTFRdeviations

Class for storing results from
[`run_methyltfr`](https://epigenomeinformatics.github.io/methylTFR/reference/run_methyltfr.md)
function.

## Details

This class inherits from
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
and most methods for that class should work for objects of this class as
well. Additionally, two accessor functions are defined for extracting
bias corrected deviations
([`deviations`](https://epigenomeinformatics.github.io/methylTFR/reference/deviations-methylTFRdeviations-method.md))
and deviation Z-scores
([`deviationZScores`](https://epigenomeinformatics.github.io/methylTFR/reference/deviationZScores-methylTFRdeviations-method.md))
