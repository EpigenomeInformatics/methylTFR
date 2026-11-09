# methylTFR: transcription factor activity from DNA methylation

Scores the methylation footprint at transcription factor binding sites
and returns one bias-corrected score per motif and sample.

## Getting started

1.  Import methylation calls with
    [`readMethylome()`](https://epigenomeinformatics.github.io/methylTFR/reference/readMethylome.md).

2.  Load annotations from methylTFRAnnotationHg38 (`getTFbindsites()`,
    `getGCfreq()`, `getGenomeGC()`).

3.  Compute scores with
    [`runMethylTFR()`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFR.md)
    or
    [`runMethylTFRFromRnBeads()`](https://epigenomeinformatics.github.io/methylTFR/reference/runMethylTFRFromRnBeads.md).

4.  Explore with
    [`deviations()`](https://epigenomeinformatics.github.io/methylTFR/reference/deviations.md),
    [`deviationZScores()`](https://epigenomeinformatics.github.io/methylTFR/reference/deviationZScores.md),
    [`computeZScoreVariability()`](https://epigenomeinformatics.github.io/methylTFR/reference/computeZScoreVariability.md),
    [`differentialDeviationTest()`](https://epigenomeinformatics.github.io/methylTFR/reference/differentialDeviationTest.md),
    [`plotMotifFootprint()`](https://epigenomeinformatics.github.io/methylTFR/reference/plotMotifFootprint.md)
    and
    [`plotExpectedFootprint()`](https://epigenomeinformatics.github.io/methylTFR/reference/plotExpectedFootprint.md).

See
[`vignette("methylTFR")`](https://epigenomeinformatics.github.io/methylTFR/articles/methylTFR.md)
for a worked example.

## See also

Useful links:

- <https://github.com/EpigenomeInformatics/methylTFR>

- <https://epigenomeinformatics.github.io/methylTFR/>

- Report bugs at <https://support.bioconductor.org/t/methylTFR>

## Author

**Maintainer**: Irem B. Gündüz <irembgunduz@gmail.com>
([ORCID](https://orcid.org/0000-0003-2641-0916))

Authors:

- Sarath Kumar Murugan <sarath.murugan@outlook.com>
  ([ORCID](https://orcid.org/0000-0001-7378-7869))

- Fabian Mueller <fabian.mueller@uni-saarland.de>
  ([ORCID](https://orcid.org/0000-0001-5809-2321))

Other contributors:

- ERA-NET Transcan-Neu III - EPILUNAR (Grant 01KT2409) \[funder\]

- Saarland University (NanoBioMed Young Investigator Grant) \[funder\]
