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

## Bundled data

The package code is MIT licensed (see `LICENSE`). The files in
`inst/extdata` have the following sources and terms.

- `example_data.rda`, `msites_sub.rda`, `gcdist_subset.rda`,
  `gcdist_BATF.rda`: CpG methylation calls, and the GC windows that
  overlap them, from ENCODE whole-genome bisulfite sequencing experiment
  ENCSR663MXB (GRCh38,
  <https://www.encodeproject.org/experiments/ENCSR663MXB/>). ENCODE data
  are released without restrictions on use
  (<https://www.encodeproject.org/help/citing-encode/>). Please
  acknowledge the ENCODE Consortium and the production laboratory, and
  cite ENCODE Project Consortium (2012) Nature 489:57-74
  (doi:10.1038/nature11247) and Luo et al. (2020) Nucleic Acids Res
  48:D882-D889 (doi:10.1093/nar/gkz1062).

- `BATF_tf_bindsites.rda`, `BATF_gcfreqs.rda`: BATF binding sites on
  chr1 and the genome-wide GC frequencies of the motif, derived from the
  JASPAR2020 BATF motif through methylTFRAnnotationHg38. JASPAR data are
  released under CC BY 4.0 (<https://jaspar.elixir.no/about/>). Cite
  Fornes et al. (2020) Nucleic Acids Res 48:D87-D92
  (doi:10.1093/nar/gkz1001).

- `immuneDeviations.rda`, `tc_mem.rda`, `tc_naive.rda`: simulated
  deviation scores. The source pseudobulk methylomes (Gunduz et al.
  2025, bioRxiv, doi:10.1101/2025.09.09.675101; Wang et al. 2025, Nat
  Genet, doi:10.1038/s41588-025-02479-6) are under controlled access in
  dbGaP (phs003204.v1.p1). Only per-cell-type motif means and SDs were
  taken from them; every sample was drawn from a normal distribution
  with those parameters. No individual-level values, methylation calls
  or donor identifiers are included. Released under the package's MIT
  license. Made by `inst/scripts/tcell_data.R` and
  `inst/scripts/vignette_data.R`.

- `*.tsv.gz`: small simulated files in each supported input format, made
  by `inst/scripts/read_methylome_data.R`.

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
