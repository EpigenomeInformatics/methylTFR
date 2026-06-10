# differential_deviation_test

Differential analysis is to test which motifs are having significant
deviations among different groups. Adapted from
chromVAR::differential_test

## Usage

``` r
differential_deviation_test(
  deviations,
  groups = NULL,
  motifs = rownames(deviations),
  alternative = c("two.sided", "less", "greater"),
  parametric = TRUE,
  padjMethod = "BH"
)
```

## Arguments

- deviations:

  a `data.frame` contains deviations for each motif

- groups:

  a character vector of group names or colnames(deviations)

- motifs:

  a character vector of motif names used in analysis or
  rownames(deviations)

- alternative:

  a character string specifying the alternative hypothesis, must be one
  of "two.sided" (default), "greater" or "less".

- parametric:

  logical, if TRUE, parametric tests are used, otherwise non-parametric
  tests are used.

- padjMethod:

  method for p-value adjustment, default is "BH"

## Value

a `data.frame` contains motifs with corresponding p-value and adjusted
p-value

## Examples

``` r
# Load example data
load(system.file("extdata",
    "tc_mem.rda",
    package = "methylTFR"
))
load(system.file("extdata",
    "tc_naive.rda",
    package = "methylTFR"
))
devs <- cbind(tc_mem, tc_naive)

# Construct group labels from sample names
get_groupname <- function(x) {
    return(unlist(strsplit(x, split = "_"))[1])
}
groups <- sub(".bedGraph", "", unlist(lapply(
    FUN = get_groupname,
    X = colnames(devs)
)))

# Perform differential analysis
tc_result <- differential_deviation_test(
    deviations = devs,
    groups = groups,
    motifs = rownames(devs),
    alternative = "two.sided",
    parametric = TRUE,
    padjMethod = "BH"
)
```
