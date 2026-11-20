# deviations

Extract bias corrected deviations from methylTFRdeviations object.

## Usage

``` r
# S4 method for class 'methylTFRdeviations'
deviations(x)
```

## Arguments

- x:

  methylTFRdeviations object.

## Value

A matrix of bias corrected deviations.

## Examples

``` r
# Load the data
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
# Get deviations
deviations(tc_mem)
#>               Tc-Mem_01   Tc-Mem_02   Tc-Mem_03   Tc-Mem_04    Tc-Mem_05
#> FOXF2        0.01494122  0.01369928 -0.01658407 -0.00226295 -0.006950625
#> FOXD1       -0.05959598  0.02830943 -0.09710344 -0.02045687 -0.053317843
#> IRF2        -0.10449205 -0.13012348 -0.11231972 -0.09466112 -0.098712731
#> MZF1(var.2) -0.06747888 -0.08052227 -0.05915465 -0.08524569 -0.086930799
#> MAX::MYC    -0.19711284 -0.20767966 -0.16780734 -0.19513375 -0.231953931
#> PPARG       -0.05230212 -0.04212760 -0.05674918 -0.07438880 -0.044912732
#> PAX6         0.01515875 -0.02251858 -0.02195243 -0.03301345 -0.033578557
#> PBX1        -0.02239313 -0.07460262 -0.05639729 -0.03780550  0.008966850
#> RORA        -0.03127471 -0.13307714 -0.06685338 -0.13247760 -0.087208509
#> RORA(var.2) -0.05525696 -0.02819860 -0.06655285 -0.05332298 -0.041200472
```
