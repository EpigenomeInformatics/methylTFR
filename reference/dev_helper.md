# dev_helper

Ratio of the mean methylation in the centre bin to the mean of the two
outer flank bins.

## Usage

``` r
dev_helper(data, breaks = c(-250, -200, -25, 25, 200, 250))
```

## Arguments

- data:

  A `data.table` with the columns `x` (distance to the site centre) and
  `avg_methyl`.

- breaks:

  Six strictly increasing break points defining five bins.

## Value

A numeric value, or `NA` if the centre bin or either outer flank bin is
empty.
