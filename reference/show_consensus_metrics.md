# Illustrates importance consensus metrics

Plots a heatmap to illustrate the behaviour of different importance
consensus metrics.

## Usage

``` r
show_consensus_metrics(
  metrics = c("min", "harmonic", "geometric", "product", "average", "l2", "max")
)
```

## Arguments

- metrics:

  Character vector of metrics to show. Should be valid values for the
  `metric` argument of
  [`consensus_importance_metric()`](https://plant-food-research-open.github.io/moiraine/reference/consensus_importance_metric.md),
  i.e. "min", "max", "average", "product", "l2", "geometric",
  "harmonic".

## Value

A ggplot.
