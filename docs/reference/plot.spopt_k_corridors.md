# Plot k-Diverse Corridors

Renders all k corridors on a single plot. Rank 1 is drawn with full
opacity and maximum line width; higher ranks fade and thin
proportionally.

## Usage

``` r
# S3 method for class 'spopt_k_corridors'
plot(x, ...)
```

## Arguments

- x:

  An `spopt_k_corridors` object from
  [`route_k_corridors`](https://walker-data.com/spopt-r/reference/route_k_corridors.md).

- ...:

  Additional arguments passed to the initial
  [`plot()`](https://rspatial.github.io/terra/reference/plot.html) call.

## Value

No return value, called for side effects (draws the corridors to the
active graphics device).
