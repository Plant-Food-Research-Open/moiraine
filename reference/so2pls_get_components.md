# Get list of latent components from sO2PLS results

Extracts the list of joint and specific latent components sO2PLS
results.

## Usage

``` r
so2pls_get_components(so2pls_res)
```

## Arguments

- so2pls_res:

  The sO2PLS results generated with
  [`get_output_so2pls()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

## Value

A list with the following elements:

- `joint`: character vector with the name of the joint latent
  components.

- `specific`: named list of length 2. Each element corresponds to a
  dataset (names of the list are the datasets name), and is a character
  vector with the name of the specific latent components for the
  corresponding dataset.
