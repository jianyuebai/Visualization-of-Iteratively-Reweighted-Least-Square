
# InteracDiagnosis

<!-- badges: start -->
<!-- badges: end -->

A Shiny tool for visualising how adding or omitting covariates affects parameter estimates, focusing on changes in the coefficient of interest, along with theoretical and bootstrap standard errors, partial (added-variable) plots, and model-comparison statistics (AIC, BIC, and likelihood-ratio tests).



## Installation

You can install the development version of InteracDiagnosis from [GitHub](https://github.com/) with:

``` r
# install.packages("InteracDiagnosis")
pak::pak("jianyuebai/InteracDiagnosis")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(InteracDiagnosis)
run_shiny_app()
[Compare models (uploaded CSV): coefficient shifts + partial plots + bootstrap SE.pdf](https://github.com/user-attachments/files/24180940/Compare.models.uploaded.CSV.coefficient.shifts.%2B.partial.plots.%2B.bootstrap.SE.pdf)


```

