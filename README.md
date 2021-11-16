
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pdbSelectOcclude

<!-- badges: start -->
<!-- badges: end -->

## Description

`pdbSelectOcclude` is an R package for visualizing relationships between
protein structural and non-spatial data, and for selecting occluding
structures in a 3D molecular visualization environment. The objective of
the package is to help improve current selecting (aka brushing) protein
structural residues techniques. Currently, most brushing techniques
(such as point-and-click) for 3D spatial data do not work well for
occluded structures (i.e. molecular structures that are insidea protein
and not just on its surfaces). The package is targeted for those who are
interested in molecular visualization.

This includes the main components: DESCRIPTION, NAMESPACE, man
subdirectory and R subdirectory. Additionally, licence, README and
subdirectories vignettes, tests, data and inst are also explored. To
develop this package, R (version 4.0.2) and Mac platform was used.

## Installation

You can install the development version of pdbSelectOcclude like so:

``` r
require("devtools")
devtools::install_github("JerrieFeng/pdbSelectOcclude", build_vignettes = TRUE)
library("pdbSelectOcclude")
```

To run the shinyApp: Under construction

## Overview

``` r
ls("package:pdbSelectOcclude")
data(package = "pdbSelectOcclude") # optional
```

`pdbSelectOcclude` contains 3 functions. The InfCriteriaCalculation
function calculates the information criteria values. Specifically,
Bayesian information criterion (BIC), Akaike information criterion (AIC)
and Integrated Complete Likelihood (ICL) are calculated. The
InfCriteriaPlot generates a plot of information criteria values.
NormFactors is a function that calculates normalization factors via
Trimmed Mean of M-values (TMM). The runTestingPackage is the function
that launches the shiny app for this package. The package also contains
RNA sequencing dataset GeneCounts. Refer to package vignettes for more
details.

``` r
browseVignettes("pdbSelectOcclude")
```

## Contributions

The author of the package is Jerrie Feng. The InfCriteriaCalculation
function makes use of map function from mclust R package to generate ICL
values. The stats R package is used for generating multinomially
distributed random number vectors. (If applies: Part of the code for
InfCriteriaCalculation function has been taken from <NamePackage> R
package. Section of the borrowed code is clearly indicated and
referenced in the InfCriteriaCalculation help file). The InfCriteriaPlot
makes use of the graphics R package. NormFactors uses Trimmed Mean of
M-values (TMM) as implemented in edgeR R package.

## References

## Acknowledgements

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinfor-matics, University of Toronto, Toronto, CANADA.
