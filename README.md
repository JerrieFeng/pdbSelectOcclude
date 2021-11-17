
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

`pdbSelectOcclude` contains 3 functions. The visPDB function is able to
visualize the protein structure given a PDB file. The interactPDB
function allows user to interact with the 3D structure of the protein.
The occludePDB function analyzes the protein sequences and enables users
to select occluding components.

``` r
browseVignettes("pdbSelectOcclude")
```

## Contributions

The author of the package is Jerrie Feng. The visPDB function is able to
visualize the protein structure given a PDB file.

## References

Skjaerven, L., Yao X.Q., Grant B.J. (2006). Getting started with Bio3D.
Grant Lab: ComputationalBiophysics & Bioinformatics.
<http://thegrantlab.org/bio3d/articles/online/intro_vignette/Bio3D_introduction.html#references-1>.

Su, W., Johnston, B. (2021). r3dmol: Create Interactive 3D
Visualizations of Molecular Data.
Github.https://github.com/swsoyee/r3dmol.

Xiao, N. (2019). protr: R package for generating various numerical
representation schemes of pro-tein sequences. ProtrWeb.
<https://cran.r-project.org/web/packages/protr/vignettes/protr.html#1_introduction>.

## Acknowledgements

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinfor-matics, University of Toronto, Toronto, CANADA.
