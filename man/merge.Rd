% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runShiny.R
\name{merge}
\alias{merge}
\title{Compile all functions into Shiny}
\usage{
merge(
  pdbFile,
  name,
  styleSel,
  codes = NULL,
  polarity = NULL,
  letter = "ALL",
  values = NULL
)
}
\arguments{
\item{pdbFile}{A PDB file, could be downloaded from PDB online}

\item{name}{the 4-letter PDB codes/identifiers of the PDB file
to be visualized.}

\item{index}{Valid index that refers to a chain found in the protein. Run 'showInfoPDB'
to get all possible chains in chain_info}
}
\value{
Returns output-text stating the chain, and it displays the 3D model
of the protein with the selected chain visible.
}
\description{
Function that takes a PDB file and visualizes it as a 3D model.
}
\examples{
# merge(pdb2, "1SI4", "m_style_stick", "VAL", NULL, "ALL", "Metal")
# runModel(pdb2, "1SI4")


}
\references{
Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D.
bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.

Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). chain.pdb: Find Possible PDB Chain Breaks.
bio3D. https://rdrr.io/cran/bio3d/man/chain.pdb.html.

Su, W. (2021). Introduction to r3dmol. r3dmol.
https://cran.r-project.org/web/packages/r3dmol/vignettes/r3dmol.html.
}
