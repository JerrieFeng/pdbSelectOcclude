#' visualize a PDB file
#'
#' Function that visualizes the given PDB file
#'
#' @param pdbfile A PDB file, could be downloaded from PDB online
#' @return nothing for now. It just shows the 3D visuals
#'
#' @references Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D. bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.
#'
#' @export
#' @import bio3d
visPDB <- function(pdbfile) {

  pdb <- read.pdb(pdbfile)
  modes <- bio3d::nma(pdb)
  cij <- bio3d::dccm(modes)
  net <- bio3d::cna(cij, cutoff.cij=0.3)
  vmd(net, pdb, launch = TRUE)

  return()
}

# [END]
