#' Show the information on chosen PDB file
#'
#' Function that takes a PDB file as input, and then extracts the AA (amino acid)
#' information, and chain information of the chosen protein.
#'
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' that was chosen as 'pdbFile' input
#'
#' @return Returns a list of containing 2 pieces information from the PDB file:
#' AA and chain information.
#'
#' @examples
#' # Using pdb2 and pdb file available with package
#' showInfoPDB(pdb, "1bm8")
#'
#' showInfoPDB(pdb2, "1si4")
#'
#' @references
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D.
#' bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.
#'
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). chain.pdb: Find Possible PDB Chain Breaks.
#' bio3D. https://rdrr.io/cran/bio3d/man/chain.pdb.html.
#'
#' @export
#' @import bio3d


showInfoPDB <- function(pdbFile, name){

  # Check if the name and given PDB file refer to the same protein
  if(any(grepl(name, pdbFile, ignore.case=TRUE))
     & nchar(name) == 4){
    cat("Ok the PDB file and given name match :) \n")
  }else{
    stop("The PDB given and the name provided don't match! Make sure you've entered
       the correct 4-letter PDB code for the PDB file. Or check if the PDB file is named correctly.")
  }

  # Obtain the AA residue
  seq <- bio3d::pdbseq(pdbFile, inds = NULL, aa1 = FALSE)
  # Unique AA found in protein
  aaList <- unique(seq)

  # Find potential chain breaks based on connective C-alpha/peptide bond atom separation
  chainBreaks <- bio3d::chain.pdb(pdbFile)
  chainIndex <- !duplicated(chainBreaks)
  # Residue number of where the new chain begins
  chainSection <- seq_along(chainBreaks)[chainIndex]
  # Get one-letter chain vector name
  chainName <- chainBreaks[chainIndex]

  # Collect all info together and return it
  infoAA <- list(code_sequence = seq, code_unique = aaList)
  infoChain <- list(residue_num = chainSection, chains = chainName)
  infoPDB <- list(name = name, AA_info = infoAA, chain_info = infoChain)

  return(infoPDB)

}

#[END]
