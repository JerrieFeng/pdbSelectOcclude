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
#' #showInfoPDB(pdb, "1bm8")
#'
#' #showInfoPDB(pdb2, "1si4")
#'
#'
#' @references
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D.
#' bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.
#'
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). chain.pdb: Find Possible PDB Chain Breaks.
#' bio3D. https://rdrr.io/cran/bio3d/man/chain.pdb.html.
#'
#' @export getPTM
#' @export getFunctionSites
#' @export downloadPDB
#' @export showInfoPDB
#'
#' @import bio3d


getPTM <- function(accession){
  if (! requireNamespace("UniprotR", quietly = TRUE)) {
    install.packages("UniprotR", repos = "http://cran.us.r-project.org")
    install.packages("Matrix", repos = "http://cran.us.r-project.org")
  }

  #pull data
  ptm <- UniprotR::GetPTM_Processing(accession)
  gly3 <- c()
  mod3 <- c()

  #pull if not empty
  if(any(!is.na(ptm$Glycosylation))){
    gly1 <- as.list(strsplit(ptm$Glycosylation, 'CARBOHYD')[[1]])
    gly2 <- gsub(";.*", "", gly1)

    for(each in gly2){
      if(grepl(".", each)){
        first <- gsub("\\..*", "", each)#first
        last <- gsub(".*\\.", "", each) #later
        together <- seq(first, last)
        gly3 <- c(gly3, together)
      }else{
        gly3 <- c(gly3, each)
      }
    }
    gly3 <- gly3[gly3 != ""]
    gly3 <- unique(gly3)

  }
  if(any(!is.na(ptm$Modified.residue))){
    mod1 <- as.list(strsplit(ptm$Modified.residue, 'MOD_RES')[[1]])
    mod2 <- gsub(";.*", "", mod1)

    for(each in mod2){
      if(grepl(".", each)){
        first <- gsub("\\..*", "", each)#first
        last <- gsub(".*\\.", "", each) #later
        together <- seq(first, last)
        mod3 <- c(mod3, together)
      }else{
        mod3 <- c(mod3, each)
      }
    }
    mod3 <- mod3[mod3 != ""]
    mod3 <- unique(mod3)
  }

  sites <- list(Glycosylation=gly3, Modified=mod3)

  return(sites)
}


getFunctionSites <- function(accession){
  if (! requireNamespace("UniprotR", quietly = TRUE)) {
    install.packages("UniprotR", repos = "http://cran.us.r-project.org")
    install.packages("Matrix", repos = "http://cran.us.r-project.org")
    install.packages("lattice", repos = "http://cran.us.r-project.org")
    install.packages("nlme", repos = "http://cran.us.r-project.org")
  }

  #pull data
  func <- UniprotR::GetProteinFunction(accession)
  active3 <- c()
  bind3 <- c()
  dna3 <- c()
  met3 <- c()
  site3 <- c()

  #pull if not empty
  if(any(!is.na(func$Active.site))){
    active1 <- as.list(strsplit(func$Active.site, 'ACT_SITE')[[1]])
    active2 <- gsub(";.*", "", active1)

    for(each in active2){
      if(grepl(".", each)){
        first <- gsub("\\..*", "", each)#first
        last <- gsub(".*\\.", "", each) #later
        together <- seq(first, last)
        active3 <- c(active3, together)
      }else{
        active3 <- c(active3, each)
      }
    }
    active3 <- active3[active3 != ""]
    active3 <- unique(active3)
  }
  if(any(!is.na(func$Binding.site))){
    bind1 <- as.list(strsplit(func$Binding.site, 'BINDING')[[1]])
    bind2 <- gsub(";.*", "", bind1)

    for(each in bind2){
      if(grepl(".", each)){
        first <- gsub("\\..*", "", each)#first
        last <- gsub(".*\\.", "", each) #later
        together <- seq(first, last)
        bind3 <- c(bind3, together)
      }else{
        bind3 <- c(bind3, each)
      }
    }
    bind3 <- bind3[bind3 != ""]
    bind3 <- unique(bind3)
  }
  if(any(!is.na(func$DNA.binding))){
    dna1 <- as.list(strsplit(func$DNA.binding, 'DNA_BIND')[[1]])
    dna2 <- gsub(";.*", "", dna1)

    for(each in dna2){
      if(grepl(".", each)){
        first <- gsub("\\..*", "", each)#first
        last <- gsub(".*\\.", "", each) #later
        together <- seq(first, last)
        dna3 <- c(dna3, together)
      }else{
        dna3 <- c(dna3, each)
      }
    }
    dna3 <- dna3[dna3 != ""]
    dna3 <- unique(dna3)
  }
  if(any(!is.na(func$Metal.binding))){
    met1 <- as.list(strsplit(func$Metal.binding, 'METAL')[[1]])
    met2 <- gsub(";.*", "", met1)

    for(each in met2){
      if(grepl(".", each)){
        first <- gsub("\\..*", "", each)#first
        last <- gsub(".*\\.", "", each) #later
        together <- seq(first, last)
        met3 <- c(met3, together)
      }else{
        met3 <- c(met3, each)
      }
    }
    met3 <- met3[met3 != ""]
    met3 <- unique(met3)
  }
  if(any(!is.na(func$Site))){
    site1 <- as.list(strsplit(func$Site, 'SITE')[[1]])
    site2 <- gsub(";.*", "", site1)

    for(each in site2){
      if(grepl(".", each)){
        first <- gsub("\\..*", "", each)#first
        last <- gsub(".*\\.", "", each) #later
        together <- seq(first, last)
        site3 <- c(site3, together)
      }else{
        site3 <- c(site3, each)
      }
    }
    site3 <- site3[site3 != ""]
    site3 <- unique(site3)
  }

  sites <- list(Active=active3, Binding=bind3, DNA=dna3, Metal=met3, Site=site3)

  return(sites)
}


downloadPDB <- function(name, path = "."){
  pdbFile <- bio3d::get.pdb(name, path)
}




showInfoPDB <- function(pdbFile, name){
  if (! requireNamespace("ptm", quietly = TRUE)) {
    BiocManager::install("muscle")
    install.packages("MASS", repos = "http://cran.us.r-project.org")
    install.packages("ptm", repos = "http://cran.us.r-project.org")
  }

  # Check if the name and given PDB file refer to the same protein
  if(any(grepl(name, pdbFile, ignore.case=TRUE))
     & nchar(name) == 4){
    cat("Ok the PDB file and given name match :) \n")
  }else{
    stop("The PDB given and the name provided don't match! Make sure you've entered
       the correct 4-letter PDB code for the PDB file. Or check if the PDB file is named correctly.
         The name should be written within quotation marks.")
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

  #Get Uniprot info
  UPNum <- ptm::pdb.uniprot(name) #convert pdb to uniprot number
  infoFunc <- getFunctionSites(UPNum)
  infoPTM <- getPTM(UPNum)

  # Collect all info together and return it
  infoAA <- list(code_sequence = seq, code_unique = aaList)
  infoChain <- list(residue_num = chainSection, chains = chainName)

  infoPDB <- list(name = name, AA_info = infoAA, chain_info = infoChain,
                  uniprot_func = infoFunc, uniprot_ptm = infoPTM)

  return(infoPDB)

}

#[END]
