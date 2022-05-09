#' Visualize Uniprot feature keys
#'
#' Function that takes a PDB file and visualizes it as a 3D model. The user to
#' select on a scale to view uniprot feature keys
#'
#' @param accession Enter accession key.
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' to be visualized.
#'
#' @return Displays the 3D model of the protein with the selected uniprot values.
#'
#' @examples
#' # Using pdb file available with package
#' # ?idpr
#' #selUniProt("DNA", pdb, "1bm8")
#'
#' @references
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D.
#' bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.
#'
#' McFadden, W.M. (2020). Charge and Hydropathy Vignette.
#' https://rdrr.io/bioc/idpr/f/vignettes/chargeHydropathy-vignette.Rmd.
#'
#' @export checkValid
#' @export selUniProt
#'
#' @import bio3d
#' @import dplyr
#' @import r3dmol


checkValid <- function(infoPDB){
  valid_func <- c()
  valid_ptm <- c()
  i <- 1
  for(each in infoPDB$uniprot_func){
    if(length(each) != 0){
      valid_func <- append(valid_func, names(infoPDB$uniprot_func[i]))
    }
    i <- i+1
  }
  i <- 1
  for(each in infoPDB$uniprot_ptm){
    if(length(each) != 0){
      valid_ptm <- append(valid_ptm, names(infoPDB$uniprot_ptm[i]))
    }
    i <- i+1
  }

  valid <- c(valid_func, valid_ptm)
  return(valid)
}



selUniProt <- function(values, pdbFile, name){

  # Get info from showInfoPDB
  infoPDB <- showInfoPDB(pdbFile, name)

  #Categorize
  valid_func <- c()
  valid_ptm <- c()
  i <- 1
  for(each in infoPDB$uniprot_func){
    if(length(each) != 0){
      valid_func <- append(valid_func, names(infoPDB$uniprot_func[i]))
    }
    i <- i+1
  }
  i <- 1
  for(each in infoPDB$uniprot_ptm){
    if(length(each) != 0){
      valid_ptm <- append(valid_ptm, names(infoPDB$uniprot_ptm[i]))
    }
    i <- i+1
  }


  #Select all chosen values by residue num
  chosen_func <- intersect(valid_func, values)
  chosen_ptm <- intersect(valid_ptm, values)
  res_func <- c()
  res_ptm <- c()
  i <- 1
  for(each in infoPDB$uniprot_func){
    if(names(infoPDB$uniprot_func[i]) %in% chosen_func){
      res_func <- append(res_func, infoPDB$uniprot_func[[i]])
    }
    i <- i+1
  }
  i <- 1
  for(each in infoPDB$uniprot_ptm){
    if(names(infoPDB$uniprot_ptm[i]) %in% chosen_ptm){
      res_ptm <- append(res_ptm, infoPDB$uniprot_ptm[[i]])
    }
    i <- i+1
  }

  #If null for selection
  if(is.null(res_func) || is.null(res_ptm)){
    #Visualize the protein
    r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
      r3dmol::m_zoom_to() %>%
      r3dmol::m_add_style(
        style = c(
          m_style_stick(colorScheme="yellowCarbon")
        ),
        sel = m_sel(resi = c(res_func, res_ptm))
      )


  }else{
    #Visualize the protein
    r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
      r3dmol::m_zoom_to() %>%
      r3dmol::m_add_style(
        style = c(
          m_style_stick(colorScheme="yellowCarbon")
        ),
        sel = m_sel(resi = res_func)
      ) %>%
      r3dmol::m_add_style(
        style = c(
          m_style_stick(colorScheme="magentaCarbon")
        ),
        sel = m_sel(resi = res_ptm)
      )  %>% r3dmol::m_add_label(
        text = paste(c("YELLOW", chosen_func, "PINK", chosen_ptm)),
        sel = m_vector3(-6.85, 0.70, 0.30),
        style = m_style_label(inFront = FALSE))
  }

}




#' Visualize polarity with Shiny
#'
#' Function that takes a PDB file and visualizes it with Shiny.
#' Specifically for polarity.
#'
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' to be visualized.
#'
#' @return None - it opens up shiny app with protein structure shown.
#'
#' @examples
#' # Using pdb file available with package
#' #selUniShiny(pdb, "1bm8")
#'
#' @references
#' Chang, W, et al. (2017). Using sliders. Shiny from Rstudio.
#' https://shiny.rstudio.com/articles/sliders.html.
#'
#' @export
#' @import bio3d
#' @import dplyr
#' @import r3dmol
#' @import shiny


selUniShiny <- function(pdbFile, name){
  if (! requireNamespace("shiny", quietly = TRUE)) {
    install.packages("shiny", repos = "http://cran.us.r-project.org")
  }

  # Get info from showInfoPDB
  infoPDB <- showInfoPDB(pdbFile, name)

  # Define UI
  ui <- fluidPage(

    # Title ----
    titlePanel("Select occluded components by Hydrophobicity"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

      # Sidebar for slider options ----
      sidebarPanel(
        # Input:
        checkboxGroupInput("uni", "Select Uniprot Function feature keys:",
                           c('Active', 'Binding', 'DNA', 'Metal', 'Site')),
        checkboxGroupInput("uni2", "Select Uniprot PTM feature keys:",
                           c('Glycosylation', 'Modified'))
      ),

      # Main panel for displaying outputs ----
      mainPanel(
        # Output: 3D model of PDB
        r3dmolOutput("model")
      )
    )
  )

  # Define server logic
  server <- function(input, output) {
    output$model <- renderR3dmol({
      selUniProt(input$uni, pdbFile, name)
      selUniProt(input$uni2, pdbFile, name)
    })
  }

  # Run the shiny server
  shinyApp(ui, server)

}

#[END]
