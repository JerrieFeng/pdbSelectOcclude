#' Visualize occluded AA of PDB
#'
#' Function that takes a PDB file and visualizes it as a 3D model. The user to select which AA to view.
#' This function allows user to interact and view occluded AA.
#'
#' @param index Valid residue index that refers to a AA code found in sequence.
#' The index can be found in showInfoPDB as the a number within the
#' range of code_unique.
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' to be visualized.
#'
#' @return Returns chosen AA, and it displays the 3D model of the protein with
#' the selected AA visible.
#'
#' @examples
#' # Using pdb file available with package
#' selAA(1, pdb, "1bm8")
#'
#'
#' @references
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D.
#' bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.
#'
#' Su, W. (2021). Introduction to r3dmol. r3dmol.
#' https://cran.r-project.org/web/packages/r3dmol/vignettes/r3dmol.html.
#'
#' @export
#' @import bio3d
#' @import dplyr
#' @import r3dmol


selAA <- function(index, pdbFile, name){

  # Get info from showInfoPDB
  infoPDB <- pdbSelectOcclude::showInfoPDB(pdbFile, name)

  #Check if given index is valid
  if( index <= length(infoPDB$AA_info$code_unique) & index > 0){
    cat( paste("Chosen AA is: ", d$AA_info$code_unique[index]) )
  } else{
    stop("The index is out of range. Please make sure the index is a number within
       the possible AA sequences. Refer to `length(infoPDB$AA_info$code_unique)` to
       find the max index that can be inputed.")
  }

  #Visualize the protein
  r3dmol() %>%
    m_add_model(data = m_fetch_pdb(name)) %>%
    m_set_style(style = m_style_cartoon()) %>%
    m_zoom_to() %>%
    m_add_style(
      style = c(
        m_style_cartoon(color="#636efa")
      ),
      sel = m_sel(resn = infoPDB$AA_info$code_unique[index])
    ) %>%
    m_add_style(
      style = c(
        m_style_cartoon(opacity = 0.4)
      ),
      sel = m_sel(resn = infoPDB$AA_info$code_unique[-index])
    ) %>% m_add_label(
      text = infoPDB$AA_info$code_unique[index],
      sel = m_vector3(-6.85, 0.70, 0.30),
      style = m_style_label(inFront = FALSE))

}






#' Visualize occluded AA of PDB with Shiny
#'
#' Function that takes a PDB file and visualizes it with Shiny. It creates a
#' slider for the user to select which AA to view. This function allows user to
#' interact and view occluded AA.
#'
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' to be visualized.
#'
#' @return None - it opens up shiny app with protein structure shown.
#'
#' @examples
#' # Using pdb file available with package
#' selAASlider(pdb, "1bm8")
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


selAASlider <- function(pdbFile, name){

  # Get info from showInfoPDB
  infoPDB <- pdbSelectOcclude::showInfoPDB(pdbFile, name)

  # Define UI
  ui <- fluidPage(

    # Title ----
    titlePanel("Select occluding Amino acid by residue index"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

      # Sidebar for slider options ----
      sidebarPanel(
        # Input:
        sliderInput("resn", "Residue Index:",
                    min = 1, max = length(infoPDB$AA_info$code_unique),
                    value = 1),
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
      pdbSelectOcclude::selAA(input$resn, pdbFile, name)
    })
  }

  # Run the shiny server
  shinyApp(ui, server)

}

#[END]
