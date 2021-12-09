#' Visualize occluded chains
#'
#' Function that takes a PDB file and visualizes it as a 3D model. The user to select
#' which chain to view. This function allows user to interact and view occluded chains.
#'
#' @param index Valid index that refers to a chain found in the protein. Run 'showInfoPDB'
#' to get all possible chains in chain_info
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' to be visualized.
#'
#' @return Returns output-text stating the chain, and it displays the 3D model
#' of the protein with the selected chain visible.
#'
#' @examples
#' # Using pdb2 file available with package
#' selChain(1, pdb2, "1SI4")
#'
#'
#' @references
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D.
#' bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.
#'
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). chain.pdb: Find Possible PDB Chain Breaks.
#' bio3D. https://rdrr.io/cran/bio3d/man/chain.pdb.html.
#'
#' Su, W. (2021). Introduction to r3dmol. r3dmol.
#' https://cran.r-project.org/web/packages/r3dmol/vignettes/r3dmol.html.
#'
#' @export
#' @import bio3d
#' @import dplyr
#' @import r3dmol


selChain <- function(index, pdbFile, name){

  # Get info from showInfoPDB
  infoPDB <- pdbSelectOcclude::showInfoPDB(pdbFile, name)

  #Check if given index is valid
  if( index <= length(infoPDB$chain_info$chains) & index > 0){
    cat( paste("Chosen chain is: ", infoPDB$chain_info$chains[index]) )
  } else{
    stop( paste("Index number out of range! Please pick one of the chains that
                exist in the protein: ", toString(infoPDB$chain_info$chains),
                ". \n Pick by their index number (1,2,3...etc).") )
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
      sel = m_sel(chain = infoPDB$chain_info$chains[index])
    )  %>%
    m_add_style(
      style = c(
        m_style_cartoon(opacity = 0.4)
      ),
      sel = m_sel(chain = infoPDB$chain_info$chains[-index])
    ) %>% m_add_label(
      text = paste("Residue Number:",
                   infoPDB$chain_info$residue_num[index], ":",
                   infoPDB$chain_info$residue_num[index+1],
                   ", Chain:", infoPDB$chain_info$chains[index]),
      sel = m_vector3(-6.89, 0.75, 0.35),
      style = m_style_label(inFront = FALSE))

}






#' Visualize occluded chains with Shiny
#'
#' Function that takes a PDB file and visualizes it with Shiny. It creates a
#' slider for the user to select which chain to view. This function allows user to
#' interact and view occluded protein chains.
#'
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' to be visualized.
#'
#' @return None - it opens up shiny app with protein structure shown.
#'
#' @examples
#' # Using pdb2 file available with package
#' selChainSlider(pdb2, "1SI4")
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


selChainSlider <- function(pdbFile, name){

  # Get info from showInfoPDB
  infoPDB <- pdbSelectOcclude::showInfoPDB(pdbFile, name)

  # Define UI
  ui <- fluidPage(

    # Title ----
    titlePanel("Select occluding Chains by corresponding index"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

      # Sidebar for slider options ----
      sidebarPanel(
        # Input:
        sliderInput("resi", "Chain index:",
                    min = 1, max = length(infoPDB$chain_info$chains),
                    value = 1,
                    step = 1),
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
      pdbSelectOcclude::selChain(input$resi, pdbFile, name)
    })
  }

  # Run the shiny server
  shinyApp(ui, server)

}

#[END]
