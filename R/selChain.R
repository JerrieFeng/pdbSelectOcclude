#' Visualize occluded chains
#'
#' Function that takes a PDB file and visualizes it as a 3D model. The user to select
#' which chain to view. This function allows user to interact and view occluded chains.
#'
#' @param letter Valid chain letter that refers to a chain found in the protein. Run 'showInfoPDB'
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
#' #selChain("A", pdb2, "1SI4", "m_style_stick")
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


selChain <- function(letter, pdbFile, name, style){

  # Get info from showInfoPDB
  infoPDB <- showInfoPDB(pdbFile, name)

  #Check if given letter is valid
  if( letter %in% infoPDB$chain_info$chains){
    cat( paste("Chosen chain is: ", letter) )
  } else{
    stop( paste("Letter out of range! Please pick one of the chains that
                exist in the protein: ", toString(infoPDB$chain_info$chains),
                ". \n Pick by their letter (A,B,C...etc).") )
  }

  if(style == "m_style_stick"){
    style=r3dmol::m_style_stick
  }else{
    style=r3dmol::m_style_sphere
  }

  #Visualize the protein
  r3dmol::r3dmol() %>%
    r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
    r3dmol::m_zoom_to() %>%
    r3dmol::m_add_style(
      style = c(
        style(colorScheme="shapely")
      ),
      sel = m_sel(chain = letter)
    )  %>%
    r3dmol::m_add_style(
      style = c(
        r3dmol::m_style_cartoon(opacity = 0.8, color="#e3e5fc")
      ),
      sel = m_sel(chain = infoPDB$chain_info$chains[-(match(letter, infoPDB$chain_info$chains))])
    ) %>% r3dmol::m_add_label(
      text = paste("Residue Number:",
                   infoPDB$chain_info$residue_num[match(letter, infoPDB$chain_info$chains)], ":",
                   infoPDB$chain_info$residue_num[match(letter, infoPDB$chain_info$chains) + 1],
                   ", Chain:", letter),
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
#' #selChainSlider(pdb2, "1SI4")
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
  infoPDB <- showInfoPDB(pdbFile, name)

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
