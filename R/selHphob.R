#' Visualize hydrophobicity
#'
#' Function that takes a PDB file and visualizes it as a 3D model. The user to
#' select on a scale to view hydrophobicity. Uses Kyte & Doolittle scale.
#'
#' @param polarity Enter 'nonpolar' to select all nonpolar AA, 'polar' to select all polar/uncharged AA,
#' 'positive' to select all positively charged AA, 'negative' to select all negatively charged AA.
#' Make sure to add quotation marks.
#' @param pdbFile A PDB file, could be downloaded from PDB online
#' @param name the 4-letter PDB codes/identifiers of the PDB file
#' to be visualized.
#'
#' @return Displays the 3D model of the protein with the selected hydrophobicity.
#'
#' @examples
#' # Using pdb file available with package
#' # ?idpr
#' #selHphob("polar", pdb, "1bm8")
#'
#' @references
#' Yao, X., G. Scarabelli, L. Skjaerven, and B. Grant (2020). Protein Structure Networks with Bio3D.
#' bio3D. http://thegrantlab.org/bio3d/articles/online/cna_vignette/cna_vignette.spin.html#references-1.
#'
#' McFadden, W.M. (2020). Charge and Hydropathy Vignette.
#' https://rdrr.io/bioc/idpr/f/vignettes/chargeHydropathy-vignette.Rmd.
#'
#' @export hphobScale
#' @export selHphob
#'
#' @import bio3d
#' @import dplyr
#' @import r3dmol



hphobScale <- function(polarity){
  nonpolar <- c("ILE", "PHE", "VAL", "LEU", "TRP", "MET", "ALA", "GLY", "PRO")
  polar <- c("CYS", "TRY", "THR", "SER", "ASN", "GLN") #also uncharge
  posCharge <- c("HIS", "LYS", "ARG")
  negCharge <- c("GLU", "ASP")

  if(polarity == "nonpolar"){
    return(nonpolar)
  }else if(polarity == "polar"){
    return(polar)
  } else if(polarity == "positive"){
    return(posCharge)
  } else if(polarity == "negative"){
    return(negCharge)
  }else{
    stop("Please enter 'nonpolar' to select all nonpolar AA, 'polar' to select all polar/uncharged AA,
         'positive' to select all positively charged AA, 'negative' to select all negatively charged AA.
         Make sure to add quotation marks.")
  }

}

selHphob <- function(polarity, pdbFile, name){

  # Get info from showInfoPDB
  infoPDB <- showInfoPDB(pdbFile, name)

  #Check if given string is valid
  valid <- c()
  for(each in polarity){
    temp <- hphobScale(each)
    valid <- append(valid, temp)
  }

  #Select all hphob in the chosen protein
  hphob <- intersect(valid, infoPDB$AA_info$code_unique)

  #Get non-hphob values
  opposite <- infoPDB$AA_info$code_unique
  opposite <- opposite[opposite %in% hphob == FALSE]

  #Visualize the protein
  r3dmol::r3dmol() %>%
    r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
    r3dmol::m_zoom_to() %>%
    r3dmol::m_add_style(
      style = c(
        m_style_stick(colorScheme="shapely")
      ),
      sel = m_sel(resn = hphob)
    ) %>%
    r3dmol::m_add_style(
      style = c(
        m_style_cartoon(color = "#e3e5fc", opacity = 1)
      ),
      sel = m_sel(resn = opposite)
    ) %>% r3dmol::m_add_label(
      text = paste(c(polarity, ":", hphob)),
      sel = m_vector3(-6.85, 0.70, 0.30),
      style = m_style_label(inFront = FALSE))

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
#' #selHphobSlider(pdb, "1bm8")
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


selHphobSlider <- function(pdbFile, name){
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
        checkboxGroupInput("polarity", "Select polarity:", c('polar', 'nonpolar', 'positive', 'negative'))
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
      selHphob(input$polarity, pdbFile, name)
    })
  }

  # Run the shiny server
  shinyApp(ui, server)

}

#[END]
