#' Compile all functions into Shiny
#'
#' Function that takes a PDB file and visualizes it as a 3D model.
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
#' # merge(pdb2, "1SI4", "m_style_stick", "VAL", NULL, "ALL", "Metal")
#' # runModel(pdb2, "1SI4")
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
#' @export merge
#' @export runModel
#'
#' @import bio3d
#' @import dplyr
#' @import r3dmol
#' @import shiny



#Merge all functions into one
merge <- function(pdbFile, name, styleSel, codes=NULL, polarity=NULL, letter="ALL", values=NULL){

  if (! requireNamespace("shiny", quietly = TRUE)) {
    install.packages("shiny", repos = "http://cran.us.r-project.org")
  }

  # Get info from showInfoPDB
  infoPDB <- showInfoPDB(pdbFile, name)

  #selection colour
  if(styleSel == "m_style_stick"){
    styleSel=r3dmol::m_style_stick
  }else{
    styleSel=r3dmol::m_style_sphere
  }


  # If no polarity or AA selected
  valid <- c()
  # POLARITY
  if(is.null(codes) && is.null(polarity)){
    res <- infoPDB$AA_info$code_unique
  }else if(is.null(codes) && !is.null(polarity)){
    #Check if given string is valid
    if(polarity != "ALL"){
      for(each in polarity){
        temp <- hphobScale(each)
        valid <- append(valid, temp)
      }
    }
    if(polarity == "ALL"){
      valid <- infoPDB$AA_info$code_unique
    }
    #Select all hphob in the chosen protein
    hphob <- intersect(valid, infoPDB$AA_info$code_unique)
    res <- hphob
  }else if(!is.null(codes) && is.null(polarity)){
    res <- codes
  }else{
    #Check if given string is valid
    if(polarity != "ALL"){
      for(each in polarity){
        temp <- hphobScale(each)
        valid <- append(valid, temp)
      }
    }
    #Select all hphob in the chosen protein
    hphob <- intersect(valid, infoPDB$AA_info$code_unique)

    # INTERSECT - (polarity INTERSECT AA) within Chains
    res <- intersect(hphob, codes)
  }


  # UNIPROT
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


  #Default option - all
  if(letter == "ALL"){
    letters <- infoPDB$chain_info$chains
  }else{
    letters <- letter
  }

  #if all default
  if(all(infoPDB$AA_info$code_unique %in% res) && letter == "ALL" && all(is.null(res_func)) && all(is.null(res_ptm)) ){
    # VIZUALIZE
    r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
      r3dmol::m_zoom_to() %>%
      r3dmol::m_add_style(
        style = c(
          r3dmol::m_style_cartoon(color="#e3e5fc")
        ),
        sel = m_sel(resn = res)
      ) %>% r3dmol::m_add_label(
        text = infoPDB$name,
        sel = m_vector3(-6.89, 0.75, 0.35),
        style = m_style_label(inFront = FALSE))

  }else if(length(res_func)!=0 && length(res_ptm)==0){
    # VIZUALIZE
    r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
      r3dmol::m_zoom_to() %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="shapely")
        ),
        sel = m_sel(chain = letters, resn = res)
      )  %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="yellowCarbon")
        ),
        sel = m_sel(chain = letters, resi = res_func)
      )  %>%
      r3dmol::m_add_style(
        style = c(
          r3dmol::m_style_cartoon(opacity = 0.8, color="#e3e5fc")
        ),
        sel = m_sel(chain = infoPDB$chain_info$chains[-(match(letters, infoPDB$chain_info$chains))])
      ) %>% r3dmol::m_add_label(
        text = infoPDB$name,
        sel = m_vector3(-6.89, 0.75, 0.35),
        style = m_style_label(inFront = FALSE))

  }else if(length(res_func)==0 && length(res_ptm)!=0){
    # VIZUALIZE
    r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
      r3dmol::m_zoom_to() %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="shapely")
        ),
        sel = m_sel(chain = letters, resn = res)
      ) %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="blueCarbon")
        ),
        sel = m_sel(chain = letters, resi = res_ptm)
      ) %>%
      r3dmol::m_add_style(
        style = c(
          r3dmol::m_style_cartoon(opacity = 0.8, color="#e3e5fc")
        ),
        sel = m_sel(chain = infoPDB$chain_info$chains[-(match(letters, infoPDB$chain_info$chains))])
      ) %>% r3dmol::m_add_label(
        text = infoPDB$name,
        sel = m_vector3(-6.89, 0.75, 0.35),
        style = m_style_label(inFront = FALSE))

  }else if(length(res_func)==0 && length(res_ptm)==0){
    # VIZUALIZE
    r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
      r3dmol::m_zoom_to() %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="shapely")
        ),
        sel = m_sel(chain = letters, resn = res)
      ) %>%
      r3dmol::m_add_style(
        style = c(
          r3dmol::m_style_cartoon(opacity = 0.8, color="#e3e5fc")
        ),
        sel = m_sel(chain = infoPDB$chain_info$chains[-(match(letters, infoPDB$chain_info$chains))])
      ) %>% r3dmol::m_add_label(
        text = infoPDB$name,
        sel = m_vector3(-6.89, 0.75, 0.35),
        style = m_style_label(inFront = FALSE))

  }else{
    # VIZUALIZE
    r3dmol::r3dmol() %>%
      r3dmol::m_add_model(data = r3dmol::m_fetch_pdb(name)) %>%
      r3dmol::m_zoom_to() %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="shapely")
        ),
        sel = m_sel(chain = letters, resn = res)
      )  %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="yellowCarbon")
        ),
        sel = m_sel(chain = letters, resi = res_func)
      )  %>%
      r3dmol::m_add_style(
        style = c(
          styleSel(colorScheme="magentaCarbon")
        ),
        sel = m_sel(chain = letters, resi = res_ptm)
      ) %>%
      r3dmol::m_add_style(
        style = c(
          r3dmol::m_style_cartoon(opacity = 0.8, color="#e3e5fc")
        ),
        sel = m_sel(chain = infoPDB$chain_info$chains[-(match(letters, infoPDB$chain_info$chains))])
      ) %>% r3dmol::m_add_label(
        text = infoPDB$name,
        sel = m_vector3(-6.89, 0.75, 0.35),
        style = m_style_label(inFront = FALSE))
  }


}








#Running the 3D model with Shiny

runModel <- function(pdbFile, name){
  infoPDB <- showInfoPDB(pdbFile, name)
  valids <- checkValid(infoPDB)

  # Define UI
  ui <- fluidPage(

    titlePanel("Hello! Select the different features "),

    # Sidebar layout with input and output definitions ----
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: 3D model of PDB
      r3dmolOutput("model")
    ),

    # Sidebar for slider options ----
    fluidRow(
      # Input:
      column(3,
             radioButtons("style", "Select a style of selection:", c("m_style_stick", "m_style_sphere"), selected="m_style_stick"),
             radioButtons("chain", "Select a Chain:", c("ALL", infoPDB$chain_info$chains), selected="ALL")
      ),
      column(4, offset = 1,
             checkboxGroupInput("polarity", "Select polarity:", c('polar', 'nonpolar', 'positive', 'negative')),
             checkboxGroupInput("aacodes", "Select Amino Acid:", infoPDB$AA_info$code_unique)
      ),
      column(4,
             checkboxGroupInput("values", "Select Uniprot feature keys (Yellow = Funtional Sites ; Blue = PTM sites):",
                                valids)
      )


    )




  )


  # Define server logic
  server <- function(input, output) {

    output$model <- renderR3dmol({
      merge(pdbFile, name, input$style, input$aacodes, input$polarity, input$chain, input$values)

    })

  }


  # Run the shiny server
  shinyApp(ui, server)
}

#[END]
