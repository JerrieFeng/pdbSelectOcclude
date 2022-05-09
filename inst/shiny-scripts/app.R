library(shiny)
library(dplyr)


merge <- function(pdbFile, name, styleSel, codes=NULL, polarity=NULL, letter="ALL" ){
  # Get info from showInfoPDB
  infoPDB <- pdbSelectOcclude::showInfoPDB(pdbFile, name)

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

  #Default option - all
  if(letter == "ALL"){
    letters <- infoPDB$chain_info$chains
  }else{
    letters <- letter
  }

  #if all default
  if(all(infoPDB$AA_info$code_unique %in% res) && letter == "ALL"){
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
          r3dmol::m_style_cartoon(opacity = 0.8, color="#e3e5fc")
        ),
        sel = m_sel(chain = infoPDB$chain_info$chains[-(match(letters, infoPDB$chain_info$chains))])
      ) %>% r3dmol::m_add_label(
        text = infoPDB$name,
        sel = m_vector3(-6.89, 0.75, 0.35),
        style = m_style_label(inFront = FALSE))
  }


}



#runn(pdb2, "1SI4")

runn <- function(pdbFile, name){
  infoPDB <- pdbSelectOcclude::showInfoPDB(pdbFile, name)

  # Define UI
  ui <- fluidPage(

    titlePanel("Select :) "),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

      # Sidebar for slider options ----
      sidebarPanel(
        # Input:
        radioButtons("style", "Select a style of selection:", c("m_style_stick", "m_style_sphere"), selected="m_style_stick"),
        radioButtons("chain", "Select a Chain:", c("ALL", infoPDB$chain_info$chains), selected="ALL"),

        checkboxGroupInput("polarity", "Select polarity:", c('polar', 'nonpolar', 'positive', 'negative')),
        checkboxGroupInput("aacodes", "Select Amino Acid:", infoPDB$AA_info$code_unique)
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
      merge(pdbFile, name, input$style, input$aacodes, input$polarity, input$chain)

    })


  }


  # Run the shiny server
  shinyApp(ui, server)
}







