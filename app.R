library(shiny)
library(shinythemes)
library(data.table)
library(shinyjs)
library(jsonlite)
load("KEGGs.RData")

options(shiny.maxRequestSize=50*1024^2)

## ui.R

ui <- shinyUI(fluidPage(theme = shinytheme("cerulean"),
                       useShinyjs(),
                       # Application title
                       titlePanel("InChI Unlock"),
                       p("Upload your InChI Keys to get PubChem CIDs, SMILES, and KEGG IDs for your compounds."),
                       
                       # Sidebar with column selection input
                       sidebarLayout(
                         sidebarPanel(
                           fileInput("file1", "Upload a csv file with InChI Keys in one column. If the app times out due to too many compounds, try splitting the file into two.",
                                     multiple = FALSE, accept = c("text/csv", ".csv", 
                                                                  "text/comma-separated-values,text/plain")),
                           selectInput("inchikey_column", "Select InchiKey Column", c("Upload data, then select column")),
                           hr(),
                           actionButton("go", "Start"),
                           downloadButton("download_data", "Download")
                         ),
                         
                         # Main panel will have a completion message or error message
                         mainPanel(textOutput("text")
                                   )
                         )
                       )
             )

## server.R

server <- shinyServer(function(input, output, session) {
  
  # Disable the two buttons
  disable("download_data")
  disable("go")
  
  # Once a csv file is uploaded, populate the inchikey column selection and enable the "go" button
  observe({
    x <- input$file1
    if(is.null(x)){
    }else{
      req(input$file1)
      enable("go")
      file_location <- input$file1$datapath
      dat1 <- data.frame(read.csv(file_location))
      
      updateSelectInput(session, "inchikey_column", 
                        label = "Select InchiKey Column", 
                        choices = colnames(dat1))

    }
  })
  
  # Reload the data in again, and select the already-specified inchikey column
  dat2 <- eventReactive(input$go,{
    showNotification("Reading dataset")
    
    req(input$file1)
    
    file_location <- input$file1$datapath
    df <- data.frame(read.csv(file_location))
    inchi_column_name <- input$inchikey_column
    inchikeys <- df[[inchi_column_name]]
    
    # Note the rows that don't have NA values
    inchikeys_not_missing <- which(!is.na(inchikeys))
    
    # Get KEGGs from the KEGGs.RData file that is part of the app
    showNotification("Getting KEGG IDs")
    keggs <- merge(data.table(inchikey = inchikeys), data.table(kegg_list), by.x = "inchikey", by.y = "InChIKey", all.x = TRUE, all.y = FALSE, sort = FALSE)[, "KEGG"]
    showNotification("KEGG ID translation complete")
      
    # Create a blank list for the PubChem CIDs
    cids <- list()
    showNotification("Getting the PubChem CIDs")
    withProgress(message = 'InChI Key to PubChem CID translation in progress', value = 0, {
        
    # If there are fewer than 100 inchikeys, put them all into PubChem ID exchange at once. If there are more than 100, do them in blocks of 100.
    # This is so that the PubChem ID exchange doesn't get overloaded and blocks the app from working. This uses the PubChem PUG REST API
    if(length(inchikeys_not_missing)<100){
        trans_block <- c(1,length(inchikeys_not_missing))
      }else{
        trans_block <- c(seq(1,length(inchikeys_not_missing),100), length(inchikeys_not_missing))
      }
        
    for(i in 2:length(trans_block)){
        cids[[i-1]] <- tryCatch({
          CIDS <- fromJSON(url(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", paste0(inchikeys[inchikeys_not_missing[trans_block[i-1]:trans_block[i]]], collapse = ","),"/property/InChIKey/JSON")))
          CIDS$PropertyTable$Properties
          
        }, error = function(er){
          print(er)
          return(NA)
        })
      }
        
    # Put the PubChem CIDs into a vector
    cids_data <- do.call('rbind', cids)
    # If there are any PubChem CIDs, merge with the inchikeys
    if(length(cids_data) > 0){
        cids_data <- cids_data[!duplicated(cids_data$InChIKey), ]
        cids <- merge(data.table(inchikey = inchikeys), data.table(cids_data), by.x = "inchikey", by.y = "InChIKey", all.x = TRUE, all.y = FALSE, sort = FALSE)[, "CID"][[1]]
        # If there are no PubChem CIDs, make a data frame of NA values 
      } else {
        cids <- data.frame(CID = character(length(inchikeys_not_missing)))
        cids$CID <- rep(NA, times = length(inchikeys_not_missing))
      }
    })
    showNotification("PubChem CID translation complete")
    
    # Create a blank list for the SMILES codes
    smiles <- list()
    showNotification("Getting the SMILES")
    withProgress(message = 'InChI Key to SMILES translation in progress', value = 0, {
    
    for(i in 2:length(trans_block)){
      smiles[[i-1]] <- tryCatch({
        SMILES <- fromJSON(url(paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", paste0(inchikeys[inchikeys_not_missing[trans_block[i-1]:trans_block[i]]], collapse = ","),"/property/CanonicalSMILES,InChIKey/JSON")))
        SMILES$PropertyTable$Properties
          
        }, error = function(er){
          print(er)
          return(NA)
        })
      }
        
    # Put the SMILES codes into a vector
    smiles_data <- do.call('rbind', smiles)
    smiles_data <- smiles_data[, !colnames(smiles_data) %in% "CID"]
    # If there are any SMILES codes, merge with the inchikeys
    if(length(smiles_data) > 0){
        smiles_data <- smiles_data[!duplicated(smiles_data$InChIKey), ]
        smiles <- merge(data.table(inchikey = inchikeys), data.table(smiles_data), by.x = "inchikey", by.y = "InChIKey", all.x = TRUE, all.y = FALSE, sort = FALSE)[, "CanonicalSMILES"][[1]]
        # If there are no SMILES codes, make a data frame of NA values 
    } else {
        smiles <- data.frame(smiles = character(length(inchikeys_not_missing)))
        smiles$smiles <- rep(NA, times = length(inchikeys_not_missing))
      }
    })
    showNotification("SMILES translation complete")
      
    # Merge the original uploaded data, the KEGGs, SMILES, and CIDs into a single data frame
    result <- cbind(df, keggs, smiles, cids)
    # Name the columns nicely
    colnames(result)[(ncol(result)-2):ncol(result)] <- c("KEGG ID", "SMILES", "PubChem CID")
    # Write the data frame to a csv file
    fwrite(result, "Translated_IDs.csv", row.names = FALSE)
    # Enable the download button for file download
    enable("download_data")
        
    return("Compound translation is complete. Please click the 'Download' button to download the csv file.")
    
})

# Put the text in line 159 on the main panel of the app
output$text <- renderText({dat2()})

# Create the zip file with the csv output file in it
output$download_data <- downloadHandler(
  filename = function() {"Translation.zip"},
  content = function(fname) {zip(fname, c("Translated_IDs.csv"))},
  contentType = "application/zip")

})  

shinyApp(ui = ui, server = server)