library(shiny)
library(DT)
library(dplyr)
library(zoo)
library(htmlwidgets)
library(jsonlite)
library(bslib)
ui <- fluidPage(

  theme = bslib::bs_theme(bootswatch="flatly", bg = "rgb(38, 38, 38)", fg = "#fff", primary = "rgb(181, 72, 98)", secondary = "rgb(82, 165, 171)"),
  titlePanel(
    #Instructions text at the top of the ui
    'CCP4 Contacts Analysis'),
  sidebarLayout(
    sidebarPanel(style="overflow-y:auto; overflow-x:auto",
    h5("Works with the following CCP4 parameters:"),
    
    #Embedded images, file in the 'www' folder
  actionButton("ccp4image", "Show parameters"),
    h5("Inter chain is shown here, but inter residue should also work."),
    
    #Instructions continued
    h2("Instructions:", style = "color:#C24450"),
    h5("1: Upload CCP4 log output in .txt format. Must be in plain text format. TextEdit users, Format -> Make Plain Text"),
    h5("2: Enter maximum distances desired for H-bond, salt bridge, and Van der Waals annotations:"),
    actionButton("instructions", "Show detailed instructions"),
    h4("'bonding_type' criteria:", style = "color:#C24450"),
    h5("- 'HBOND' (hydrogen bond) if the following hydrogen donor and acceptor is present in both 'donor_atom' and 'acceptor_atom' AND 'contact_distance' < input:"),
    actionButton("hbondtable", "Show table"),
    h5("- 'Salt' (salt bridge) if the following is present in 'donor_atom' AND in 'acceptor_atom' AND 'contact_distance' < input"),
    actionButton("salttable", "Show table"),  
    h5("- 'VDW' (Van der Waals) if C is present in 'donor_atom' AND C is present in 'acceptor_atom' AND 'contact_distance' < input"),
    tags$a(href="http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/#MLformula", "Click here for protein atom nomenclature.", target="_blank"),
    h5("DNA atom nomenclature is shown below:"),
    actionButton("DNAatomnames", "Show nomenclature"),
    
    h4("Exporting selected results:", style = "color:#C24450"),
    h5("- Important: Type desired file name prior to selecting rows or filtering. Entering or changing the file name after the fact resets the table."),
    h5("- Select relevant rows and click on the desired format in the top left."),
    h5("- If nothing is selected, the function exports the entire table."),
    h4('Legend:', style = "color:#C24450"),
    actionButton("tablelegend", "Show legend")),
    
    mainPanel(
    h2("Define parameters:", style = "color:#C24450"),
    #Textbox to rename file if exporting as csv/xslx
    textInput("filename","If exporting as CSV/XSLX, enter name for file"),
    numericInput(inputId = "hbonddistanceinput", label = "Define maximum hydrogen bond distance", 
                 value = 3.50),
    numericInput(inputId = "saltdistanceinput", label = "Define maximum salt bridge distance", 
                 value = 6),
    numericInput(inputId = "vdwdistanceinput", label = "Define maximum VDW distance", 
                 value = 4.50),
    #File input box, only .txt files accepted. Uploaded file is internally named 'file1'
    h2("Upload your data:", style = "color:#C24450"),
    options(shiny.maxRequestSize = 500*1024^2),
    fileInput("file1", "Upload TXT File",         
              multiple = FALSE,
              accept = c("text",
                         ".txt")),
   
    #Code to display the DT table in the UI    
     fluidRow(
      DT::DTOutput("table1")),
    shinyjs::useShinyjs(),
    tags$style(".dtsp-name {
    color: black;}")
    )))

server <- function(input, output) {
  
  #Code that manipulates the data. 'datasetInput' is what the final curated data frame is named
  datasetInput <- eventReactive(input$file1, {
    #Names the input file (file1) file2
    req(input$hbonddistanceinput)
    req(input$saltdistanceinput)
    req(input$vdwdistanceinput)
    file2 <- input$file1 
    file3 <- readr::read_fwf(file2$datapath, readr::fwf_widths(c(6, 4, 2, 6, 5, 6, 5, 2, 5, 8, 5, 5, 100000000000000004))) #Converts the uploaded file to a data frame, input treated as a fixed width format file, list indicates the width that the input data should be separated by
    file3$X1[file3$X1 == ""] <- NA #Converts empty cells in Column 1 to NA. 'read_fwf' names columns as X1, X2, etc. initially
    file3$X1 <- na.locf(file3$X1) #Fills in cells in X1(Column1) that contain NA with the contents of the cell directly above them. 
    file3$X2[file3$X2 == ""] <- NA
    file3$X2 <- na.locf(file3$X2)
    file3$X3[file3$X3 == ""] <- NA
    file3$X3 <- na.locf(file3$X3)
    file3$X4[file3$X4 == ""] <- NA
    file3$X4 <- na.locf(file3$X4)
    file3$X13 <- NULL #In case there's a lot of text above/below the data, which would create extra columns, convert those columns to NULL which would not be displayed
    file3$X14 <- NULL
    
    #Searches the uploaded file for specific characters in specific columns and renames Column10 (X10) with the annotated bond type if criteria met
    
    file3$X10[grepl("C", file3$X4) & grepl("C", file3$X9) & file3$X11 < as.numeric(input$vdwdistanceinput)] <- "VDW"
    
    file3$X10[grepl("Arg", file3$X1) & grepl("NE|NH1|NH2|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Asn", file3$X1) & grepl("OD1|ND2|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Asn", file3$X6) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("ND2|OD1|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Asp", file3$X1) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("OD1|OD2|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) &  file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Gln", file3$X1) & grepl("NE2|OE1|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Gln", file3$X6) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("NE2|OE1|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Glu", file3$X1) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1|", file3$X4) & grepl("OE1|OE2|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("His", file3$X1) & grepl("ND1|NE2|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("His", file3$X6) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("ND1|NE2|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Lys", file3$X1) & grepl("NZ|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Ser", file3$X1) & grepl("OG|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Ser", file3$X6) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("OG|CG2|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Thr", file3$X1) & grepl("OG1|CG2|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1|CG2", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Thr", file3$X6) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("OG1|CG2|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Trp", file3$X1) & grepl("NE1|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Tyr", file3$X1) & grepl("OH|\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Tyr", file3$X6) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("OH|\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    
    #backbone hbond acceptors/donors
    file3$X10[grepl("\\bN\\b", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1", file3$X9) & !grepl("\\bP\\b", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1", file3$X4) & grepl("\\bO\\b", file3$X9) & !grepl("\\bP\\b", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    
    
    #DNA/RNA
    #all donors
    "NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|N|O2|O2*|OD1|OE1|N6|N1|N4|N3"
    
    #all acceptors
    "OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|O|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1|N1|N3|N7|O6|O2|O4|O5|O3"
    
    file3$X10[grepl("Da", file3$X1) & grepl("N6|OP1|OP2", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1|N1|N3|N7|O6|O2|O4|O5|O3", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Da", file3$X6) & !grepl("Asp|Glu", file3$X1) & !grepl("\\bO\\b", file3$X4) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1|N6|N1|N4|N3", file3$X4) & grepl("N1|N3|N7", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Dt", file3$X1) & grepl("N3|OP1|OP2", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1|N1|N3|N7|O6|O2|O4|O5|O3", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Dt", file3$X6) & !grepl("Asp|Glu", file3$X1) & !grepl("\\bO\\b", file3$X4) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1|N6|N1|N4|N3", file3$X4) & grepl("O2|O4", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Dc", file3$X1) & grepl("N4|OP1|OP2", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1|N1|N3|N7|O6|O2|O4|O5|O3", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Dc", file3$X6) & !grepl("Asp|Glu", file3$X1) & !grepl("\\bO\\b", file3$X4) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1|N6|N1|N4|N3", file3$X4) & grepl("N3|O2", file3$X9) & !grepl("C", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Dg", file3$X1) & grepl("N1|OP1|OP2", file3$X4) & grepl("OD1|OD2|OE1|OE2|ND1|ND2|NE2|OG|OG1|OH|\\bO\\b|O4'|O2'|O4*|O2*|OP1|OP2|OG|OG1|N1|N3|N7|O6|O2|O4|O5|O3", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    file3$X10[grepl("Dg", file3$X6) & !grepl("Asp|Glu", file3$X1) & !grepl("\\bO\\b", file3$X4) & grepl("NE|NH1|NH2|ND2|ND1|NE2|NZ|OG|OG1|NE1|OH|\\bN\\b|OD1|OE1|N6|N1|N4|N3", file3$X4) & grepl("O6|N7|N3", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- "HBOND"
    #exceptions
    file3$X10[grepl("Dg", file3$X6) & grepl("N1", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- NA
    file3$X10[grepl("Da", file3$X1) & grepl("N1", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- NA
    
    file3$X10[grepl("Dt", file3$X6) & grepl("N3", file3$X9) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- NA
    file3$X10[grepl("Da|Dg|Dc", file3$X1) & grepl("N3", file3$X4) & file3$X11 <= as.numeric(input$hbonddistanceinput)] <- NA
    
    
    file3$X10[grepl("OE1|OE2|OD1|OD2|OP1|OP2|\\bP\\b", file3$X4) & grepl("Da|Dt|Dc|Dg|Glu|Asp", file3$X1) & grepl("NZ|NH1|NH2|NE", file3$X9) & grepl("Lys|Arg", file3$X6) & file3$X11 < as.numeric(input$saltdistanceinput)] <- "Salt"
    file3$X10[grepl("NZ|NH1|NH2|NE", file3$X4) & grepl("Lys|Arg", file3$X1) & grepl("OE1|OE2|OD1|OD2|OP1|OP2|\\bP\\b", file3$X9) & grepl("Da|Dt|Dc|Dg|Glu|Asp", file3$X6) & file3$X11 < as.numeric(input$saltdistanceinput)] <- "Salt"
    
    #Deletes cells in these columns that are not numeric
    file3$X2 <- as.numeric(gsub('[a-zA-Z]', '', file3$X2))
    file3$X7 <- as.numeric(gsub('[a-zA-Z]', '', file3$X7))
    file3$X11 <- as.numeric(gsub('[a-zA-Z]', '', file3$X11)) 
    file3
  })
  
  
  #DT package-specific code
  output$table1 <- DT::renderDT(
    datasetInput(), #Our dataframe name 
    
    #Rename columns
    colnames = c('donor_mol', 'donor_pos_range', 'donor_chainID', 'donor_atom', ' ', 'acceptor_mol', 'acceptor_pos_range', 'acceptor_chainID','acceptor_atom', 'bonding_type', 'contact_distance', 'angle'), #Names columns
    
    #DT extensions (https://rstudio.github.io/DT/extensions.html), plus options for each extension that must all be placed into 'options = list()'
    extensions = c('Buttons', 'Select', 'SearchPanes', 'Scroller'), 
    server = FALSE, 
    options = list(
      search = list(regex = TRUE, caseInsensitive = FALSE, search = '\\.\\.\\.'),
      ordering = F,
      deferRender = TRUE,
      pageLength = 100, dom = 'Bfrltip',
      select = list(style = 'multiple', items = 'row'),
      buttons = list('searchPanes', 'copy', list(extend = 'csv', filename = input$filename), list(extend = 'excel', filename = input$filename), 'selectRows', 'selectNone'),
      columnDefs = list(list(searchPanes = list(show = TRUE), targets = 1:5))),
    rownames = FALSE,
    selection = 'none')
  #UI functions
  observeEvent(input$instructions, {
    showModal(modalDialog(
      img(src='instructions.png', height="100%", width="100%"),
      title='Instructions', footer=NULL, size='l', easyClose=TRUE, fade=TRUE))})
  observeEvent(input$ccp4image, {
    showModal(modalDialog(
      img(src='ccp4input.png', height="100%", width="100%"),
      title='Parameters', footer=NULL, size='l', easyClose=TRUE, fade=TRUE))})
  observeEvent(input$hbondtable, {
    showModal(modalDialog(
      img(src='hbondtable.png', height="100%", width="100%"),
      title='Table of interactions', footer='Atoms in red are to capture hydrogen bonds due to model ambiguity.', size='l', easyClose=TRUE, fade=TRUE))})
  observeEvent(input$salttable, {
    showModal(modalDialog(
      img(src='salttable.png', height="100%", width="100%"),
      title='Parameters', footer=NULL, size='l', easyClose=TRUE, fade=TRUE))})
  observeEvent(input$DNAatomnames, {
    showModal(modalDialog(
      img(src='nucleicAcidAtoms.png', height="100%", width="100%"),
      title='DNA atom', footer=NULL, size='l', easyClose=TRUE, fade=TRUE))})
  observeEvent(input$tablelegend, {
    showModal(modalDialog(
      img(src='colnames.png', height="100%", width="100%"),
      title='Column legend', footer=NULL, size='l', easyClose=TRUE, fade=TRUE))})
  
}

shinyApp(ui, server)
