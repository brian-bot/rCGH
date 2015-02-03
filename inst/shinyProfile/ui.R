# require(synapseClient)
# geneList <- read.delim(synGet("syn2347502")@filePath, stringsAsFactors=FALSE)$Gene

shinyUI(
    pageWithSidebar(
            
      headerPanel("Interactive CGH Viewer"),
      
      sidebarPanel(
#        includeCSS("shinySafir.css"),
        tags$head( tags$link(rel="stylesheet", type="text/css", href="cghViewer.css") ),

        
        h4('Gene symbol'),
        textInput("geneSymb", '', 'NONE'),
        tags$hr(),
                      
        h4('Show chromosome'),
        selectInput(inputId = "chr", label = "", choices = c('All', 1:23), selected = 'All'),
        tags$hr(),
        
        sliderInput("center", "Recenter profile", min=-1.5, max=1.5, value=0, step = .1),
        sliderInput("Ymax", "Rescale max(y)", min=.1, max=1, value=1, step=.1),
        sliderInput("Ymin", "Rescale min(y)", min=.1, max=1, value=1, step=.1),
        sliderInput("gain", "Gain threshold (Log2ratio)", min=0, max=2, value=.5, step = .25),
        sliderInput("loss", "Loss threshold (Log2ratio)", min=-2, max=0, value=-.5, step = .25),
        tags$hr(),
        p(),
        withTags(div(class='row-fluid', style="margin-top: 20px;", align="left",
                     a("@Contact us", style="font-size: 14px;", href="mailto:frederic.commo@gustaveroussy.fr?Subject=rCGH%20Viewer", target="_top")
                    )
                )
        ),

      mainPanel(
        tabsetPanel(
          tabPanel("CGH profile",
                   plotOutput("Profile", width = "100%", height = "100%"),
                   tags$hr(),
                   tableOutput("geneSummary")
                   ),
          tabPanel("CGH table",
                  h4(textOutput("tableTitle1"), align="center"),
                  h4(textOutput("tableTitle2"), align="center"),
                  tags$hr(),
                  dataTableOutput("fullTable")
                   )
          )
        )
      )
    )