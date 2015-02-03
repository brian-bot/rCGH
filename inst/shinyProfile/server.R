###########################
###########################
require(ggplot2)
require(grid)

source('helpers.R')
cat("Loading files...")
load(file.path(getwd(), "extdata/hg19.rda"))
load(file.path(getwd(), "extdata/geneDB.rda"))
segTable <- readRDS(file.path(getwd(), "data/st.rds"))
geneTable <- readRDS(file.path(getwd(), "data/bg.rds"))
#hg19 <- readRDS(file.path(getwd(), "extdata/hg19.rds"))
#geneDB <- readRDS(file.path(getwd(), "extdata/geneDB.rds"))
cat("Done.\n")
###########################

shinyServer(function(input, output) {
  
  gene <- reactiveValues(symbol=character())
  observe({gene$symbol <- toupper(input$geneSymb)})
      
  reCenterSeg <- reactive({
    seg <- segTable
    seg$seg.med <- seg$seg.med + input$center
    return(seg)
    })

  reCenterGenes <- reactive({
    geneTable$Log2Ratio <- geneTable$Log2Ratio + input$center
    return(geneTable)
  })

  createCGHplot <- reactive({
    seg <- reCenterSeg()
    gPlot <- .mainPlot(seg)
    gPlot <- .addSegments(gPlot, seg, input$chr, input$gain, input$loss)
    gPlot <- .updateYscale(gPlot, input$Ymin, input$Ymax)
    gPlot <- .updateXscale(gPlot, input$chr, hg19)
    gPlot <- .addChr(gPlot, input$chr, hg19)
    gPlot <- .addTitle(gPlot, unique(segTable$ID), input$gain, input$loss)

    if(!gene$symbol %in% c('NONE', '')){
      geneAnnot <- try(.geneOfInt(seg, gene$symbol, geneDB), silent = TRUE)
      if(class(geneAnnot)[1] != 'try-error' & !is.null(geneAnnot))
        gPlot <- .addTag(gPlot, geneAnnot, input$Yexpand, input$gain, input$loss)
        }
    print(gPlot)
    })
  
  createSummary <- reactive({
    if(gene$symbol %in% c('NONE', ''))
      return(NULL)
    geneAnnot <- try(.geneOfInt(reCenterSeg(), gene$symbol, geneDB), silent = TRUE)
    if(class(geneAnnot)[1] != 'try-error' & !is.null(geneAnnot)){
      selected <- geneAnnot[,c("symbol", "fullName", "cytoband", "entrezgeneId", "Log2Ratio", "segNum", "segLength(kb)")]
#      colnames(selected)[ncol(selected)] <- "segment.len(kb)"
      selected$Log2Ratio <- round(selected$Log2Ratio, 3)
    } else{
      selected <- data.frame(message=sprintf("\"%s\" does not seem to be an official symbol.", gene$symbol))
    }
    return(selected)
  })

  createFullTable <- reactive({
    gt <- reCenterGenes()
    gt$Log2Ratio <- round(gt$Log2Ratio, 3)
    if(input$chr=="All") chr <- 1:23
    else chr <- as.numeric(input$chr)
    gt <- gt[gt$chr %in% chr,]
    greater <- as.numeric(input$gain)
    lower <- (-abs(as.numeric(input$loss)))
    gt <- gt[gt$Log2Ratio>=greater | gt$Log2Ratio<=lower,]
    if(nrow(gt)>0){
      gt[, "entrezgeneId"] <- .renderLink(gt[, "entrezgeneId"])
      return(gt[,c("symbol", "fullName", "chr", "cytoband", "entrezgeneId", "Log2Ratio", "segNum", "segLength(kb)")])
    }
    else
      return(NULL)
  })
  
  createTitle1 <- reactive({
    return( unique(segTable$ID) )
  })
  createTitle2 <- reactive({
    hi <- input$gain
    lo <- input$loss
    return( sprintf("Gain threshold: %s, Loss threshold: %s", hi, lo) )
  })
  
  output$Profile <- renderPlot({ createCGHplot() }, res=120, width=2000, height=1058)
  output$geneSummary <- renderTable({ createSummary() })
  output$tableTitle1 <- renderText({ createTitle1() })
  output$tableTitle2 <- renderText({ createTitle2() })
  output$fullTable <- renderDataTable({ createFullTable() }, options=list(lengthMenu=c(25, 50, 100)))
})
