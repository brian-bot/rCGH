###########################
###########################
.getName <- function(segTable){
  return( unique(as.character(segTable$ID)) )
}
.selectChr <- function(Choice){
  if(Choice == 'All')
    return(1:23)
  else
    return(as.numeric(Choice))
}
.mainPlot <- function(segTable, s=5){
    #  w <- 16
    if(nrow(segTable)<1)
        return(NULL)
    X <- lapply(1:nrow(segTable), function(i){
        #n <- ceiling(segTable$num.mark[i]/w)
    n <- max(25, segTable$num.mark[i])
    x <- seq(segTable$loc.start[i], segTable$loc.end[i], len=n)
    y <- rnorm(n, segTable$seg.med[i], segTable$probes.Sd[i]/s)
    return(cbind(loc=x, l2r=y))
    })
    X <- as.data.frame(do.call(rbind, X))
  
  gPlot <- ggplot(data=X, aes(x=loc, y=l2r)) +
    geom_point(pch = 19, cex = 0.15, col = 'grey50') +
    geom_hline(yintercept = 0) +
    xlab('Genomic position (bp)') +
    ylab('Log2(Ratio)') +
    theme_bw() +
    theme(  panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text=element_text(size=18),
            axis.title=element_text(size=22),
            axis.title.x=element_text(vjust=-.5),
            axis.title.y=element_text(hjust=.5, vjust = .25),
            plot.title = element_text(lineheight=1.1, size = 25, vjust = 2, face="bold")
    )
  return(gPlot)
}
.updateYscale <- function(gPlot, Ymin, Ymax){
  ymin <- min(gPlot$data$l2r, na.rm=TRUE)*Ymin
  ymax <- max(gPlot$data$l2r, na.rm=TRUE)*Ymax
  gPlot <- gPlot + coord_cartesian(ylim = range(min(-1, ymin)-.25, max(1, ymax)+.25))
  return(gPlot)
}
.updateXscale <- function(gPlot, chr, hg19){
  if(chr=="All")
    return(gPlot)
  else {
    chr <- as.numeric(chr)
    cumLen <- cumsum(as.numeric(hg19$length))
    xmin <- ifelse(chr==1, 0, cumLen[chr-1])
    xmax <- cumLen[chr]
    gPlot <- gPlot + coord_cartesian(xlim = range(xmin, xmax))
    return(gPlot)
  }
}
.addSegments <- function(gPlot, segTable, chr, gain, loss){
  if(chr=="All") chr <- 1:23
  else chr <- as.numeric(chr)
  
  myBlue <- rgb(0, 0.45, 1, 1)
  
  segTable <- segTable[which(segTable$chrom %in% chr),]
  idx <- which(segTable$seg.med<= loss | segTable$seg.med>= gain)
  if(length(idx)>0){
    subTable <- segTable[idx,]
    GLcolors <- ifelse(subTable$seg.med<= loss, 'red3',
                       ifelse(subTable$seg.med>= gain, myBlue, "black"))
    gPlot <- gPlot+
      geom_segment(data=subTable,
                   aes(x=loc.start, xend=loc.end, y=seg.med, yend=seg.med),
                   colour=GLcolors, size=2)
  }
  return(gPlot)   
}
.addChr <- function(gPlot, chr, hg19){
  if(chr=="All") chr <- as.numeric(1:23)
  else chr <- as.numeric(chr)
  
  ylim <- gPlot$coordinates$limits$y
  if(is.null(ylim))
    ylim <- range(gPlot$data$l2r)
  cumCentr <- 1/2*hg19$length+hg19$cumlen
  gPlot <- gPlot+
    geom_vline(xintercept = hg19$cumlen[chr], color = 'grey30', linetype = 2, size = 0.25) +
    annotate('text', x=cumCentr[chr], y=rep(max(ylim, na.rm=TRUE)*.95, length(chr)),
             label=chr, size = 4, colour = 'grey40')
  return(gPlot)     
}
.addTitle <- function(gPlot, sampleName, gain, loss){
  Title = paste(sampleName, '\nGain threshold: ', round(gain, 3), ' Loss threshold:', round(loss, 3))
  gPlot <- gPlot + ggtitle(Title)
  return(gPlot)  
}
.addTag <- function(gPlot, geneAnnot, Yexpand, gain, loss){
  tagExp <- 2.5
  myBlue <- rgb(0, 0.45, 1, 1)
  ylim <- gPlot$coordinates$limits$y
  ymin <- min(ylim); ymax <- max(ylim)
  symbol <- as.character(geneAnnot$symbol)
  lr <- geneAnnot$Log2Ratio
  xLabel <- ifelse((lr+tagExp)<ymax, geneAnnot$genomStart, geneAnnot$genomStart-2e8)
  yLabel <- ifelse((lr+tagExp)<ymax, lr+tagExp, lr-1)
  
  if(is.na(lr))
    return(gPlot)
  Col <- ifelse(lr<= loss, 'red3', ifelse(lr>=gain, myBlue, 'grey40'))
  Loc <- data.frame(xstart=xLabel,
                    xend=geneAnnot$genomStart,
                    ystart=ifelse((lr+tagExp)<ymax, lr+1.75, lr),
                    yend=ifelse((lr+tagExp)<ymax, lr+.25, lr),
                    Col=Col)
  gPlot <- gPlot +
    annotate("text", x=max(xLabel, 2e8), y=yLabel,
             label = paste0(symbol, '\n(Log2R = ', round(lr, 3), ')'), cex = 7, colour=as.character(Loc$Col)) +
    geom_segment(data=Loc, aes(x=xstart, xend=xend, y=ystart, yend=yend),
                 colour = "black", arrow=arrow(angle=30, length=unit(.3, "cm")), size = 1.25)
  return(gPlot)
}
.geneOfInt <- function(segTable, symbol, geneDB){
  symbol <- toupper(symbol)
  tmp <- geneDB[which(geneDB$symbol == symbol),]
  if(nrow(tmp)==0){
    return(NULL)
  }
  if(!tmp$chr %in% segTable$chrom)
    return(cbind(tmp, Log2Ratio = NA, "segNum"=NA, "segLength(kb)"=NA))
  nPrev <- sum(segTable$chrom<tmp$chr)
  segTable <- segTable[segTable$chrom==tmp$chr,]
  geneStart <- tmp$genomStart
  geneEnd <- tmp$genomEnd
  bound1 <- max(which(segTable$loc.start <= geneStart))
  bound2 <- min(which(segTable$loc.end >= geneEnd))
  bounds <- unique(c(bound1, bound2))
  if(is.na(geneStart) | is.na(geneEnd)){
    tmp <- cbind(tmp, Log2Ratio=NA, "segNum"=NA, "segLength(kb)"=NA)
  } else if(length(bounds) > 0){
    geneLR <- unique(segTable$seg.med[bounds])
    tmp <- cbind(tmp, Log2Ratio=geneLR, "segNum"=nPrev+bounds, "segLength(kb)"=abs(segTable$loc.end[bounds] - segTable$loc.start[bounds])/1e3)
  } else{
    tmp <- cbind(tmp, Log2Ratio=NA, "segNum"=NA, "segLength(kb)"=NA)
  }
  tmp$entrezgeneId <- as.factor(tmp$entrezgeneId)
  return( as.data.frame(tmp) )
}
.renderLink <- function(uid){
  sprintf("<a href=\"http://www.ncbi.nlm.nih.gov/gene/?term=%s[uid]\" target=\"_blank\" style=\"font-size:18px; \">%s</a>", uid, uid)
}

# End helper functions
###########################
###########################
