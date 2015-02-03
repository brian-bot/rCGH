
#########################
# SET FUNCTION
#########################
setMethod('setInfo', signature('cghObj'), function(object, item = NULL, value = NULL) {
  if (is.null(item)){
    cat("No item specified.\n")
  } 
  else if (is.null(value)){
    cat("No value specified.\n")
  } else 
    object@info[item] <- value
  return(object)
}
)

#########################
# ACCESSOR FUNCTIONS
#########################
setMethod('getInfo', signature('cghObj'), function(object, item = NULL) {
  info <- object@info
  if (is.null(item)){
    return (info)
  } else if (item %in% names(info)){
    as.character(info[item])
  } else
    cat(paste('\'', item, '\'', sep = ''), 'is not an available item.\n')
}
)
setMethod('getByGene', 'cghObj', function(object) object@byGene)
setMethod('getCNset', 'cghObj', function(object) object@cnSet)
setMethod('getDensity', 'cghObj', function(object, ...) {
  if(is.null(Title)){
    correct <- getParam(object)$correctionValue
    Title <- sprintf("%s\nEM centralization: correction.factor = %s", getInfo(object, 'sampleName'), round(correct, 5))
  }
  dplot <- object@probesDensity
  update(dplot, ylim=dplot$y.limits*.9, main=list(label=Title, cex=2), ...)
})
setMethod('getParam', 'cghObj', function(object) object@param)
setMethod('getProfile', 'cghObj',
          function(object, gene=NULL, gain=.5, loss=(-.5), Title=NULL, s=5, hg=hg19,...) {
            
            myBlue <- rgb(0, 0.45, 1, 1)
            
            segTable <- getSegTable(object)
            segTable <- segTable[which(segTable$chrom != 24),]
            hg <- hg[1:23,]
            segTable <- .convertLoc(segTable, hg)
            cumLen <- hg$cumlen
            
            args <- list(...)
            yexist <- "ylim" %in% names(args)
            if(!yexist){
              miny <- min(-2, min(segTable$seg.med, na.rm=TRUE)) - .25
              maxy <- max(2, max(segTable$seg.med, na.rm=TRUE)) + .25
              ylim <- range(miny, maxy)
            } else{
              ylim <- range(min(args[["ylim"]]), max(args[["ylim"]]))
            }
            cat("ylim", ylim, "\n")
            
            idx <- which(segTable$seg.med<= loss | segTable$seg.med>= gain)
            subTable <- segTable[idx,]
            
            GLcolors <- ifelse(subTable$seg.med<= loss, 'red3',
                               ifelse(subTable$seg.med>= gain, myBlue, NA))
            
            if(is.null(Title))
              Title = paste(getInfo(object, 'sampleName'), '-',
                            getInfo(object, 'analyseDate'),
                            '\nGain threshold: ', round(gain, 3), ' Loss threshold:', round(loss, 3))
            platform <- getInfo(object, "platform")
            
            w <- ifelse(grepl("Affymetrix", platform), 90, 10)
            X <- lapply(1:nrow(segTable), function(i){
              n <- ceiling(segTable$num.mark[i]/w)
              n <- max(15, n)
              x <- seq(segTable$loc.start[i], segTable$loc.end[i], len=n)
              y <- rnorm(n, segTable$seg.med[i], segTable$probes.Sd[i]/s)
              return(cbind(chr=segTable$chrom[i], loc=x, l2r=y))
            })
            X <- as.data.frame(do.call(rbind, X))
            
            cumCentr <- 1/2*hg$length + cumLen
            
            gPlot <- ggplot(data = X, aes(x=loc, y=l2r)) +
              geom_point(pch = 19, cex = 0.1, col = 'grey50') +
              geom_hline(yintercept = 0) +
              geom_vline(xintercept = cumLen[1:23], color = 'grey30', linetype = 2, size = 0.25) +
              ggtitle(Title) +
              xlab('Genomic position (bp)') +
              ylab('Log2(Ratio)') +
              theme_bw() +
              theme(  panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      plot.title = element_text(lineheight=.8, size = rel(2.0), face="bold"),
                      axis.title.x = element_text(size = rel(1.8), angle = 00),
                      axis.text.x = element_text(size = rel(1.5)),
                      axis.title.y = element_text(size = rel(1.8), angle = 90),
                      axis.text.y = element_text(size = rel(1.5))
              )
            
            gPlot <- gPlot+
              geom_segment(data = subTable,
                           aes(x = loc.start, xend = loc.end, y = seg.med, yend = seg.med),
                           colour = GLcolors, size = 2)+
              coord_cartesian(ylim = ylim)+
              annotate('text', x = c(-1e8, cumCentr[1:23]), y = rep(max(ylim)-0.2, 24),
                       label = c("Chr", seq(1, 23)), size = 4, colour = 'grey40')
            
            if(!is.null(gene)){
              gene <- toupper(gene)
              bg <- getByGene(object)
              v <- bg[grep(sprintf("^%s$", gene), bg$symbol),]
              if(nrow(v)>0){
                gPlot <- gPlot +
                  geom_point(x=v$genomStart, y=v$Log2Ratio, cex=6, pch=1) +
                  annotate('text', x = v$genomStart - 2e8, y = v$Log2Ratio + .15,
                           label = sprintf("%s\n%s", gene, format(v$Log2Ratio, digits=3)),
                           size = 7, colour = 'grey25')
              }
            }
            
            return(gPlot)
          }
)
setMethod('getSegTable', 'cghObj', function(object) object@segTable)
setMethod('view', 'cghObj',
          function(object, hg=hg19,...) {
            # path <- "~/Documents/myProjects/cgh_workflow_paper/rCGH_backup/inst/shinyProfile"
            #  path <- file.path(.libPaths()[1], "cghViewer")
            path <- system.file("inst", "shinyProfile", package="rCGH")
            platform <- getInfo(object, "platform")
            if(grepl("Affymetrix", platform)){
              w <- 60
            } else{
              w <- 10
            }
            segTable <- .convertLoc(getSegTable(object), hg)
            segTable$num.mark <- round(segTable$num.mark/w)
            saveRDS(segTable, file.path(path, "data/st.rds"))
            saveRDS(getByGene(object), file.path(path, "data/bg.rds"))
            runApp(path, ...)
          })


#########################
# MAIN FUNCTIONS
#########################
adjustSignal <- function(object, Scale=TRUE, Cy=TRUE, Ref="cy3", GC=TRUE, gcDB=agilentDB, preSeg=TRUE, suppOutliers=FALSE){
  
  cnSet <- getCNset(object)
  if (grepl("Agilent", getInfo(object, 'platform'))){
    if(Cy){
      cnSet <- .CyAdjust(cnSet, Ref)			# helper function
      object@param$CyAdjusted = TRUE
    }
    if(GC){
      cnSet <- .GCadjust(cnSet, gcDB)						# helper function
      object@param$GCAdjusted = TRUE
    }
  } else{
    object@param$CyAdjusted = FALSE
    object@param$GCAdjusted = FALSE
  }
  
  object@param$dLRs <- .dlrs(cnSet$Log2Ratio)
  object@param$MAD <- .MAD(cnSet$Log2Ratio)
  if(Scale){
    cat("Scaling...")
    cnSet$Log2Ratio <- scale(cnSet$Log2Ratio, center=FALSE)
    if (grepl("Affymetrix", getInfo(object, 'platform')))
      cnSet$Log2Ratio <- cnSet$Log2Ratio*1.25
    cat("Done.\n")
  }
  if(suppOutliers)
    cnSet$Log2Ratio <- .supprOutliers(cnSet$Log2Ratio,
                                      n = ifelse(grepl("Agilent", getInfo(object, "platform")), 5, 25))
  if(preSeg){
    cnSet <- .preSeg(cnSet, getInfo(object, "sampleName"), getParam(object))
    }
  object@cnSet <- cnSet
  cat('dLRs:', object@param$dLRs, '\tMAD:', object@param$MAD, '\n')
  return(object)
}
EMnormalize <- function(object, cut=c(.1, .9), G=3:6, B=100, peakThresh=0.5, ksmooth=801, useN=2e3, mergeVal=0, Title=NA){
            
            # object:   			an object of class cghObj
            # cut:					Quantiles thresholds in EM centralization
            # G:						Number of groups to consider in EM centralization
            # peakThresh: 		Proportion of the maximum density value to consider a peak as the central reference
            # MergePeaks:		Allow to merge two peaks if there distance is lower than Mergeval
            # MergeVal:			Minimum distance to consider two peaks as different.
            # method:			Define the method to choose the reference peak. Left = left major peak, Zero = major peak closed to Zero, Right = right major peak.
            # Plot:	 				Edit a density plot with EM peaks
            # Expand:				A graphic parameter to expand the x axis
            # Save:				To save automatically the density plot. The current default folder is Root/~/Safir CentralizeProfiles
            # Root: 				The root to save the density plot.
            
            '
            Use a EM algorithm to model the LogRatio density as a mixture of gaussian distributions
            '
            cnSet <- getCNset(object)
            testLR <- cnSet$Log2Ratio
            cut <- as.numeric(quantile(testLR, c(cut[1], cut[2])))
            runLR <- .smoothLR(testLR, getInfo(object, 'platform'), cut, K=ksmooth)		# helper function
            CHR <- cnSet$ChrNum[runLR$index]
            runLR <- runLR$runLR
            
            cat("Analyzing mixture:")
            EMmodels <- lapply(1:B, function(b){cat("."); .buildEMmodel(runLR, G, Len = useN)})	
            nG <- do.call(c, lapply(EMmodels, function(m) m$nG))
            tnG <- table(nG)
            mednG <- as.numeric(names(tnG)[which.max(tnG)])
            
            models <- EMmodels[nG==mednG]
            m <- do.call(rbind, lapply(models, function(m) m$m))
            p <- do.call(rbind, lapply(models, function(m) m$p))
            s <- do.call(rbind, lapply(models, function(m) m$s))
            cat('Done.\n')
            
            m <- apply(m, 2, median, na.rm=TRUE)
            s <- apply(s, 2, median, na.rm=TRUE)
            p <- apply(p, 2, median, na.rm=TRUE)
            
            if(mergeVal>0){
              mergedPars <- .mergePeaks(mednG, length(runLR), m, s, p, mergeVal)
              m <- mergedPars$m
              s <- mergedPars$s
              p <- mergedPars$p
            }

            cat('Computing densities...')
            computeD <- .computeDensities(length(runLR), m, p, s)								# helper function
            dList <- computeD$dList
            peaks <- computeD$peaks
            cat('Done.\n')
            
            cat('Gaussian mixture:\n')
            cat("n.peaks = ", mednG, '\n') 
            bestPeak <- .chooseBestPeak(peaks, m, peakThresh)
            correct <- m[bestPeak]
            
            cat ('\nGroup parameters:\n')
            for (grp in 1:length(m)){
              cat('Grp', grp, '\nprop:', p[grp], ':\tmean:', m[grp],
                  '\tSd:', sqrt(s[grp]), '\tCV:', signif(sqrt(s[grp])/abs(m[grp]), 4), "\n")
            }
            cat("\nEM correction factor = ", correct, "\n")
            
            if(is.na(Title))
              Title <- sprintf("%s\nEM centralization: correction.factor = %s", getInfo(object, 'sampleName'), round(correct, 5))
            EMplot <- .plotEMmodel(runLR, dList, m, bestPeak, cut, Title)
            
            cnSet$Log2Ratio = cnSet$Log2Ratio - correct
            object@cnSet = cnSet
            object@param$EMcentralized = TRUE
            object@param$nPeak = median(nG)
            object@param$PeakProp = as.numeric(p)
            object@param$PeakValues = as.numeric(m)
            object@param$SigmaSq = s
            object@param$centralPeak = as.numeric(bestPeak)
            object@param$correctionValue = as.numeric(correct)
            object@probesDensity = EMplot
            cat('Use getDensity(cghObject) to visualize the probes density plot.\n')
            return(object)
          }
SegmentCGH <- function(object, Smooth=TRUE, UndoSD=NULL, minMarks=8){
            cnSet <- getCNset(object)
            cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
            
            params = getParam(object)
            if(is.null(UndoSD))
              params$UndoSD <- .25 + sqrt(max(getParam(object3)$SigmaSq))
            else params$UndoSD <- UndoSD
            
            L2R <- cnSet$Log2Ratio
            Chr <- cnSet$ChrNum
            Pos <- cnSet$ChrStart
            sampleName <- getInfo(object, 'sampleName')
            if(is.na(sampleName))
              sampleName <- "sample_x"
            
            cat('Computing standard segmentation...\n')
            segTable <- .computeSegmentation(L2R, Chr, Pos, sampleName, params, Smooth)
            segTable <- .smoothSeg(segTable, minMarks)
            segTable <- .computeMedSegm(segTable, L2R)    # helper function
            segTable <- .mergeLevels(segTable)
            probeValues.left <- .probeSegValue(segTable, use.medians = TRUE) #cnSet,
            cat("Number of segments:", nrow(segTable), "\n")
            
            object@param <- params
            object@segTable = segTable
            object@cnSet = cbind.data.frame(cnSet, Segm = probeValues.left)
            
            return(object)
          }
byGeneTable <- function(object, DB=geneDB, HG=hg19){
            segTable <- getSegTable(object)
            object@byGene <- .ByGene(segTable, DB, HG)
            return(object)
          }
