################################################
# HELPER FUNCTIONS (ALL INTERNAL TO THE PACKAGE)
################################################

###############################################
# Loading required files
###############################################
# cat("Loading required files...")
# path <- system.file("inst/extdata", package = "rCGH")
# hg19 <- readRDS(file.path(path, "human.chrom.info.hg19.FC.rds"))
# geneDB <- readRDS(file.path(path, "geneAnnots_2013_10_28.rds"))
# agilentDB <- readRDS(file.path(path, "022060_D_F_20111015.rds"))
# cat("Done.\n)
###############################################
# Preprocessings called in constructors
###############################################
# These internal functions are in buildData.R
#
###############################################
# Preprocessings called in adjustSignal.R
###############################################
.CyAdjust <- function(cnSet, Ref){
  '
	Called by adjustSignal(object)
  If Fract (prop of tumor in the sample), adjust the tumor signal (Cy5)
  Compute and adjust the Log2(Cy3/Cy5)
  Return the data.frame with a supplementary column: Log2Ratio
  '
  cat('Cy effect adjustment...\t')
  if(Ref=="cy3"){
    ref <- log2(cnSet$gMedianSignal)					# Ref in Cy3 (g)
    test <- log2(cnSet$rMedianSignal)					# Test in Cy5 (r)
  } else{
    ref <- log2(cnSet$rMedianSignal)					# Ref in Cy5 (r)
    test <- log2(cnSet$gMedianSignal)					# Test in Cy3 (g)
  }
  M <- test - ref
  A <- (test + ref)/2
  Loess <- loessFit(M, A)$fitted
  LR <- M - Loess
  cnSet$Log2Ratio <- LR
  cnSet$Log2Ratio - median(cnSet$Log2Ratio, na.rm = T)
  cat('Done.\n')
  return (cnSet)
}
.GCadjust <- function(cnSet, gcDB){
  '
  Called by adjustSignal(object)
  Adjust the GC%
  '
  cat('GC% adjustment...')
  cnSet <- cnSet[order(cnSet$ProbeName),]
  gcDB <- gcDB[order(gcDB$ProbeID),]
  
  idx <- match(cnSet$ProbeName, gcDB$ProbeID)
  cnSet <- cnSet[which(!is.na(idx)), ]
  gcDB <- gcDB[idx[!is.na(idx)],]
  
  if(!all(as.character(cnSet$ProbeName) == as.character(gcDB$ProbeID)))
    stop("An error occured in GCadjust: probeIds do not match with Agilent grid.\n")
  
  lr = cnSet$Log2Ratio
  GC <- gcDB$GCpercent
  adjLr <- lr - loessFit(lr, GC)$fitted
  cnSet$Log2Ratio <- adjLr
  cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
  cat('\tDone.\n')
  
  return(cnSet)
}
.dlrs <- function(x){
  nx <- length(x)
  if (nx<3) {
    stop("Vector length>2 needed for computation")
  }
  tmp <- embed(x,2)
  diffs <- tmp[,2]-tmp[,1]
  dlrs <- IQR(diffs, na.rm = TRUE)/(sqrt(2)*1.34)
  return(dlrs)
}
.MAD <- function(LR){
  tmp <- abs(LR - median(LR, na.rm = TRUE))
  return(median(tmp, na.rm = TRUE))
}
.supprOutliers <- function(x, n){
  # x: vector of values
  # n: window width
  '
  Called by adjustSignal(object)
  Replace the outliers by the medians of the +/- n neighbours
  '
  cat('Suppressing outliers...')
  nx <- length(x)
  tmp <- lapply(1:(nx-2*n), function(i) x[i:(i+2*n)])
  Q <- sapply(tmp, function(x) IQR(x[-(n+1)], na.rm = T))
  M <- sapply(tmp, function(x) median(x[-(n+1)], na.rm = T))
  tmp <- do.call(rbind, tmp)
  idx <- which(tmp[,(n+1)] < M - 1.5*Q | tmp[,(n+1)] > M + 1.5*Q)
  tmp[idx,(n+1)] <- M[idx]
  cat(length(idx), 'probes adjusted.\n')
  newx <- c(x[1:n], tmp[,(n+1)], x[(nx-n+1):nx])
  return(newx)
}
.preSeg <- function(cnSet, sampleName, params){
  cat('Computing pre-segmentation...\n')
  cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
  params$UndoSD <- .1
  
  L2R <- cnSet$Log2Ratio
  Chr <- cnSet$ChrNum
  Pos <- cnSet$ChrStart
  if(is.na(sampleName) | is.null(sampleName))
    sampleName <- "sample_x"
  
  segTable <- .computeSegmentation(L2R, Chr, Pos, sampleName, params, Smooth=11)
  segTable <- .smoothSeg(segTable, minMarks=3)
  segTable <- .computeMedSegm(segTable, L2R)    # helper function
  segTable <- .mergeLevels(segTable, thresMin=.1)
  
  cnSet <- .newLog(cnSet, segTable)
  return(cnSet)
}
.newLog <- function(cn, st){
  for(idx in 1:nrow(st)){
    ii <- which(cn$ChrNum==st$chrom[idx] & cn$ChrStart==st$loc.start[idx])
    jj <- which(cn$ChrNum==st$chrom[idx] & cn$ChrStart==st$loc.end[idx])
    cn$Log2Ratio[ii:jj] <- rnorm(length(ii:jj), st$seg.med[idx], st$probes.Sd[idx])
  }
  return(cn)
}

###############################################
# Processings called in EMnormalize.R
###############################################
.smoothLR <- function(LR, Platform, cut, K){
  '
	Called by EMnormalize(object)
  Smoothing the LR vector to improve the EM classification.
  '
  if(any(is.na(LR)))
    LR <- LR[!is.na(LR)]
  if(grepl("Affymetrix", Platform)){
    ii <- seq(1, length(LR), by = 6)
    LR <- LR[ii]
  } else {
    ii <- seq(1, length(LR))
  }
  runLR <- runmed(LR, k = K)	
  q1 <- cut[1]; q2 <- cut[2]
  index <- which(runLR>=q1 & runLR<=q2)
  return (list(runLR=runLR[index], index=ii[index]))
}
.buildEMmodel <- function(LR, G, Len){
  '
	Called by EMnormalize(object)
  Model the distribution as a gaussian mixture.
  '
	if(is.na(Len))
    Len <- floor(length(LR)*.25)
	model <- Mclust(LR[sample(1:length(LR), size = Len)], G = G)
  nG <- model$G
	p <- model$parameters$pro
	m <- model$parameters$mean
	s <- model$parameters$variance$sigmasq
	if(length(s)<length(m)) s <- rep(s, length(m))
	p <- p[order(m)]
	s <- s[order(m)]
	m <- m[order(m)]
	return(list(nG = nG, m = m, p = p, s =s))
}
.mergePeaks <- function(nG, n, m, s, p, mergeVal){
  cat("Merging peaks closer than", mergeVal, "...\n")
  Raw <- c(1, 1)
  while(length(Raw)!=0){
    Mdist <- matrix(0, nG, nG)
    for(i in 1:nG)
      for(j in 1:nG){
        Mdist[i, j] <- abs(m[i] - m[j])
      }
    diag(Mdist) <- NA
    Raw <- ceiling(which(Mdist<mergeVal)/nG)
    
    if(length(Raw)!=0){
      C1 <- Raw[1]
      C2 <- Raw[2]
      mu1 <- m[C1]; mu2 <- m[C2]
      s1 <- s[C1]; s2 <- s[C2]
      p1 <- p[C1]; p2 <- p[C2]
      newmu <- (p1*mu1 + p2*mu2)/(p1 + p2)
      news <- (p1*(s1+mu1^2) + p2*(s2+mu2^2))/(p1 + p2) - newmu^2
      #        news <- ((p1*n - 1)*s1 + (p2*n - 1)*s2)/(p1*n + p2*n - 2)
      newp <- p1 + p2
      m[C1] <- newmu; s[C1] <- news; p[C1] <- newp
      m <- m[-C2]; s <- s[-C2]; p <- p[-C2]
      nG <- length(m)
      cat("Merged parameters:\n")
      cat("means:", m, "\nVar:", s, "\nprops:", p, "\n\n")
    }
  }
  return(list(nG = nG, m = m, s = s, p = p))
}
.computeDensities <- function(n, m, p, s){
  '
	Called in EMnormalize(object)
  Simulates the mixture model according to the returned EM paramaters.
  '
	dList = list()
	peaks <- c()
	for(i in 1:length(m)){
		tmp <- rnorm(n*p[i], m[i], sqrt(s[i]))
		tmpD <- density(tmp, na.rm = T)
		tmpD$y = tmpD$y *p[i]
		dList[[i]] <- tmpD
		peaks <- c(peaks, max(tmpD$y))
		}
	return(list(dList = dList, peaks = peaks))
}
.chooseBestPeak <- function(peaks, m, peakThresh){
	best <- which(peaks>=max(peaks)*peakThresh)
  if(length(best)>0){
	  cat(length(best), 'peak(s) above', sprintf("%s%s", peakThresh*100, "%"), 'of max peak.\n')
  }
  else{
    cat('No peak above threshold:', peakThresh, '\n')
  }
  bestPeak <- best[1]
	cat('Peak at', m[bestPeak], 'has been chosen.\n')
	return(bestPeak)
}
.plotEMmodel <- function(LR, dList, m, bestPeak, cut, Title){
  dLR <- density(LR)
  currentPlot = xyplot(dLR$y ~ dLR$x, type = "n",
                       xlab=list(
                         label=expression(Log[2](ratio)),
                         cex=1.5),
                       xlim = cut,
                       ylim = range(0, max(dLR$y)*1.5),
                       ylab=list(
                         label="Density",
                         cex=1.5),
                       scales=list(x=list(cex=1.25), y=list(cex=1.25)),
                       panel = function(x, y){
                         lpolygon(dLR$x, dLR$y, col = 'grey90')
                         n = length(LR)
                         nG = length(dList)
                         for (i in 1:nG){
                           tmp <- dList[[i]]
                           ymin <- min(max(tmp$y*1.25, na.rm = TRUE), max(dLR$y)*1.25)
                           llines(tmp$x, tmp$y, lwd = 1, col = rgb(i/nG, 0.2, (nG-i)/nG, 0.75))
                           lpolygon(tmp$x, tmp$y, col = rgb(i/nG, 0.2, (nG-i)/nG, 0.25))
                           ltext( x = mean(tmp$x), y = max(0.25, ymin),
                                  labels = round(m[i], 3), cex = ifelse(i == bestPeak, 1.5, 1.25),
                                  font =  ifelse(i == bestPeak, 2, 1))
                         }
                       }
  )
  return(currentPlot)
}

###############################################
# Processings called in segmentCGH.R
###############################################
.computeSegmentation <- function(LR, Chr, Pos, sampleId, params, Smooth){

  Ksmooth <- params$Ksmooth
  Kmax <- params$Kmax
  Nmin <- params$Nmin
  Mwidth <- params$Mwidth
  Alpha <- params$Alpha
  UndoSD <- params$UndoSD
  
  cna.obj <- CNA(LR, Chr, Pos, data.type = "logratio", sampleid = sampleId)
  
  if(Smooth){
    smooth.cna.obj <- smooth.CNA(cna.obj, smooth.region = Ksmooth)
    seg.cna.obj <- segment(smooth.cna.obj, undo.splits = "sdundo", undo.SD = UndoSD, alpha = Alpha, kmax = Kmax, nmin = Nmin, min.width = Mwidth)
  }
  else seg.cna.obj <- segment(cna.obj, undo.splits = "sdundo", undo.SD = UndoSD, alpha = Alpha, kmax = Kmax, nmin = Nmin, min.width = Mwidth)
  cat('UndoSD:', UndoSD, '\n')
  return(seg.cna.obj$output)
}
.smoothSeg <- function(segTable, minMarks){
  splitSegTables <- split(segTable, segTable$chrom)
  
  cat("Removing segments shorter than", minMarks, "markers.\n")
  adjustedLocs <- lapply(splitSegTables, function(sst){
    while(any(sst$num.mark<minMarks)){
      i <- which(sst$num.mark<minMarks)[1]
      j <- .getCloser(sst, i)
      sst <- .mergeSegments(sst, i, j)
      sst <- sst[-i,]
    }
    return(sst)
  }
  )
  adjustedLocs <- as.data.frame(do.call(rbind, adjustedLocs))
  rownames(adjustedLocs) <- seq(1, nrow(adjustedLocs))
  return(adjustedLocs)
}
.getCloser <- function(segTable, idx){
  if(idx==1)
    return(idx+1)
  else if (idx==nrow(segTable))
    return(idx-1)
  else{
    delta <- abs(segTable$seg.mean[c(idx-1, idx+1)] - segTable$seg.mean[idx])
    i <- ifelse(which.min(delta)==1, idx-1, idx+1)
    return(i)
  }
}
.mergeSegments <- function(segTable, i, j){
  if(j<i){
    segTable$loc.end[j] <- segTable$loc.end[i]
  } else {
    segTable$loc.start[j] <- segTable$loc.start[i]
  }
  segTable$num.mark[j] <- segTable$num.mark[j] + segTable$num.mark[i]
  return(segTable)
}
.computeMedSegm <- function(segTable, L2R){
  '
	Called by SegmentCGH()
  '
  nMark <- segTable$num.mark
  e = 0
  seg.med <- Sd <- c()
  for(i in 1:length(nMark)){
    s = e + 1
    e = e + nMark[i]
    tmpL2R <- L2R[s:e]
    tmpMed <- tukey.biweight(tmpL2R[!is.na(tmpL2R)])
    tmpSd <- sd(tmpL2R[!is.na(tmpL2R)], na.rm=TRUE)
    seg.med <- c(seg.med, tmpMed)
    Sd <- c(Sd, tmpSd)
  }
  segTable <- cbind.data.frame(segTable, seg.med = seg.med, probes.Sd = Sd)
  return(segTable)
}
.mergeLevels <- function(st, thresMin=0.1, ...){
  vObs <- st$seg.mean
  vPred <- st$seg.med
  vFit <- mergeLevels(vObs, vPred, thresMin=thresMin, ...)
  st$seg.med <- vFit$vecMerged
  return(st)
}
.probeSegValue <- function(segTable, use.medians){  
  '
	Called by SegmentCGH()
	'
  segVal <- segTable$seg.mean
  if(use.medians) segVal <- segTable$seg.med
  nMark <- segTable$num.mark
  output <- lapply(1:length(segVal), function(i){rep(segVal[i], nMark[i])})
  return(do.call(c, output))
}

###############################################
# Processings called in byGeneTable.R
###############################################
.ByGene <- function(segTable, geneDB, hg19){
  cat("Creating byGene table...")
  genelist <- lapply(1:nrow(segTable), function(i){
    .getGeneList(i, segTable[i,], geneDB)
  })
  sampleId <- unique(segTable$ID)
  cmValues <- .cmValues(segTable, hg19)
  bygene <- do.call(rbind, genelist)
  bygene$relativeLog <- .relativeLog(bygene, cmValues, hg19)
  genelist <- cbind.data.frame(patientId=.getPatientId(sampleId), sampleId=sampleId, bygene)
  cat("Done\n")
  return(genelist)
}
.getGeneList <- function(i, seg, geneDB, removeConflicts=TRUE){
  selecItems <- c("symbol", "fullName", "chr", "cytoband", "chrStart", "genomStart", "entrezgeneId")
  segLen <- abs(seg$loc.end - seg$loc.start)
  tmpDB <- geneDB[geneDB$chr==seg$chrom,]
  idx <- which(seg$loc.start<=tmpDB$chrStart & tmpDB$chrStart<=seg$loc.end |
                 seg$loc.start<=tmpDB$chrEnd & tmpDB$chrEnd<=seg$loc.end)
  if(length(idx)>0){
    bygene <- cbind.data.frame(tmpDB[idx, selecItems],
                               Log2Ratio=seg$seg.med,
                               num.mark=seg$num.mark,
                               segNum=i,
                               "segLength(kb)"=segLen/1e3)
    if(removeConflicts)
      bygene <- .conflicts(bygene)
    
    return(bygene)
  }
}
.cmValues <- function(segTable, hg19){
  cmLocs <- .locateCM(segTable, hg19)
  lapply(cmLocs, function(locs) c(segTable$seg.med[locs[1]], segTable$seg.med[locs[2]]))
}
.locateCM <- function(segTable, hg19){
  chrs <- unique(segTable$chrom)
  cmLocs <- lapply(chrs, function(chr){
    cStart <- hg19$centromerStart[hg19$chrom==chr]
    cEnd <- hg19$centromerEnd[hg19$chrom==chr]
    tmp <- segTable[segTable$chrom==chr,]
    locStart <- which(tmp$loc.start<=cStart & cStart<=tmp$loc.end)
    locEnd <- which(tmp$loc.start<=cEnd & cEnd<=tmp$loc.end)
    if(length(locStart)==0)
      locStart <- which.min(abs(tmp$loc.end - cStart))
    if(length(locEnd)==0)
      locEnd <- which.min(abs(tmp$loc.start - cEnd))
    as.numeric(rownames(tmp)[c(locStart, locEnd)])
  })
  return(cmLocs)
}
.relativeLog <- function(bygene, cmValues, hg19){
  relativeLog <- lapply(1:length(cmValues), function(chr){
    tmp <- bygene[bygene$chr==chr, ]
    ii <- which(tmp$chrStart<hg19$centromerStart[hg19$chrom==chr])
    jj <- which(tmp$chrStart>hg19$centromerEnd[hg19$chrom==chr])
    rl <- rep(NA, nrow(tmp))
    rl[ii] <- tmp$Log2Ratio[ii] - cmValues[[chr]][1]
    rl[jj] <- tmp$Log2Ratio[jj] - cmValues[[chr]][2]
    return(rl)
  })
  return(do.call(c, relativeLog))
}
.conflicts <- function(bygene){
  dup <- intersect(which(duplicated(bygene$symbol)), which(duplicated(bygene$Log2Ratio)))
  if(length(dup)>0)
    bygene <- bygene[-dup,]
  return(bygene)
}
.getPatientId <- function(sampleId){
  gsub("(.*)_(.*)_(.*)+", "\\2", sampleId)
}

###############################################
# Processings called in view.R
###############################################

.convertLoc <- function(segTable, hg){
  ss <- split(segTable, segTable$chrom)
  sconv <- lapply(ss, function(tmp){
    chr <- unique(tmp$chrom)
    tmp$loc.start <- tmp$loc.start + hg$cumlen[chr]
    tmp$loc.end <- tmp$loc.end + hg$cumlen[chr]
    return(tmp)
  })
  return(as.data.frame(do.call(rbind, sconv)))
}

###############################################
# NOT USED
###############################################

# 
# .chrCumLen <- function(hg){
#   cumLen = cumsum(as.numeric(hg$length))
#   cumLen = c(0, cumLen[-length(cumLen)])
#   return(cumLen)
# }
# 
# .AdjustCent <- function(eset, N, K, thresh){
#   # N: how many centromer probes to consider
#   .getProbes <- function(Values, Loc, D, N){
#     if(length(D) >= N) return(Values[D[1:N]])
#   }
#   if(!exists('hg19')){
#     ent <- synGet('syn2141399')
#     hg19 <- read.csv(ent@filePath, header = TRUE, sep = '\t')
#   }
#   Cmers <- hg19[hg19$chrom == unique(eset$ChrNum),]
#   probesLoc = eset$ChrStart
#   L2R = eset$Log2Ratio
#   if(any(is.na(L2R))) L2R[is.na(L2R)] <- median(L2R, na.rm = TRUE)
#   chr <- eset$ChrNum
#   Rmed <- runmed(L2R, k = K)
#   DL <- order(which(probesLoc < Cmers$centromerStart), decreasing = TRUE)
#   DR <- which(probesLoc > Cmers$centromerEnd)
#   probesLeft <- probesRight <- rnorm(N, tukey.biweight(Rmed), sd(Rmed))
#   if(length(DL) >= N)
#     probesLeft <- .getProbes(Rmed, probesLoc, DL, N)
#   if(length(DR) >= N)
#     probesRight <- .getProbes(Rmed, probesLoc, DR, N)
#   p <- t.test(probesLeft, probesRight)$p.value
#   test <- cbind.data.frame(chr = Cmers$chrom,
#                            leftMean = mean(probesLeft, na.rm = TRUE),
#                            rightMean = mean(probesRight, na.rm = TRUE),
#                            pvalue = p)
#   if(p > thresh){
#     allProbes <- c(probesLeft, probesRight)
#     L2R <- L2R - mean(allProbes, trim = .05)
#   }
#   test <- cbind.data.frame(test, adjusted = ifelse(test$pvalue > thresh, "yes", 'no'))
#   return(list(L2R = L2R, test = test))
# }
# 
# 
# 
# dLRsd <- function(LR){
# 	'
# 	Not used yet
# 	'
# 	n <- length(LR)
# 	V1 <- LR[-1]
# 	V2 <- LR[-n]
# 	dLR <- V2-V1
# 	q1 <- quantile(dLR, 0.25, na.rm = TRUE)
# 	q3 <- quantile(dLR, 0.75, na.rm = TRUE)
# 	s <- sd(dLR[which(dLR > q1 & dLR < q3)], na.rm = TRUE)/sqrt(2)
# 	return(s)
# 	}

#########################
# Gene tables functions
#########################


# .addGenomicPos <- function(cnSet){
#   cat('Adding hg19 genomic positions...\t')
#   if(!exists('hg19')){
#     hg19 <- readRDS(synGet('syn2342423')@filePath)
#   }
#   cumLen <- hg19$cumlen
#   cnSet = cnSet[order(cnSet$ChrNum, cnSet$ChrStart, cnSet$ProbeName), ]
#   chrTable = table(cnSet$ChrNum, useNA = 'ifany')
#   genomic = c()
#   for(i in 1:length(chrTable)){
#     n = chrTable[i]
#     genomic = c(genomic, rep(cumLen[i], n))
#   }
#   genomic = ifelse(!is.na(cnSet$ChrStart), genomic + cnSet$ChrStart, NA)
#   cnSet <- cbind.data.frame(cnSet[,c("ProbeName", "ChrNum", "ChrStart", "ChrEnd")],
#                             genomicPos = genomic,
#                             Log2Ratio = cnSet[,'Log2Ratio'],
#                             stringsAsFactors=FALSE)
#   .checkOverlaps(cnSet)
#   return(cnSet)
# }
# .checkOverlaps <- function(cnSet){
#   overlaps <- sapply(seq(2, 23), function(chr){
#     last <- max(cnSet$genomicPos[cnSet$ChrNum == chr-1], na.rm = TRUE)
#     first <- min(cnSet$genomicPos[cnSet$ChrNum == chr], na.rm = TRUE)
#     last > first
#   })
#   if(any(overlaps)) cat('overlap between', which(overlaps)[-1], "and", which(overlaps), '\n')
#   else cat('No overlap\n')
# }

###############################################
# QCs
###############################################



# .Iqr <- function(x, n){
#   '
#   Called by .supprOutliers
#   If xn is an outlier, then replace xn by the median of the +/-n neighbours (computed without xn itself)
#   '
#   xprim = x[-(n+1)]
#   q = quantile(xprim, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
#   iqr = q[3] - q[1]
#   if(is.na(x[n+1]) | (x[n+1] < q[2]-1.5*iqr | x[n+1] > q[2]+1.5*iqr))
#     x[n+1] <- q[2]
#   return(x)
# }


###############################################





#################
# To Do: QCSegm
#################

# .getSD <- function(object, n=1000, q=.95, B=10000){
#   cn <- getCNset(object)
#   x <- cn$Log2Ratio
#   N <- length(x)
#   SD <- sapply(sample(N, B, replace = TRUE), function(i){
#     if(i>N-n+1)
#       return(NA)
#     return(sd(x[i:(i+n)]))
#   })
#   quantile(SD, q, na.rm=TRUE)
# }
# .armSeg <- function(arm, ...){
#   if(nrow(arm)==0)
#     return(NULL)
#   lrr <- as.vector(as.vector(arm$Log2Ratio))
#   cp <- try(cpt.meanvar(lrr, method="BinSeg", penalty="Asymptotic", pen.value = 1e-12), silent=TRUE)
#   if(class(cp)[1]=="try-error")
#     cp <- cpt.meanvar(lrr, method="BinSeg", penalty="SIC", pen.value = 1e-12)
#   cpts <- cp@cpts
#   idx <- c(1, cpts[-length(cpts)])
#   seg.med <- mapply(function(ii, jj) median(lrr[ii:jj], na.rm=TRUE), idx, cpts)
#   segs <- data.frame(chrom=rep(unique(arm$ChrNum), length(idx)),
#                      loc.start=arm$ChrStart[idx],
#                      loc.end=arm$ChrStart[cpts-1],
#                      num.mark=cpts-1-idx,
#                      seg.mean=cp@param.est$mean,
#                      seg.med=seg.med,
#                      probes.Sd=sqrt(cp@param.est$variance)
#   )
#   return(segs)
# }
# .binSeg <- function(cnSet, hg19, minMarks=10, ...){
#   Chrs <- unique(cnSet$ChrNum)
#   segs <- lapply(Chrs, function(chr){
#     cat("Segmenting", chr, "...")
#     tmp <- cnSet[cnSet$ChrNum==chr,]
#     parm <- tmp[tmp$ChrStart<=hg19$centromerStart[chr],]
#     qarm <- tmp[tmp$ChrStart>=hg19$centromerEnd[chr],]
#     seg <- rbind(.armSeg(parm,...), .armSeg(qarm,...))
#     cat(nrow(seg), "segments.\n")
#     return(seg)
#   })
#   st <- as.data.frame(do.call(rbind, segs))
#   st <- .smoothSeg(st, minMarks)
#   st <- .mergeLevels(st, ...)
#   return(st)
# }
