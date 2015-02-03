
################################
# Build a Agilent object
################################

buildAgilent <- function(filePath, Ref="cy3", sampleName=NA, labName=NA, supFlags=TRUE){
  # Load cghData from a synapse entity and build an AgilentObject.
  fileName <- gsub("(.*)/", "", filePath)
  object <- new("cghObj", info = c(fileName=fileName,
                                sampleName=sampleName,
                                labName=labName,
                                platform='Agilent'
                                )
                       )

  object@info <- c(object@info, .readAgilentInfo(filePath))
  object@cnSet <- .readAgilentMatrix(filePath)
  if(supFlags)
    object <- .suppressFlags(object)
  object <- .suppressDuplic(object)
  object <- .preset(object)
  return (object)
}
.readAgilentInfo <- function(filePath){
  cat('Reading information...')
  arrayInfo <- read.csv(filePath, header = F, fill = T, skip = 1, nrows = 8, sep = "\t")
  tmpNames <- as.vector(arrayInfo[1,])
  barCode = as.character(arrayInfo[2, which(tmpNames == "FeatureExtractor_Barcode")])
  gridName = as.character(arrayInfo[2, which(tmpNames == "Grid_Name")])
  scanDate = as.character(arrayInfo[2, which(tmpNames == "Scan_Date")])
  scanDate <- gsub("(.*)-(.*)-(.*) (.*)+", "\\3-\\1-\\2", scanDate)
  programVersion = as.character(arrayInfo[2, which(tmpNames == "Protocol_Name")])
  gridGenomicBuild = as.character(arrayInfo[2, which(tmpNames == "Grid_GenomicBuild")])
  ref = 'Dual color hybridization'
  cat('\tDone.\n')
  return(c(barCode = barCode,
           gridName = gridName,
           scanDate = as.character(scanDate),
           programVersion = programVersion,
           gridGenomicBuild = gridGenomicBuild,
           reference = ref,
           analyseDate = format(Sys.Date(), "%Y-%m-%d")
           )
         )
}
.readAgilentMatrix <- function(filePath){
  cat('Reading values...')
  arrayInfo <- readLines(filePath, n = 15)
  startAt <- grep("^FEATURES", arrayInfo)
  cnSet <- read.csv(filePath, header = T, skip = startAt-1, sep = "\t", stringsAsFactors = FALSE)
  cat('\tDone.\n')
  cnSet <- .curateAgilentCnSet(cnSet)
  return(cnSet)
}
.getRFlags <- function(cnSet){
  flags <- which(cnSet$rIsSaturated == 1 |     			  # 1 = non valid rIsSaturated
                   cnSet$rIsFeatNonUnifOL == 1 | 			# 1 = non valid rIsFeatureNonUnifOL
                   cnSet$rIsWellAboveBG == 0)					# 0 = non valid rIsWellAboveBG
  cat(length(flags), 'flagged probes on chromosome', unique(cnSet$ChrNum), '\n')
  return(flags)
}
.getGFlags <- function(cnSet){
  flags <- which(cnSet$gIsSaturated == 1 |            # 1 = non valid gIsSaturated
                   cnSet$gIsFeatNonUnifOL == 1 | 			# 1 = non valid gIsFeatureNonUnifOL
                   cnSet$gIsWellAboveBG == 0) 				# 0 = non valid gIsWellAboveBG
  cat(length(flags), 'flagged probes on chromosome', unique(cnSet$ChrNum), '\n')
  return(flags)
}
.medFlag <- function(values, flagged, minpos, maxpos){
  sapply(flagged, function(f){
    ii <- max(minpos, f-8)
    jj <- min(maxpos, f+8)
    median(values[ii:jj], na.rm=TRUE)
  })
}
.replaceFlags <- function(cnSet){
  S <- split(cnSet, cnSet$ChrNum)

  cat("\nRed channel:\n")
  rflags <- sapply(S, function(subset) .getRFlags(subset))
  
  cat("\nGreen channel:\n")
  gflags <- sapply(S, function(subset) .getGFlags(subset))
  
  newR <- lapply(names(rflags), function(chr){
    chr <- as.numeric(chr)
    flagged <- rflags[[chr]]
    tmp <- S[[chr]]
    tmp$rMedianSignal[flagged] <- .medFlag(tmp$rMedianSignal, flagged, 1, nrow(tmp))
    as.numeric(tmp$rMedianSignal)
  })
  
  newG <- lapply(names(gflags), function(chr){
    chr <- as.numeric(chr)
    flagged <- gflags[[chr]]
    tmp <- S[[chr]]
    tmp$gMedianSignal[flagged] <- .medFlag(tmp$gMedianSignal, flagged, 1, nrow(tmp))
    as.numeric(tmp$gMedianSignal)
  })
  
  cnSet$rMedianSignal <- do.call(c, newR)
  cnSet$gMedianSignal <- do.call(c, newG)
  
  return(cnSet)
}
.suppressFlags <- function(object){
  if(grepl("Agilent", getInfo(object, 'platform'))){
    cat('Suppressing flagged probes...\n')
    cnSet <- getCNset(object)
    cnSet <- .replaceFlags(cnSet)
    flagNames <- c('gIsSaturated', 'rIsSaturated', 'gIsFeatNonUnifOL', 'rIsFeatNonUnifOL', 'gIsWellAboveBG', 'rIsWellAboveBG')
    cnSet <- cnSet[,-which(colnames(cnSet) %in% flagNames)]
    object@cnSet <- cnSet
  }
  cat("\n")
  return(object)
}
.suppressDuplic <- function(object){
  cnSet <- getCNset(object)
  cnSet <- cnSet[order(cnSet$ProbeName),]
  if (!any(colnames(cnSet) == 'ProbeName')) stop('None of the columns can be identifed as ProbeNames')
  dup <- duplicated(cnSet$ProbeName)
  if(any(dup)){
    cat('Suppressing duplicated probes...')
    duplicProbes <- as.character(unique(cnSet$ProbeName[dup]))
    duplicSet <- subset(cnSet, cnSet$ProbeName %in% duplicProbes)
    medianSet <- ddply(.data = duplicSet, .variables=.(ProbeName, SystematicName, ChrNum, ChrStart, ChrEnd), summarize,
                       rMedianSignal = median(rMedianSignal, na.rm=TRUE), gMedianSignal = median(gMedianSignal, na.rm=TRUE))
    cnSet <- rbind.data.frame(cnSet[-which(cnSet$ProbeName %in% duplicProbes),], medianSet)
  }
  cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
  rownames(cnSet) <- seq(1, nrow(cnSet))
  object@cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
  cat('Done.\n')
  return(object)
}
.preset <- function(object, Ksmooth=1000, Kmax=20, Nmin=Kmax*8, Mwidth=2, UndoSD=0.75, Alpha=1e-10){
  cat('Adding presettings...')
  if(grepl("Affymetrix", getInfo(object, 'platform'))){
    Ksmooth=7500; UndoSD=1
  }
  param <- list(Ksmooth = Ksmooth, Kmax = Kmax, Nmin = Nmin, Mwidth = Mwidth, UndoSD = UndoSD, Alpha = Alpha)
  object@param <- param
  cat('Done.\n')
  return (object)
}
.curateAgilentCnSet <- function(cnSet){
  keepCol <- which(as.character(colnames(cnSet)) %in%
                     c(	"ProbeName", "SystematicName",
                        "gMedianSignal", "rMedianSignal",
                        "gIsSaturated", "rIsSaturated",
                        "gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
                        "gIsWellAboveBG", "rIsWellAboveBG")
  )															
  cat('Filtering control probes...')
  isChr = grep('^chr[^Mrandom]*$', cnSet$SystematicName)
  cnSet <- cnSet[isChr, keepCol]
  cat('Done.\n')
  
  cat('Checking chr nums...')
  systNames <- cnSet$SystematicName
  chr <- gsub(":(.*)", "", systNames)
  chrNum <- gsub("(chr)(.*):(\\d+)-(\\d+)", "\\2", systNames)
  chrNum[chrNum=="X"] <- 23
  chrNum[chrNum=="Y"] <- 24
  chrNum <- as.numeric(chrNum)
  chrStart <- as.numeric(gsub("(chr)(.*):(\\d+)-(\\d+)", "\\3", systNames))
  chrEnd <- as.numeric(gsub("(chr)(.*):(\\d+)-(\\d+)", "\\4", systNames))
  cnSet <- cbind.data.frame(ProbeName = cnSet$ProbeName,
                            SystematicName = cnSet$SystematicName,
                            ChrNum=chrNum,
                            ChrStart=chrStart,
                            ChrEnd=chrEnd,
                            cnSet[,-c(1:2)]
                            )
  cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
  cat('Done.\n')
  
  return(cnSet)
}

################################
# Build a AffyCytoScan object
################################

buildAffyCytoScan <- function(filePath, sampleName=NA, labName=NA, useSNP=FALSE){
  fileName <- gsub("(.*)/", "", filePath)
  object <- new("cghObj", info = c(fileName=fileName,
                             sampleName=sampleName,
                             labName=labName,
                             platform='Affymetrix_CytoScanHD'
                             )
                    )
  affyData <- .readCytoScan(filePath, useSNP)
  object@info <- c(object@info, affyData$infos)
  object@cnSet <- affyData$values
  object <- .preset(object)
  return (object)
}
.getTagValue <- function(arrayInfos, tag){
  x <- arrayInfos[grep(tag, arrayInfos)]
  return(unlist(strsplit(x, '='))[2])
}
.readCytoScan <- function(filePath, useSNP){
  cat('Reading information...')
  fileName <- gsub("(.*)/", "", filePath)

  arrayInfos <- readLines(filePath, n = 1000)
  oldVersion <- any(grepl("#%affymetrix-array-type", arrayInfos))
  
  if(!oldVersion){
    startAt=1
    arrayType <- barCode <- gridName <- scanDate <- programVersion <- ucsc <- ensembl <- gridGenomicBuild <- ref <- NA    
  } else{
      arrayType <- .getTagValue(arrayInfos, "#%affymetrix-array-type")
      barCode <- .getTagValue(arrayInfos, "#%affymetrix-array-barcode")
      gridName <- .getTagValue(arrayInfos, "#%affymetrix-algorithm-param-state-annotation-file")
      scanDate <- .getTagValue(arrayInfos, "#%affymetrix-scan-date")
      programVersion <- .getTagValue(arrayInfos, "#%affymetrix-algorithm-version")
      ucsc <- .getTagValue(arrayInfos, "genome-version-ucsc")
      ensembl <- .getTagValue(arrayInfos, "genome-version-ensembl")
      gridGenomicBuild <- paste(ucsc, ensembl, sep = '/')
      ref <- .getTagValue(arrayInfos, "#%affymetrix-algorithm-param-state-reference-file")
      startAt <- grep("ProbeSetName", arrayInfos)
    }
  cat('\tDone.\n')
  infos <- c(barCode=barCode, gridName=gridName,
             scanDate=format(as.Date(scanDate), "%Y-%m-%d"), programVersion=programVersion,
             gridGenomicBuild=gridGenomicBuild, reference=ref,
             analyseDate=format(Sys.Date(), "%Y-%m-%d")
             )
  values <- .readCytoScanMatrix(filePath, oldVersion, startAt, useSNP)
  return(list(infos=infos, values=values))
}
.readCytoScanMatrix <- function(filePath, oldVersion, startAt, useSNP){
  cat('Reading values...')
  cnSet <- read.csv(filePath, header=TRUE, skip=startAt-1, sep="\t", stringsAsFactors=FALSE)
  colnames(cnSet) <- gsub("\\..*", "", colnames(cnSet))
  colnames(cnSet)[1:3] <- c("ProbeName", "ChrNum", "ChrStart")
  if(useSNP){
    cnSet <- cnSet[grep("S-\\d", cnSet$ProbeName),]
  } else{
    cnSet <- cnSet[grep("C-\\d", cnSet$ProbeName),]
  }
  cnSet$ChrNum <- .renameChr(cnSet$ChrNum)
#   if(!oldVersion){
#     if(any(cnSet$ChrNum == 24))
#       cnSet$ChrNum[cnSet$ChrNum == 24] <- 23
#     if(any(cnSet$ChrNum == 25))
#       cnSet$ChrNum[cnSet$ChrNum == 25] <- 24    
#   } else{
#     if(any(cnSet$ChrNum == "X"))
#       cnSet$ChrNum[cnSet$ChrNum == "X"] <- 23
#     if(any(cnSet$ChrNum == "Y"))
#       cnSet$ChrNum[cnSet$ChrNum == "Y"] <- 24
#   }
  for(i in 2:ncol(cnSet)) cnSet[,i] <- as.numeric(cnSet[,i])
  cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
  idx <- which(is.na(cnSet$ChrNum) | is.na(cnSet$ChrStart) | cnSet$WeightedLog2Ratio==0 | is.na(cnSet$SmoothSignal))
  cnSet <- cnSet[-idx,]
  cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
  cat('\tDone.\n')
  return(cnSet)
}
.renameChr <- function(ChrNum){
  if(any(ChrNum == "X"))
    ChrNum[ChrNum == "X"] <- 23
  if(any(ChrNum == "Y"))
    ChrNum[ChrNum == "Y"] <- 24
  # On new Affy ChAS version, chr23 and chr 24 are coded 24 and 25, resp.
  if(any(ChrNum == 25)){
    ChrNum[ChrNum == 24] <- 23
    ChrNum[ChrNum == 25] <- 24
  }
  return( as.numeric(ChrNum) )
}
############################
# Read SNP6 from local dir
############################
buildAffySNP6 <- function(filePath, sampleName=NA, labName=NA, useSNP=FALSE){
  fileName <- gsub("(.*)/", "", filePath)
  object <- new("cghObj", info = c(fileName=fileName,
                             sampleName=sampleName,
                             labName=labName,
                             synapseId=NA,
                             platform='Affymetrix_snp6'
                             )
                )
  affyData <- .readSNP6(filePath, useSNP)
  object@info <- c(object@info, affyData$infos)
  object@cnSet <- affyData$values
  object <- .preset(object)
  return (object)
}
.readSNP6 <- function(filePath, useSNP){
  cat('Reading information...')
#   if(any(grepl("bz2", filePath))){
#     destname <- gsub("[.]bz2$", "", filePath, ignore.case=TRUE)
#     arrayInfos <- readLines(gunzip(filePath, destname=destname, overwrite=TRUE, remove=FALSE), n = 750)    
#   } else{
#   }
  arrayInfos <- readLines(filePath, n = 750)
  
  arrayType <- .getTagValue(arrayInfos, "#ArraySet")
  barCode = NA
  gridName <- .getTagValue(arrayInfos, "#state-annotation-file")
  Date <- .getTagValue(arrayInfos, "#state-time-start")
  Date <- unlist(strsplit(Date, ' '))
  scanDate = paste(Date[5], Date[2], Date[3])#, sep = '-')
  programVersion <- .getTagValue(arrayInfos, "#option-program-version")
  ucsc <- .getTagValue(arrayInfos, "#genome-version-ucsc")
  ensembl <- .getTagValue(arrayInfos, "#genome-version-ncbi")
  gridGenomicBuild <- paste(ucsc, ensembl, sep = '/')
  ref <- .getTagValue(arrayInfos, "#state-reference-file")
  cat('\tDone.\n')
  infos <- c(barCode=barCode, gridName=gridName,
             scanDate=scanDate, programVersion=programVersion,
             gridGenomicBuild=gridGenomicBuild, reference=ref,
             analyseDate=format(Sys.Date(), "%Y-%m-%d")
             )
  startAt <- grep("ProbeSet", arrayInfos)
  values <- .readSNP6Matrix(filePath, startAt, useSNP)
  return(list(infos=infos, values=values))
}
.readSNP6Matrix <- function(filePath, startAt, useSNP){
  cat('Reading values...')
  cnSet <- read.csv(filePath, header=TRUE, skip=startAt-1, sep="\t", stringsAsFactors=FALSE)
  colnames(cnSet)[1:3] <- c("ProbeName", "ChrNum", "ChrStart")
  if(useSNP){
    cnSet <- cnSet[grep("^SNP_A-\\d+", cnSet$ProbeName),]
  } else{
    cnSet <- cnSet[grep("CN_\\d+", cnSet$ProbeName),]
  }
  cnSet$ChrNum <- .renameChr(cnSet$ChrNum)
#   cnSet$ChrNum[cnSet$ChrNum == "X"] <- 23
#   cnSet$ChrNum[cnSet$ChrNum == "Y"] <- 24
#  cnSet$ChrNum <- as.numeric(cnSet$ChrNum)
  cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
  cat('\tDone.\n')
  return(cnSet)
}
