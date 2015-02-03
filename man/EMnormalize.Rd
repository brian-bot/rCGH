\name{EMnormalize}
\alias{EMnormalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EMnormalize(object, cut = c(0.1, 0.9), G = 3:6, B = 100, peakThresh = 0.5, ksmooth = 801, useN = 2000, mergeVal = 0, Title = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{cut}{
%%     ~~Describe \code{cut} here~~
}
  \item{G}{
%%     ~~Describe \code{G} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{peakThresh}{
%%     ~~Describe \code{peakThresh} here~~
}
  \item{ksmooth}{
%%     ~~Describe \code{ksmooth} here~~
}
  \item{useN}{
%%     ~~Describe \code{useN} here~~
}
  \item{mergeVal}{
%%     ~~Describe \code{mergeVal} here~~
}
  \item{Title}{
%%     ~~Describe \code{Title} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, cut = c(0.1, 0.9), G = 3:6, B = 100, peakThresh = 0.5, 
    ksmooth = 801, useN = 2000, mergeVal = 0, Title = NA) 
{
    "\n            Use a EM algorithm to model the LogRatio density as a mixture of gaussian distributions\n            "
    cnSet <- getCNset(object)
    testLR <- cnSet$Log2Ratio
    cut <- as.numeric(quantile(testLR, c(cut[1], cut[2])))
    runLR <- .smoothLR(testLR, getInfo(object, "platform"), cut, 
        K = ksmooth)
    CHR <- cnSet$ChrNum[runLR$index]
    runLR <- runLR$runLR
    cat("Analyzing mixture:")
    EMmodels <- lapply(1:B, function(b) {
        cat(".")
        .buildEMmodel(runLR, G, Len = useN)
    })
    nG <- do.call(c, lapply(EMmodels, function(m) m$nG))
    tnG <- table(nG)
    mednG <- as.numeric(names(tnG)[which.max(tnG)])
    models <- EMmodels[nG == mednG]
    m <- do.call(rbind, lapply(models, function(m) m$m))
    p <- do.call(rbind, lapply(models, function(m) m$p))
    s <- do.call(rbind, lapply(models, function(m) m$s))
    cat("Done.\n")
    m <- apply(m, 2, median, na.rm = TRUE)
    s <- apply(s, 2, median, na.rm = TRUE)
    p <- apply(p, 2, median, na.rm = TRUE)
    if (mergeVal > 0) {
        mergedPars <- .mergePeaks(mednG, length(runLR), m, s, 
            p, mergeVal)
        m <- mergedPars$m
        s <- mergedPars$s
        p <- mergedPars$p
    }
    cat("Computing densities...")
    computeD <- .computeDensities(length(runLR), m, p, s)
    dList <- computeD$dList
    peaks <- computeD$peaks
    cat("Done.\n")
    cat("Gaussian mixture:\n")
    cat("n.peaks = ", mednG, "\n")
    bestPeak <- .chooseBestPeak(peaks, m, peakThresh)
    correct <- m[bestPeak]
    cat("\nGroup parameters:\n")
    for (grp in 1:length(m)) {
        cat("Grp", grp, "\nprop:", p[grp], ":\tmean:", m[grp], 
            "\tSd:", sqrt(s[grp]), "\tCV:", signif(sqrt(s[grp])/abs(m[grp]), 
                4), "\n")
    }
    cat("\nEM correction factor = ", correct, "\n")
    if (is.na(Title)) 
        Title <- sprintf("\%s\nEM centralization: correction.factor = \%s", 
            getInfo(object, "sampleName"), round(correct, 5))
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
    cat("Use getDensity(cghObject) to visualize the probes density plot.\n")
    return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
