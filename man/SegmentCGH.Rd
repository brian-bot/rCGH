\name{SegmentCGH}
\alias{SegmentCGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SegmentCGH(object, Smooth = TRUE, UndoSD = NULL, minMarks = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{Smooth}{
%%     ~~Describe \code{Smooth} here~~
}
  \item{UndoSD}{
%%     ~~Describe \code{UndoSD} here~~
}
  \item{minMarks}{
%%     ~~Describe \code{minMarks} here~~
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
function (object, Smooth = TRUE, UndoSD = NULL, minMarks = 8) 
{
    cnSet <- getCNset(object)
    cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
    params = getParam(object)
    if (is.null(UndoSD)) 
        params$UndoSD <- 0.25 + sqrt(max(getParam(object3)$SigmaSq))
    else params$UndoSD <- UndoSD
    L2R <- cnSet$Log2Ratio
    Chr <- cnSet$ChrNum
    Pos <- cnSet$ChrStart
    sampleName <- getInfo(object, "sampleName")
    if (is.na(sampleName)) 
        sampleName <- "sample_x"
    cat("Computing standard segmentation...\n")
    segTable <- .computeSegmentation(L2R, Chr, Pos, sampleName, 
        params, Smooth)
    segTable <- .smoothSeg(segTable, minMarks)
    segTable <- .computeMedSegm(segTable, L2R)
    segTable <- .mergeLevels(segTable)
    probeValues.left <- .probeSegValue(segTable, use.medians = TRUE)
    cat("Number of segments:", nrow(segTable), "\n")
    object@param <- params
    object@segTable = segTable
    object@cnSet = cbind.data.frame(cnSet, Segm = probeValues.left)
    return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
