\name{adjustSignal}
\alias{adjustSignal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
adjustSignal(object, Scale = TRUE, Cy = TRUE, Ref = "cy3", GC = TRUE, gcDB = agilentDB, preSeg = TRUE, suppOutliers = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{Scale}{
%%     ~~Describe \code{Scale} here~~
}
  \item{Cy}{
%%     ~~Describe \code{Cy} here~~
}
  \item{Ref}{
%%     ~~Describe \code{Ref} here~~
}
  \item{GC}{
%%     ~~Describe \code{GC} here~~
}
  \item{gcDB}{
%%     ~~Describe \code{gcDB} here~~
}
  \item{preSeg}{
%%     ~~Describe \code{preSeg} here~~
}
  \item{suppOutliers}{
%%     ~~Describe \code{suppOutliers} here~~
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
function (object, Scale = TRUE, Cy = TRUE, Ref = "cy3", GC = TRUE, 
    gcDB = agilentDB, preSeg = TRUE, suppOutliers = FALSE) 
{
    cnSet <- getCNset(object)
    if (grepl("Agilent", getInfo(object, "platform"))) {
        if (Cy) {
            cnSet <- .CyAdjust(cnSet, Ref)
            object@param$CyAdjusted = TRUE
        }
        if (GC) {
            cnSet <- .GCadjust(cnSet, gcDB)
            object@param$GCAdjusted = TRUE
        }
    }
    else {
        object@param$CyAdjusted = FALSE
        object@param$GCAdjusted = FALSE
    }
    object@param$dLRs <- .dlrs(cnSet$Log2Ratio)
    object@param$MAD <- .MAD(cnSet$Log2Ratio)
    if (Scale) {
        cat("Scaling...")
        cnSet$Log2Ratio <- scale(cnSet$Log2Ratio, center = FALSE)
        if (grepl("Affymetrix", getInfo(object, "platform"))) 
            cnSet$Log2Ratio <- cnSet$Log2Ratio * 1.25
        cat("Done.\n")
    }
    if (suppOutliers) 
        cnSet$Log2Ratio <- .supprOutliers(cnSet$Log2Ratio, n = ifelse(grepl("Agilent", 
            getInfo(object, "platform")), 5, 25))
    if (preSeg) {
        cnSet <- .preSeg(cnSet, getInfo(object, "sampleName"), 
            getParam(object))
    }
    object@cnSet <- cnSet
    cat("dLRs:", object@param$dLRs, "\tMAD:", object@param$MAD, 
        "\n")
    return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
