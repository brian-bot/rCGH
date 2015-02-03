\name{buildAgilent}
\alias{buildAgilent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buildAgilent(filePath, Ref = "cy3", sampleName = NA, labName = NA, supFlags = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filePath}{
%%     ~~Describe \code{filePath} here~~
}
  \item{Ref}{
%%     ~~Describe \code{Ref} here~~
}
  \item{sampleName}{
%%     ~~Describe \code{sampleName} here~~
}
  \item{labName}{
%%     ~~Describe \code{labName} here~~
}
  \item{supFlags}{
%%     ~~Describe \code{supFlags} here~~
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
function (filePath, Ref = "cy3", sampleName = NA, labName = NA, 
    supFlags = TRUE) 
{
    fileName <- gsub("(.*)/", "", filePath)
    object <- new("cghObj", info = c(fileName = fileName, sampleName = sampleName, 
        labName = labName, platform = "Agilent"))
    object@info <- c(object@info, .readAgilentInfo(filePath))
    object@cnSet <- .readAgilentMatrix(filePath)
    if (supFlags) 
        object <- .suppressFlags(object)
    object <- .suppressDuplic(object)
    object <- .preset(object)
    return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
