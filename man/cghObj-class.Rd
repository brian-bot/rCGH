\name{cghObj-class}
\Rdversion{1.1}
\docType{class}
\alias{cghObj-class}
\alias{getByGene,cghObj-method}
\alias{getCNset,cghObj-method}
\alias{getDensity,cghObj-method}
\alias{getInfo,cghObj-method}
\alias{getParam,cghObj-method}
\alias{getProfile,cghObj-method}
\alias{getSegTable,cghObj-method}
\alias{setInfo,cghObj-method}
\alias{show,cghObj-method}
\alias{view,cghObj-method}

\title{Class \code{"cghObj"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cghObj", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{info}:}{Object of class \code{"character"} ~~ }
    \item{\code{cnSet}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{param}:}{Object of class \code{"vector"} ~~ }
    \item{\code{segTable}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{byGene}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{probesDensity}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{gProfile}:}{Object of class \code{"ANY"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getByGene}{\code{signature(object = "cghObj")}: ... }
    \item{getCNset}{\code{signature(object = "cghObj")}: ... }
    \item{getDensity}{\code{signature(object = "cghObj")}: ... }
    \item{getInfo}{\code{signature(object = "cghObj")}: ... }
    \item{getParam}{\code{signature(object = "cghObj")}: ... }
    \item{getProfile}{\code{signature(object = "cghObj")}: ... }
    \item{getSegTable}{\code{signature(object = "cghObj")}: ... }
    \item{setInfo}{\code{signature(object = "cghObj")}: ... }
    \item{show}{\code{signature(object = "cghObj")}: ... }
    \item{view}{\code{signature(object = "cghObj")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("cghObj")
}
\keyword{classes}
