\name{Path2PPI-class}
\Rdversion{1.1}
\docType{class}
\alias{Path2PPI-class}
\alias{Path2PPI}
\alias{addReference,Path2PPI-method}
\alias{getHybridNetwork,Path2PPI-method}
\alias{getPPI,Path2PPI-method}
\alias{initialize,Path2PPI-method}
\alias{plot.Path2PPI,Path2PPI-method}
\alias{predictPPI,Path2PPI-method}
\alias{removeReference,Path2PPI-method}
\alias{show,Path2PPI-method}
\alias{showInteraction,Path2PPI-method}
\alias{showReferences,Path2PPI-method}


\title{Class \code{"Path2PPI"}}
\description{
An instance of the class \code{Path2PPI} is the major object in the Path2PPI 
package. It manages all reference species and the target species. 
The prediction algorithm is implemented in this class as well.
}
\usage{
Path2PPI(...)
}

\arguments{
  \item{...}{Argument list (see Note below).}
}


\section{Slots}{
  \describe{
    \item{\code{pathway}:}{Object of class \code{"character"} }
    \item{\code{targetSpecies}:}{Object of class \code{".TargetSpecies"} }
    \item{\code{referenceContainer}:}{Object of class 
    \code{".ReferenceContainer"} }
    \item{\code{h.thresh}:}{Object of class \code{"numeric"} }
    \item{\code{h.range}:}{Object of class \code{"numeric"} }
    \item{\code{i.thresh}:}{Object of class \code{"numeric"} }
    \item{\code{consider.complexes}:}{Object of class \code{"logical"} }
    \item{\code{max.complex.size}:}{Object of class \code{"numeric"} }
    \item{\code{raw.ppi}:}{Object of class \code{"data.frame"} }
    \item{\code{ppi}:}{Object of class \code{"data.frame"} }
  }
}

\section{Methods}{
  \describe{
    \item{\code{\link{addReference}}}{\code{signature(path2ppi = "Path2PPI")} }
    \item{\code{\link{getHybridNetwork}}}{
    \code{signature(path2ppi = "Path2PPI")} }
    \item{\code{\link{getPPI}}}{\code{signature(path2ppi = "Path2PPI")} }
    \item{\code{\link{initialize}}}{\code{signature(.Object = "Path2PPI")} }
    \item{\code{\link{plot.Path2PPI}}}{\code{signature(x = "Path2PPI")} }
    \item{\code{\link{predictPPI}}}{\code{signature(path2ppi = "Path2PPI")} }
    \item{\code{\link{removeReference}}}{
    \code{signature(path2ppi = "Path2PPI")} }
    \item{\code{\link{show}}}{\code{signature(object = "Path2PPI")} }
    \item{\code{\link{showInteraction}}}{
    \code{signature(path2ppi = "Path2PPI")} }
    \item{\code{\link{showReferences}}}{
    \code{signature(path2ppi = "Path2PPI")} }
   }
}


\value{
An instance of the class \code{Path2PPI}.
}


\note{
Arguments to \code{Path2PPI()} and the \code{\link{new}} method are obligatory 
and must be named if they differ from this order:
\tabular{ll}{
  \strong{pathway} \tab A character string with the name of the pathway which 
  has to be predicted.\cr
  \strong{targetName} \tab A character string giving the taxonomy name of the 
  target species.\cr
  \strong{targetId} \tab A character string giving the taxonomy identifier of 
  the target species.\cr
}
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}

\examples{
ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")
ppi
}
