\name{getPPI}
\alias{getPPI}

\title{
Get predicted PPI
}
\description{
Get the predicted PPI of an \code{Path2PPI} object consisting of each 
predicted interaction and protein in the target species.
}
\usage{
getPPI(path2ppi, raw=FALSE, igraph=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path2ppi}{
An object of the class \code{Path2PPI}.
}
  \item{raw}{
Logical; if \code{TRUE} then the detailed view of the predicted PPI will be 
returned. That means that each predicted interaction deduced from each 
reference species is given. In contrast, \code{FALSE} leads to the actually 
predicted and combined PPI where no redundancies occur.
}
  \item{igraph}{
Logical; if \code{TRUE} then the returned PPI is given as igraph-object. 
Otherwise a data frame with each predicted interaction will be returned.
}
}

\value{
See \code{igraph} argument.
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}


\seealso{
\code{\link{getHybridNetwork}}
}
\examples{
data(ai) #Load test data set

ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")

ppi <- addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                    human.ai.irefindex, pa2human.ai.homologs)
ppi <- addReference(ppi, "Saccharomyces cerevisiae (S288c)", "559292", 
                    yeast.ai.proteins, yeast.ai.irefindex, 
                    pa2yeast.ai.homologs)

ppi <- predictPPI(ppi)

#Get the predicted PPI as data frame.
network <- getPPI(ppi)

#Get the detailed predicted PPI as data frame.
network.raw <- getPPI(ppi,raw=TRUE)
}

