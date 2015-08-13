\name{getHybridNetwork}
\alias{getHybridNetwork}

\title{
Get hybrid network of the predicted PPI
}
\description{
Get the hybrid network of the previously predicted PPI. The hybrid network 
consists of all relevant interactions from the reference species, the 
predicted interactions in the target species and all relevant homologous 
relationships.
}
\usage{
getHybridNetwork(path2ppi, igraph = FALSE)
}

\arguments{
    \item{path2ppi}{
         An object of the class \code{Path2PPI}.
         }
    \item{igraph}{
         Logical; if \code{TRUE} then the hybrid network is given as 
         igraph-object. Otherwise a data frame, consisting of each 
         interaction and homologous relationship, will be returned.
    }
}

\value{
See \code{igraph} argument.
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}


\seealso{
\code{\link{getPPI}}
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

#Return the hybrid network as data frame
hybrid <- getHybridNetwork(ppi)

#Return the hybrid network as igraph object
hybrid <- getHybridNetwork(ppi,igraph=TRUE)
}

