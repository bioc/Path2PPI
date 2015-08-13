\name{removeReference}
\alias{removeReference}

\title{
Remove reference species
}
\description{
Remove reference species previously attached to an object from the class 
\code{Path2PPI}.
}
\usage{
removeReference(path2ppi, species)
}

\arguments{
  \item{path2ppi}{
An object from the class \code{Path2PPI}.
}
  \item{species}{
Either a number between 1 and the number of stored reference species or a 
character string with the taxonomy id of the reference species to remove.
}
}

\value{
An object of the class \code{Path2PPI} with removed reference species 
\code{species}.
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}

\seealso{
\code{\link{showReferences}}, \code{\link{addReference}}
}
\examples{
data(ai) #Load test data set

ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")

ppi <- addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                    human.ai.irefindex, pa2human.ai.homologs)
ppi <- addReference(ppi, "Saccharomyces cerevisiae (S288c)", "559292", 
yeast.ai.proteins, yeast.ai.irefindex, pa2yeast.ai.homologs) 

#Remove second reference species
ppi <- removeReference(ppi,2)

#Remove reference species with taxonomy id "9606"
ppi <- removeReference(ppi,"9606") 

}
