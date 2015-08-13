\name{showReferences}
\alias{showReferences}

\title{
Information about reference species
}
\description{
Get information about the currently stored reference species. If indicated by 
\code{returnValue} a data frame - containing information about each protein or 
interaction - is provided as well.
}
\usage{
showReferences(path2ppi, species = NA, returnValue = NA)
}

\arguments{
\item{path2ppi}{
An object from the class \code{Path2PPI}.
}
\item{species}{
Either a number between 1 and the number of stored reference species or a 
character string with the taxonomy id. If no value for \code{species} is given 
then information about each stored reference species is provided.

}
\item{returnValue}{
Character value indicating whether to return a value. "proteins": a data frame 
containing the proteins associated with the pathway of interest in the  
corresponding reference species. "interactions": a data frame containing all 
processed, relevant and non-redundant interactions. "irefindex": a data frame 
containing all relevant interactions in the raw irefindex format. Is only 
reasonable if \code{species} is defined. If no value for \code{returnValue} is 
given then only general information is provided.
}
}

\value{
See description for \code{returnValue}
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}

\seealso{
\code{\link{addReference}}, \code{\link{removeReference}}, 
\code{\link{showInteraction}}
}
\examples{
data(ai) #Load test data set

ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")

ppi <- addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                    human.ai.irefindex, pa2human.ai.homologs)
ppi <- addReference(ppi, "Saccharomyces cerevisiae (S288c)", "559292", 
                    yeast.ai.proteins, yeast.ai.irefindex, 
                    pa2yeast.ai.homologs) 

#Get general information about each stored reference species
showReferences(ppi)

#Get general information about reference species with the taxonomy id "9606"
showReferences(ppi, species="9606")

#Get all proteins associated with the pathway of interest 
#and previously given by the user
proteins <- showReferences(ppi, species="9606", returnValue="proteins")

#Get all processed and non-redundant interactions previously 
#determined to be relevant for the pathway of interest
interactions <- showReferences(ppi, species="9606", returnValue="interactions")

#Get all relevant interactions in the detailed irefindex format
irefindex <- showReferences(ppi, species="9606", returnValue="irefindex")


}

