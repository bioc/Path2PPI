\name{showInteraction}
\alias{showInteraction}
\alias{showInteraction-methods}

\title{
Information about an interaction
}
\description{
Use \code{showInteraction} to get detailed information about one 
interaction of the predicted PPI.
}
\usage{
showInteraction(path2ppi, interaction, mode="default", verbose=TRUE)
}
\arguments{
\item{path2ppi}{
An object from the class \code{Path2PPI}.
}
\item{interaction}{
Character vector consisting of the identifiers of the two interactors.
}
\item{mode}{
Character string. Which information of this interaction is requested. 
"default": only the predicted interaction and some major information are 
provided. "detailed": all interactions deduced from each reference species 
with this interaction is provided. "references": each reference interaction 
of the current interaction with some major information. "references.detailed": 
each reference interaction of the current interaction with all available 
information (extracted from the corresponding iRefIndex data set).
}
\item{verbose}{
Logical. \code{FALSE} hides messages in the output.
}
}

\value{
Data frame with the requested information defined in \code{mode}.
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}

\seealso{
\code{\link{plot.Path2PPI}},\code{\link{showReferences}}
}
\examples{
data(ai) #Load test data set

ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")
ppi <- addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                    human.ai.irefindex, pa2human.ai.homologs)
ppi <- addReference(ppi, "Saccharomyces cerevisiae (S288c)", "559292", 
                    yeast.ai.proteins, yeast.ai.irefindex, 
                    pa2yeast.ai.homologs)
                    
ppi <- predictPPI(ppi,h.range=c(1e-60,1e-20))

interaction <- showInteraction(ppi,interaction=c("B2AT71","B2AE79"),
mode="detailed")
interaction
}
