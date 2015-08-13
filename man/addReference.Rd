\name{addReference}
\alias{addReference}

\title{
Add reference species
}
\description{
Adds reference species to an object from the class \code{Path2PPI}.
}
\usage{
addReference(path2ppi, taxName, taxId, proteins, irefindex, homologs)
}

\arguments{
\item{path2ppi}{
An object of the class \code{Path2PPI}.
}
\item{taxName}{
A character string giving the taxonomy name.
}
\item{taxId}{
A character string giving the taxonomy identifier.
}
\item{proteins}{
Either a character vector with the identifiers of the proteins which are 
involved in the corresponding pathway or a character vector with the protein 
names or aliases, respectively, named by the protein identifiers.
}
\item{irefindex}{
Either a data frame, representing the iRefIndex table of the current reference 
species, e.g. loaded previously via \code{read.table}, or the corresponding 
file name of the iRefIndex file.
}
\item{homologs}{
Either a data frame representing the results of the BLAST search (e.g. loaded 
previously via \code{read.table}) or the corresponding file name of the BLAST 
result file.
}
}
\details{
This method searches for all relevant interactions in the data frame or file 
defined in iRefIndex. There are different and often ambiguous protein 
identifiers defined in an iRefIndex file, and the putative "major" identifiers 
are not necessarily those defined in the corresponding "major" columns "uidA" 
and "uidB". Furthermore, iRefIndex also contains protein complexes. Hence, 
\code{Path2PPI} applies an advanced search algorithm to automatically find  
relevant interactions associated with the pathway or the proteins of interest, 
respectively. The user does not have to predefine the identifiers types 
(Uniprot, Swissprot, Ensemble etc.), since these types are often unambiguously 
assigned. The algorithm searches for each identifier in 10 columns where any 
type of identifier or accession number is defined ("uidA", "altA", 
"OriginalReferenceA", "FinalReferenceA", "aliasA", "uidB", "altB", 
"OriginalReferenceB", "FinalReferenceB" and "aliasB"). Additionally, it 
searches for each complex which contains one or more of the predefined 
proteins. Subsequently, each homologous relationship which is not relevant for 
the previously found interactions is declined. The results of these searches 
are centralized in the \code{Path2PPI} object and can be visualized using the 
appropriate methods (e.g. \code{\link{showReferences}})
}
\value{
An object from the class \code{Path2PPI} with attached reference species.
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}


\seealso{
\code{\link{showReferences}}, \code{\link{removeReference}}
}
\examples{
data(ai) #Load test data set

ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")
ppi

ppi <- addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                    human.ai.irefindex, pa2human.ai.homologs)
ppi <- addReference(ppi, "Saccharomyces cerevisiae (S288c)", "559292", 
                    yeast.ai.proteins, yeast.ai.irefindex, 
                    pa2yeast.ai.homologs) 
ppi
}

