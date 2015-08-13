\name{predictPPI}
\alias{predictPPI}

\title{
Prediction of the PPI
}
\description{
Major method of the \code{Path2PPI} class to predict the final PPI in the 
target species using the information available from the stored reference 
species. Different values for the arguments of this method can lead to 
different PPI networks, differing in the degree of reliability and strictness.
}
\usage{
predictPPI(path2ppi, mode="both", h.thresh=1e-05, 
           h.range=c(1e-100, 1e-20), i.thresh=0.7, 
           consider.complexes=FALSE, max.complex.size=5,
           decline.self.interaction.ref=FALSE,
           decline.self.interaction.tar=TRUE,
           verbose=TRUE)
}

\arguments{
  \item{path2ppi}{
An object of the class \code{Path2PPI}.
}
  \item{mode}{
Which interaction from the reference species should be taken into account. 
"both": both interactors of an interaction has to be in the initial protein 
list previously inserted by the user (recommended if it is a large network or 
many proteins were initially defined, respectively). "one": only one of the 
interactors of each reference interaction has to be in the initial protein 
list (may lead to very large networks).
}
  \item{h.thresh}{
E-value cutoff at which each homologous relationship definitely will be 
declined (see also \code{h.range} argument).
}
  \item{h.range}{
Numeric vector consisting of two values. The first value indicates the lower 
border (smallest E-value). Each E-value which is equal or less than this 
border is scored with 1 (best). The second value indicates the upper border 
(biggest E-value). Each E-value which is equal or greater than this border is 
scored with 0 (worst). 
}
  \item{i.thresh}{
Numeric. Threshold for accepted interactions. If the computed prediction score 
for an interaction is less than \code{i.thresh} it will be declined.
}
  \item{consider.complexes}{
Logical. If \code{TRUE} then interactions are also considered which actually 
indicate an association of the current protein to one bigger protein complex. 
This may lead to very large networks if \code{mode="one"} since all other 
proteins of this complex are considered as well, i.e., each protein in such 
complexes are considered to interact with each other protein of this complex. 
If \code{mode="both"} then each protein of an complex has to be in the initial 
protein list to consider each interaction (see details).
}
  \item{max.complex.size}{
Numeric. Is only considered if \code{consider.complexes=TRUE}. The maximum 
size of complexes to be considered. 
}
  \item{decline.self.interaction.ref}{
Logical. If \code{TRUE} then all self interactions from reference species are 
declined.
}
  \item{decline.self.interaction.tar}{
Logical. If \code{TRUE} then all predicted self interactions in target species 
are declined.
}
  \item{verbose}{
Logical. \code{FALSE} hides messages in the output.
}
}
\details{
Difference of \code{h.thresh} and \code{h.range}: If only one protein in the 
target species was found to be homologous to a current reference species 
protein and this homology was rated with an E-value which is equal or smaller 
than \code{h.thresh} it is scored with 1 (even if the E-value is larger than 
the upper border of \code{h.range}). See package vignette for more details.

Use the complex arguments with care, since each complex may lead to a vast 
amount of interactions, i.e., each protein is considered to interact with each 
other of this complex; e.g. if there are 10 proteins involved in one complex, 
this would lead to 10 over 2 = 45 interactions.
}
\value{
An object of the class \code{Path2PPI} with predicted PPI.
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}

\seealso{
\code{\link{plot.Path2PPI}}, \code{\link{homologyScore}}
}
\examples{
data(ai) #Load test data set

ppi <- Path2PPI("Autophagy induction", "Podospora anserina", "5145")
ppi <- addReference(ppi, "Homo sapiens", "9606", human.ai.proteins, 
                    human.ai.irefindex, pa2human.ai.homologs)
ppi <- addReference(ppi, "Saccharomyces cerevisiae (S288c)", "559292", 
                    yeast.ai.proteins, yeast.ai.irefindex, 
                    pa2yeast.ai.homologs)

#Using the default settings leads to 8 predicted interactions in the 
#target species
ppi <- predictPPI(ppi)

#Consider complexes where each complex is allowed to be up to 10 proteins 
#large. For this smaller pathway only one more interaction was predicted when 
#considering larger complexes.
ppi <- predictPPI(ppi,consider.complexes=TRUE,max.complex.size=10)

#We can be less strict and decrease h.range what obviously increases the 
#number of predicted interactions to 13
ppi <- predictPPI(ppi,h.range=c(1e-60,1e-20))
}
