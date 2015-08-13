\name{plot.Path2PPI}
\alias{plot.Path2PPI}

\title{
Plots the predicted PPI
}
\description{
Plots the predicted PPI in three different ways. Depending on the \code{type} 
argument it manages the specific layout settings and finally uses the plot 
function of the \code{\link{igraph}} package.
}
\usage{
\method{plot}{Path2PPI}(x, type = "ppi", multiple.edges = FALSE, 
              scores = FALSE, species.colors = c(), 
              vertices.opacity=0.8, use.identifiers=FALSE,
              protein.labels = NA, show.legend = TRUE, 
              vertices.coordinates = NA, return.coordinates = FALSE, 
              tkplot=FALSE,...)
}
\arguments{
\item{x}{
An object from the class \code{Path2PPI} where the PPI network already has 
been predicted.
}
\item{type}{
Character string. Which graph type to plot. "ppi": plots only the predicted 
PPI. "hybrid": plots the hybrid network which consists of all relevant 
interactions from the reference species, the predicted interactions in the 
target species and all relevant homologous relationships.
}
\item{multiple.edges}{
Logical. Is only considered if \code{type="ppi"}. If \code{TRUE} then each 
reference interaction is depicted in the species-specific color (raw mode), 
in contrast, if set to \code{FALSE} only the finalized / combined interactions 
are depicted.
}
\item{scores}{
Logical. If \code{TRUE} the edge scores will be shown.
}
\item{species.colors}{
Named vector, to specify the species colors. If no value is given then default 
colors are used.
}
  \item{vertices.opacity}{
Numeric value between 0 and 1 defining the opacity of the vertices.
}
\item{use.identifiers}{
Logical. If \code{TRUE} then only the proteins identifiers are used as the 
vertex labels.
}
\item{protein.labels}{
Named vector to define the labels of the vertices. If no value is given then 
the protein identifiers are used. The vector does not have to be complete, 
i.e. not each protein has to be defined.
}
\item{show.legend}{
Logical. If \code{TRUE} then a legend is depicted.
}
\item{vertices.coordinates}{
Data frame containing the coordinates of the vertices. If no value is given 
then coordinates are computed using the \code{\link{layout.auto}} function.
}
\item{return.coordinates}{
Logical. If \code{TRUE} the coordinates of the vertices are returned.
}
\item{tkplot}{
Logical. If \code{TRUE} the graph is drawn in the interactive graph drawing
facility \code{\link{tkplot}}.
}
\item{\dots}{
Additional plotting parameters.
}
}
\details{
The argument \code{return.coordinates} only works correctly if 
\code{tkplot=FALSE}. If you want to get the coordinates of the tkplot device 
use \code{\link{tkplot.getcoords}}.
}

\note{
If you want to export the plotted graph to postscript you have to consider
that the default font family is set to \emph{sans} for vertex and edge labels.
Please change the default font family of postscript to \emph{sans} before you 
call the plot method: \code{ps.options(family="sans")}.
Additionally, you have to consider that the default value for 
\code{vertices.opacity} is set to \emph{0.8} in order to enhance the 
visibility of the graph, since some edges may be hidden by the vertices. 
Postscript does not support semi-transparencies. Hence, please change the 
\code{vertices.opacity} argument to \emph{1} if you want to export the graph 
using postscript.
}

\value{
If \code{return.coordinates=TRUE} the coordinates of the vertices are 
returned.
}

\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}

\seealso{
\code{\link{predictPPI}}, \code{\link{igraph}} for other plotting parameters
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

#Plot the predicted PPI with the default settings and return 
#the coordinates of the vertices
set.seed(12)
coordinates <- plot(ppi, return.coordinates=TRUE)

#Plot the predicted PPI and show each underlying reference interaction. 
#Use different species specific colors. To compare both graphs, 
#use the coordinates computed before
plot(ppi,multiple.edges=TRUE,vertices.coordinates=coordinates)

#Plot the corresponding hybrid network with predefined species colors.
#Also define some labels for the proteins of the target species.
#Keep in mind: You can not use the data in "coordinates" since 
#the hybrid network consists of more vertices than the default PPI
set.seed(40)
target.labels<-c("B2AE79"="PaTOR","B2AXK6"="PaATG1",
                "B2AUW3"="PaATG17","B2AM44"="PaATG11",
                "B2AQV0"="PaATG13","B2B5M3"="PaVAC8")
species.colors <- c("5145"="red","9606"="blue","559292"="green")
plot(ppi,type="hybrid",species.colors=species.colors,
protein.labels=target.labels)
}

