\name{ai}
\alias{ai}
\alias{human.ai.irefindex}
\alias{human.ai.proteins}
\alias{pa2human.ai.homologs}
\alias{yeast.ai.irefindex}
\alias{yeast.ai.proteins}
\alias{pa2yeast.ai.homologs}
\docType{data}
\title{
Data set to predict autophagy induction in \emph{Podospora anserina}
}
\description{
This data set consists of all data files necessary to predict the putative 
interactions of the induction step of autophagy in \emph{Podospora anserina} 
by means of the corresponding PPIs in human and yeast.
}
\usage{data("ai")}
\format{
    human.ai.irefindex: Data frames with 1694 observations of 54 variables.
    yeast.ai.irefindex: Data frames with 3840 observations of 54 variables.
    pa2human.ai.homologs: Data frames with 261 observations of 12 variables.
    pa2yeast.ai.homologs: Data frames with 98 observations of 12 variables.
    human.ai.proteins: Named character vector with 5 elements.
    yeast.ai.proteins: Named character vector with 7 elements.
}
\details{
    Data frames \code{human.ai.irefindex} and \code{yeast.ai.irefindex} 
    consists of all relevant interactions of the corresponding iRefIndex 
    files. The two data frames \code{pa2human.ai.homologs} and 
    \code{pa2yeast.ai.homologs} are the necessary parts of the result files 
    from the BLAST searches of the \emph{P. anserina} proteom against the 
    proteoms of human and yeast. The named character vectors 
    \code{human.ai.proteins} and \code{yeast.ai.proteins} consists of the 
    proteins involved in the induction process of autophagy in human and 
    yeast.
}

\references{
Camacho, C. et al. (2009). BLAST+: architecture and applications. BMC 
Bioinformatics, 10(1), 421.

Razick, S. et al. (2008). iRefIndex: a consolidated protein interaction 
database with provenance. BMC Bioinformatics, 9(1), 405.
}
\examples{
data(ai)
}
\value{
Four data frames and two named character vectors (see above).
}

\keyword{datasets}
