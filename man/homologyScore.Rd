\name{homologyScore}
\alias{homologyScore}

\title{
Computes homology scores
}
\description{
Computes the homology scores based on the BLAST E-value. This function is used 
by the \code{predictPPI} method to compute homology scores to decide whether 
an interaction in a reference species is adopted to the target species 
(see package vignette for a detailed description). It can be used to test 
which E-values lead to which scores given a predefined E-value range.
}
\usage{
homologyScore(e.value, h.range)
}

\arguments{
    \item{e.value}{
    One BLAST E-value or a numeric vector with different BLAST E-values
    }
    \item{h.range}{
    Numeric vector consisting of two values. The first value indicates the 
    lower bound (smallest E-value). Each E-value which is equal or less 
    than this bound is scored with 1. The second value indicates the 
    upper bound (biggest E-value). Each E-value which is equal or greater than 
    this bound is scored with 0.
    }
}
\details{
Uses a linear function to map the E-value \eqn{v} to the range \eqn{[l,u]} 
where \eqn{l} is the lower and \eqn{u} the upper bound:


\deqn{s(v)=|m\log_{10}(v)+b|}{s(v)=|m*log_10(v)+b|}
\deqn{m=\frac{1}{\log_{10}(l)-\log_{10}(u)}}{m=1/(log_10(l)-log_10(u))}
\deqn{b=-(m\log_{10}(u))}{b=-(m*log_10(u))}


}
\value{
Numeric vector containing the scores.
}
\author{
Oliver Philipp \email{MolBI-software@bioinformatik.uni-frankfurt.de}
}

\seealso{
\code{\link{predictPPI}}
}
\examples{
l <- 1e-100                                   #lower bound
u <- 1e-20                                    #upper bound
h.range <- c(l,u)                             #define range
e.values <- c(1e-20,1e-40,1e-60,1e-80,1e-100) #some BLAST E-values

homologyScore(e.values,h.range)
}


