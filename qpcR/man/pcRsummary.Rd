\name{pcR.summary}
\alias{pcR.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pcR.summary
%%  ~~function to do ... ~~
}
\description{The most important function for the qpcR package. Input .csv must be formatted properly. Correct usage of the function will generate a table that contains the mean Ct values and standard deviation for all genes in your dataset. Careful: This function writes a .csv to your current working directory.
}
\usage{
pcR.summary(file, replicates)
}
\arguments{
  \item{file}{File path for the .csv file that contains your raw Ct values. Data MUST be composed of the following columns (in any order): Gene, Sample, Ct, Replicate.
}
  \item{replicates}{A numeric argument that specifies the number of replicates you are using in your dataset.
}
}
\details{The input file MUST be saved as .csv! To maintain stability, use one-word character strings for your Gene and Sample.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Alexander Cunha
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##pcR.summary(file = 'file_path.csv', replicates = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
