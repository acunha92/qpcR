\name{pcR.expression}
\alias{pcR.expression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pcR.expression
%%  ~~function to do ... ~~
}
\description{Creates a .csv that contains the deltaCT values as well as Relative Expression of your genes of interest in relation to your house-keeping gene of choice. Input file is generated from the pcr.summary() function. Careful: This function writes a .csv to your current working directory.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pcR.expression(summary.file, goi, hkg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summary.file}{The .csv file that is generated using the pcr.summary() function. Please do not use any other file.
%%     ~~Describe \code{x} here~~
}
  \item{goi}{Either a character string or a vector of character strings indicating your genes of interest to be tested against a house-keeping gene.
}
  \item{hkg}{A character string indicating the house-keeping gene you want to test your gene(s) of interest against. Please only run the function with one house-keeping gene at a time. The output filename will denote which house-keeping gene you chose.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##pcr.expression(summary.file = 'file_path_summary.csv', goi = c('GENE1','GENE2','GENE3'), hkg = 'HKG1'))
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
