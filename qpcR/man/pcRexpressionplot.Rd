\name{pcR.expressionplot}
\alias{pcR.expressionplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pcR.expressionplot
%%  ~~function to do ... ~~
}
\description{This function creates a bar plot PDF displaying the relative expression values of your genes of interest in relation to a house-keeping gene. The input file is generated from the pcr.expression() function. Careful: This function writes a .pdf in your current working directory.
}
\usage{
pcR.expressionplot(expression.file, title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression.file}{The .csv file that is generated using the pcr.expression() function. Please do not use any other file.
}
  \item{title}{Optional argument that provides your plot with a custom title. If there is no given title, the default title will be the filename of your expression.file.}
}
\details{
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
\author{Alexander Cunha}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##pcr.expressionplot(expression.file = 'file_path_expression.csv', title = 'Custom Plot Title')
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
