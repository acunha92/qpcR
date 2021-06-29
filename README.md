# qpcR
Repository for the qpcR package. Work in progress. 

Code written by Alexander Cunha


To install: 

	install.packages('devtools')

	library(devtools)

	devtools::install_github('acunha92/qpcR', subdir = 'qpcR')
	
	library(qpcR)




Requires: tidyverse, reshape



	pcR.summary(file, replicates),

	pcR.meanplot(summary.file, title),

	pcR.expression(summary.file, goi, hkg),

	pcR.expressionplot(expression.file, title)

