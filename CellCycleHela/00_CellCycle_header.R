#' # Load CCprofiler package and set working directory
library(devtools)
options(warn=-1)

install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
library('CCprofiler')
library('data.table')
library('ggplot2')
library('fitdistrplus')

setwd("/Users/isabell/Desktop/projects/ProteoformProject/Results/")

