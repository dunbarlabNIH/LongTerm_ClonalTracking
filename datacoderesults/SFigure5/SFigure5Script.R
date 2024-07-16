#Rohan Hosuru and Jack Yang (based off of Samson Koelle's Code)
#Copywrite Cynthia Dunbar Group
#Modified version of R code for Blood paper (Koelle et al. 2017)
#This script loads data and calls source code to recreate results from the paper.
#You should be able to copy and paste it into your R console to get results (after installing R packages and setting working directories)

#First, we set some parameters
ppi = 300
size = 8
set.seed(1)

#set the following variables
workingdir = "/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/results"
sourcedir = '/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/source'

#download archived version of DiversitySampler package
require(devtools)
install_version("DiversitySampler", version = "2.1", repos = "http://cran.us.r-project.org")

#load required packages
library(pheatmap) 
library(DiversitySampler)
library(RColorBrewer)
library(nnet)
library(foreach)
library(stringr)
library(biclust)
library(scales)

#Load functions. Function descriptions are available within the functions themselves
setwd(sourcedir)
#thresholding functions
source('basethresholdinsamples.R')
source('thresholdinsamples.R')

#functions for Figure 2
source('GFPtrack.R')
source('diversitytrack.R')
source('clonenumbertrack.R')

#functions for Figure 3
source("BCLH v1.3.R")

#functions for Figure 4
source("BCLH v1.3.R")
source("cor_analysisv1.3.R")

#functions for Figure 5
source('biastracking3.R')
source('rollingvarianceofbias.R')

#Load Data
setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/loaders")
source("zj31 loader.R")
source("zh19 loader.R")
source("zj38 loader.R")

#Make Figures
setwd(workingdir)
dirname = paste("Figure S5", Sys.time())
dir.create(dirname)
setwd(dirname)

#Top 10 Clones across all four lineages
BCLH(data = zh19alldata[,as.vector(zh19indexmatrixnoNA)], names = zh19allnames[as.vector(zh19indexmatrixnoNA)],n_clones = 10,folder = paste("ZH19"))
BCLH(data = zj31alldata[,as.vector(zj31indexmatrixnoNA)], names = zj31allnames[as.vector(zj31indexmatrixnoNA)],n_clones = 10,folder = paste("ZJ31"))
BCLH(data = zj38alldata[,as.vector(zj38indexmatrixnoNA)], names = zj38allnames[as.vector(zj38indexmatrixnoNA)],n_clones = 10,folder = paste("ZJ38"))

#Correlation Graphs
cor_analysis(data = zh19alldata, indexmatrix = zh19indexmatrix, timepoints = zh19timepoints,folder = paste("ZH19"))
cor_analysis(data = zj31alldata, indexmatrix = zj31indexmatrix, timepoints = zj31timepoints,folder = paste("ZJ31"))
cor_analysis(data = zj38alldata, indexmatrix = zj38indexmatrix, timepoints = zj38timepoints,folder = paste("ZJ38"))

setwd(workingdir)