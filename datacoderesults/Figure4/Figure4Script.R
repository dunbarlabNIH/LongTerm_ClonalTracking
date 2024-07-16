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
source("zh33 loader.R")
source("zg66 loader.R")

#Make Figures
setwd(workingdir)
dirname = paste("Figure 4 Images", Sys.time())
dir.create(dirname)
setwd(dirname)

#Top 10 Clones per lineage
BCLH(data = zh33alldata[,as.vector(zh33indexmatrixnoNA)], names = zh33allnames[as.vector(zh33indexmatrixnoNA)], n_clones = 10, folder = paste("ZH33"))
BCLH(data = zg66alldata[,as.vector(zg66indexmatrixnoNA)], names = zg66allnames[as.vector(zg66indexmatrixnoNA)], n_clones = 10,folder = paste("ZG66"))

#Correlation Graphs
cor_analysis(data = zh33alldata, indexmatrix = zh33indexmatrix, timepoints = zh33timepoints, folder = "ZH33" , celltypenames = c("T", "B","Mono",'Gr'))
cor_analysis(data = zg66alldata, indexmatrix = zg66indexmatrix, timepoints = zg66timepoints, folder = "ZG66" , celltypenames = c("T", "B","Mono",'Gr'))

setwd(workingdir)