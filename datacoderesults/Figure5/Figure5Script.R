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

#Make Figures
setwd(workingdir)
dirname = paste("Figure 5 Images", Sys.time())
dir.create(dirname)
setwd(dirname)
#make bias tracking plots for ZH33
biastracking(data = zh33alldata, uppersamples = zh33indexmatrix[,2],lowersamples =  zh33indexmatrix[,4], timepointnames = zh33timepointnames,folder = "ZH33 B v Gr",upperbiaslabel = "B", lowerbiaslabel = 'Gr')
biastracking(data = zh33alldata, uppersamples = zh33indexmatrix[,1],lowersamples =  zh33indexmatrix[,2], timepointnames = zh33timepointnames,folder = "ZH33 T v B",upperbiaslabel = "T", lowerbiaslabel = 'B')
biastracking(data = zh33alldata, uppersamples = zh33indexmatrix[,3], lowersamples =  zh33indexmatrix[,4], timepointnames = zh33timepointnames,folder = "ZH33 Mono v Gr",upperbiaslabel = "Mo", lowerbiaslabel = 'Gr')
#make rolling variance of bias plots for ZH33
rollingvarianceofbias(data = zh33alldata, uppersamples = zh33indexmatrixnoNA[,3],lowersamples = zh33indexmatrixnoNA[,4], folder = "ZH33 Mono v Gr", timepoints = zh33timepoints)
rollingvarianceofbias(data = zh33alldata, uppersamples = zh33indexmatrixnoNA[,1],lowersamples = zh33indexmatrixnoNA[,2], folder = "ZH33 T v B", timepoints = zh33timepoints)
rollingvarianceofbias(data = zh33alldata, uppersamples = zh33indexmatrixnoNA[,2],lowersamples = zh33indexmatrixnoNA[,4], folder = "ZH33 B v Gr", timepoints = zh33timepoints)
setwd(workingdir)