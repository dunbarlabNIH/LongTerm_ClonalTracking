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
source("zj31 loader.R")
source("zh19 loader.R")
source("zj38 loader.R")

#Create Images for Figure 3
setwd(workingdir)
dirname = paste("Figure 3 Images", Sys.time())
dir.create(dirname)
setwd(dirname)
#make top clone heatmaps for Gr
BCLH(data = zh33alldata[,as.vector(zh33indexmatrixnoNA[,4])], names = zh33timepoints,n_clones = 100, folder = "ZH33 Gr")
BCLH(data = zg66alldata[,as.vector(zg66indexmatrixnoNA[,4])], names = zg66timepoints, n_clones = 100, folder = "ZG66 Gr")
BCLH(data = zh19alldata[,as.vector(zh19indexmatrixnoNA[,4])], names = zh19timepoints,n_clones = 100, folder = "ZH19 Gr")
BCLH(data = zj31alldata[,as.vector(zj31indexmatrixnoNA[,4])], names = zj31timepoints, n_clones = 100, folder = "ZJ31 Gr")
BCLH(data = zj38alldata[,as.vector(zj38indexmatrixnoNA[,4])], names = zj38timepoints, n_clones = 100, folder = "ZJ38 Gr")
setwd(workingdir)