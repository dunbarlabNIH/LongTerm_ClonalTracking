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

#Make Figures
setwd(workingdir)
#Create images for Figure 2 - MAKE SURE TO USE DATA/LOADERS WITHOUT MONO
dirname = paste("Figure 2 Images", Sys.time())
dir.create(dirname)
setwd(dirname)
#Compare Shannon Diversities
diversitytrack(data = zh33alldata, indexmatrix = zh33indexmatrixnoMono, timepoints = zh33timepoints, folder = "ZH33")
diversitytrack(data = zg66alldata, indexmatrix = zg66indexmatrixnoMono, timepoints = zg66timepoints, folder = "ZG66")
diversitytrack(data = zh19alldata, indexmatrix = zh19indexmatrixnoMono, timepoints = zh19timepoints, folder = "ZH19")
diversitytrack(data = zj31alldata, indexmatrix = zj31indexmatrixnoMono, timepoints = zj31timepoints, folder = "ZJ31")
diversitytrack(data = zj38alldata, indexmatrix = zj38indexmatrixnoMono, timepoints = zj38timepoints, folder = "ZJ38")

#GFP tracking
GFPtrack(data = zh33GFP, timepoints = zh33GFPtimepoints, folder = "ZH33")
GFPtrack(data = zg66GFP, timepoints = zg66GFPtimepoints, folder = "ZG66")
GFPtrack(data = zh19GFP, timepoints = zh19GFPtimepoints, folder = "ZH19")
GFPtrack(data = zj31GFP, timepoints = zj31GFPtimepoints, folder = "ZJ31")
GFPtrack(data = zj38GFP, timepoints = zj38GFPtimepoints, folder = "ZJ38")

#Compare clone numbers
clonenumbertrack(data = zh33alldata, indexmatrixnoNA = zh33indexmatrixnoNAnoMono, timepoints = zh33timepoints, folder = "ZH33")
clonenumbertrack(data = zg66alldata, indexmatrixnoNA = zg66indexmatrixnoNAnoMono, timepoints = zg66timepoints, folder = "ZG66")
clonenumbertrack(data = zh19alldata, indexmatrixnoNA = zh19indexmatrixnoNAnoMono, timepoints = zh19timepoints, folder = "ZH19")
clonenumbertrack(data = zj31alldata, indexmatrixnoNA = zj31indexmatrixnoNAnoMono, timepoints = zj31timepoints, folder = "ZJ31")
clonenumbertrack(data = zj38alldata, indexmatrixnoNA = zj38indexmatrixnoNAnoMono, timepoints = zj38timepoints, folder = "ZJ38")
setwd(workingdir)