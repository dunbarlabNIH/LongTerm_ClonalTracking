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
source("zg66 loader.R")
source("zj31 loader.R")
source("zh19 loader bias.R")
source("zj38 loader bias.R")

#Make Figures
setwd(workingdir)
dirname = paste("Figure S7", Sys.time())
dir.create(dirname)
setwd(dirname)
#make bias tracking plots for all other monkeys
biastracking(data = zg66alldata, uppersamples = zg66indexmatrix[,2],lowersamples =  zg66indexmatrix[,4], timepointnames = zg66timepointnames,folder = "ZG66 B v Gr",upperbiaslabel = "B", lowerbiaslabel = 'Gr')
biastracking(data = zg66alldata, uppersamples = zg66indexmatrix[,1],lowersamples =  zg66indexmatrix[,2], timepointnames = zg66timepointnames,folder = "ZG66 T v B",upperbiaslabel = "T", lowerbiaslabel = 'B')
biastracking(data = zg66alldata, uppersamples = zg66indexmatrix[,3], lowersamples =  zg66indexmatrix[,4], timepointnames = zg66timepointnames,folder = "ZG66 Mono v Gr",upperbiaslabel = "Mo", lowerbiaslabel = 'Gr')

biastracking(data = zh19alldata, uppersamples = zh19indexmatrix[,2],lowersamples =  zh19indexmatrix[,4], timepointnames = zh19timepointnames,folder = "ZH19 B v Gr",upperbiaslabel = "B", lowerbiaslabel = 'Gr')
biastracking(data = zh19alldata, uppersamples = zh19indexmatrix[,1],lowersamples =  zh19indexmatrix[,2], timepointnames = zh19timepointnames,folder = "ZH19 T v B",upperbiaslabel = "T", lowerbiaslabel = 'B')
biastracking(data = zh19alldata, uppersamples = zh19indexmatrix[,3], lowersamples =  zh19indexmatrix[,4], timepointnames = zh19timepointnames,folder = "ZH19 Mono v Gr",upperbiaslabel = "Mo", lowerbiaslabel = 'Gr')

biastracking(data = zj31alldata, uppersamples = zj31indexmatrix[,2],lowersamples =  zj31indexmatrix[,4], timepointnames = zj31timepointnames,folder = "ZJ31 B v Gr",upperbiaslabel = "B", lowerbiaslabel = 'Gr')
biastracking(data = zj31alldata, uppersamples = zj31indexmatrix[,1],lowersamples =  zj31indexmatrix[,2], timepointnames = zj31timepointnames,folder = "ZJ31 T v B",upperbiaslabel = "T", lowerbiaslabel = 'B')
biastracking(data = zj31alldata, uppersamples = zj31indexmatrix[,3], lowersamples =  zj31indexmatrix[,4], timepointnames = zj31timepointnames,folder = "ZJ31 Mono v Gr",upperbiaslabel = "Mo", lowerbiaslabel = 'Gr')
#
biastracking(data = zj38alldata, uppersamples = zj38indexmatrix[,2],lowersamples =  zj38indexmatrix[,4], timepointnames = zj38timepointnames,folder = "ZJ38 B v Gr",upperbiaslabel = "B", lowerbiaslabel = 'Gr')
biastracking(data = zj38alldata, uppersamples = zj38indexmatrix[,1],lowersamples =  zj38indexmatrix[,2], timepointnames = zj38timepointnames,folder = "ZJ38 T v B",upperbiaslabel = "T", lowerbiaslabel = 'B')
biastracking(data = zj38alldata, uppersamples = zj38indexmatrix[,3], lowersamples =  zj38indexmatrix[,4], timepointnames = zj38timepointnames,folder = "ZJ38 Mono v Gr",upperbiaslabel = "Mo", lowerbiaslabel = 'Gr')

#make rolling variance of bias plots for all other monkeys
rollingvarianceofbias(data = zg66alldata, uppersamples = zg66indexmatrixnoNA[,3],lowersamples = zg66indexmatrixnoNA[,4], folder = "ZG66 Mono v Gr", timepoints = zg66timepoints)
rollingvarianceofbias(data = zg66alldata, uppersamples = zg66indexmatrixnoNA[,1],lowersamples = zg66indexmatrixnoNA[,2], folder = "ZG66 T v B", timepoints = zg66timepoints)
rollingvarianceofbias(data = zg66alldata, uppersamples = zg66indexmatrixnoNA[,2],lowersamples = zg66indexmatrixnoNA[,4], folder = "ZG66 B v Gr", timepoints = zg66timepoints)

rollingvarianceofbias(data = zh19alldata, uppersamples = zh19indexmatrixnoNA[,3],lowersamples = zh19indexmatrixnoNA[,4], folder = "ZH19 Mono v Gr", timepoints = zh19timepoints)
rollingvarianceofbias(data = zh19alldata, uppersamples = zh19indexmatrixnoNA[,1],lowersamples = zh19indexmatrixnoNA[,2], folder = "ZH19 T v B", timepoints = zh19timepoints)
rollingvarianceofbias(data = zh19alldata, uppersamples = zh19indexmatrixnoNA[,2],lowersamples = zh19indexmatrixnoNA[,4], folder = "ZH19 B v Gr", timepoints = zh19timepoints)

rollingvarianceofbias(data = zj31alldata, uppersamples = zj31indexmatrixnoNA[,3],lowersamples = zj31indexmatrixnoNA[,4], folder = "ZJ31 Mono v Gr", timepoints = zj31timepoints)
rollingvarianceofbias(data = zj31alldata, uppersamples = zj31indexmatrixnoNA[,1],lowersamples = zj31indexmatrixnoNA[,2], folder = "ZJ31 T v B", timepoints = zj31timepoints)
rollingvarianceofbias(data = zj31alldata, uppersamples = zj31indexmatrixnoNA[,2],lowersamples = zj31indexmatrixnoNA[,4], folder = "ZJ31 B v Gr", timepoints = zj31timepoints)

rollingvarianceofbias(data = zj38alldata, uppersamples = zj38indexmatrixnoNA[,3],lowersamples = zj38indexmatrixnoNA[,4], folder = "ZJ38 Mono v Gr", timepoints = zj38timepoints)
rollingvarianceofbias(data = zj38alldata, uppersamples = zj38indexmatrixnoNA[,1],lowersamples = zj38indexmatrixnoNA[,2], folder = "ZJ38 T v B", timepoints = zj38timepoints)
rollingvarianceofbias(data = zj38alldata, uppersamples = zj38indexmatrixnoNA[,2],lowersamples = zj38indexmatrixnoNA[,4], folder = "ZJ38 B v Gr", timepoints = zj38timepoints)

setwd(workingdir)