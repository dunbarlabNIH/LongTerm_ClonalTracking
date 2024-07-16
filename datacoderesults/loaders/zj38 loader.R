#Rohan Hosuru and Jack Yang (based off Samson Koelle's Code)
#copywrite Cynthia Dunbar Lab
#11-14-23
#This script loads data from animal ZJ38

wkdirsave = getwd()
#set this directory to data containing folder
setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/data")

#read main data
zj38 = read.delim("ZJ38_combined_20240620_counts.txt")

#read key file containing proper names and sampleIDs
zj38keyfile = read.delim("ZJ38_key_062024.txt")

#create vector of timepoints
zj38timepoints = sort(unique(zj38keyfile$MONTH))
zj38timepointnames = as.character(zj38timepoints)

#order data
tomatch = match(zj38keyfile$FILENAME, colnames(zj38))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zj38ordered = zj38[,tomatch]

#add null column to data
zj38ordered = cbind(zj38ordered, list(rep(0, length(zj38[[1]]))))

#index configuration
zj38indexmatrix = matrix(nrow = length(zj38timepoints), ncol = 4,data = c(24,NA,25:29,NA,30:31,
                                                                          1,NA,2:4,NA,5,NA,6:7,
                                                                          NA,16:21,NA,22:23,
                                                                          NA,8:14,NA,15))
zj38indexmatrixnoMono = matrix(nrow = length(zj38timepoints), ncol = 3, data = c(24,NA,25:29,NA,30:31,
                                                                                 1,NA,2:4,NA,5,NA,6:7,
                                                                                 NA,8:14,NA,15))
#index configuration for generating Figure 2
# zj38indexmatrix = matrix(nrow = length(zj38timepoints), ncol = 3, data = c(NA,21:24,NA,25:30,1:5,NA,6:8,NA,9:10,11:12,NA,13,NA,14:20))

#this one is for NA removed (bias lines)
# zj38indexmatrix = matrix(nrow = length(zj38timepoints), ncol = 4, data = c(19:24,1:6,13:18,7:12))


#create index matrix with 0 column instead of NA
zj38indexmatrixnoNA = zj38indexmatrix
zj38indexmatrixnoNA[is.na(zj38indexmatrix)] = nrow(zj38keyfile) + 1

zj38indexmatrixnoNAnoMono = zj38indexmatrixnoMono
zj38indexmatrixnoNAnoMono[is.na(zj38indexmatrixnoMono)] = nrow(zj38keyfile) + 1

#save data with no threshold applied
zj38alldatanothresh <- zj38ordered

#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zj38alldatanothresh = basethresholdinsamples(data = zj38alldatanothresh, samples = zj38indexmatrixnoNA, thresh = 100)

#apply threshold
zj38alldata = thresholdinsamples(data = zj38alldatanothresh, samples = zj38indexmatrixnoNA, thresh = 2000)

#get names
zj38allnames = zj38keyfile$GIVENNAME

#import GFP levels
zj38GFPraw = read.delim("ZJ38_GFP.txt")
rownames(zj38GFPraw) = zj38GFPraw[[1]]
zj38GFP = (zj38GFPraw[c(1:13),c(4,5,2,3)])
colnames(zj38GFP) = c("T","B",'Mono',"Gr")
zj38GFPtimepoints = sort(unique(zj38GFPraw$Months))

setwd(wkdirsave)

