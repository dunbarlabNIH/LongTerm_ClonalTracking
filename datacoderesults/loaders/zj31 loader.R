#Rohan Hosuru and Jack Yang (based off Samson Koelle's Code)
#copywrite Cynthia Dunbar Lab
#07-12-23
#This script loads data from animal ZJ31

wkdirsave = getwd()
#set this directory to data containing folder
setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/data")

#read main data
zj31 = read.delim("ZJ31_combined_20220512_counts.txt")

#read key file containing proper names and sampleIDs
zj31keyfile = read.delim("zj31_key_031824.txt")
# zj31keyfile = read.delim("zj31_key_031824_fig2.txt")
# zj31keyfile = read.delim("zj31_key_edited_NAremoved.txt")

#create vector of timepoints
zj31timepoints = sort(unique(zj31keyfile$MONTH))
zj31timepointnames = as.character(zj31timepoints)

#order data
tomatch = match(zj31keyfile$FILENAME, colnames(zj31))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zj31ordered = zj31[,tomatch]

#add null column to data
zj31ordered = cbind(zj31ordered, list(rep(0, length(zj31[[1]]))))

#index configuration
zj31indexmatrix = matrix(nrow = length(zj31timepoints), ncol = 4, data = c(31:40,1:10,21:30,11:20))
zj31indexmatrixnoMono = matrix(nrow = length(zj31timepoints), ncol = 3, data = c(31:40,1:10,11:20))

#index configuration for generating Figure 2
# zj31indexmatrix = matrix(nrow = length(zj31timepoints), ncol = 3, data = c(29:42,1:14,15:28))

#create index matrix with 0 column instead of NA
zj31indexmatrixnoNA = zj31indexmatrix
zj31indexmatrixnoNA[is.na(zj31indexmatrix)] = nrow(zj31keyfile) + 1

zj31indexmatrixnoNAnoMono = zj31indexmatrixnoMono
zj31indexmatrixnoNAnoMono[is.na(zj31indexmatrixnoMono)] = nrow(zj31keyfile) + 1

#save data with no threshold applied
zj31alldatanothresh <- zj31ordered

#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zj31alldatanothresh = basethresholdinsamples(data = zj31alldatanothresh, samples = zj31indexmatrixnoNA, thresh = 100)

#apply threshold
zj31alldata = thresholdinsamples(data = zj31alldatanothresh, samples = zj31indexmatrixnoNA, thresh = 2000)

#get names
zj31allnames = zj31keyfile$GIVENNAME

#import GFP levels
zj31GFPraw = read.delim("ZJ31_GFP.txt")
rownames(zj31GFPraw) = zj31GFPraw[[1]]
zj31GFP = (zj31GFPraw[c(1:20),c(4,5,2,3)])
colnames(zj31GFP) = c("T","B",'Mono',"Gr")
zj31GFPtimepoints = sort(unique(zj31GFPraw$Months))

setwd(wkdirsave)

