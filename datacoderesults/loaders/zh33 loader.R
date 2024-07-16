#Rohan Hosuru and Jack Yang (based off Samson Koelle's Code)
#copywrite Cynthia Dunbar Lab
#07-12-23
#This script loads data from animal ZH33

wkdirsave = getwd()
#set this directory to data containing folder
setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/data")

#read main data
zh33 = read.delim("ZH33_combined_20231019_counts.txt")

#read key file containing proper names and sampleIDs
zh33keyfile = read.delim("ZH33_key_031824.txt")
# zh33keyfile = read.delim("ZH33_key_031824_fig2.txt")


#create vector of timepoints
zh33timepoints = sort(unique(zh33keyfile$MONTH))
zh33timepointnames = as.character(zh33timepoints)

#order data
tomatch = match(zh33keyfile$FILENAME, colnames(zh33))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zh33ordered = zh33[,tomatch]

#add null column to data
zh33ordered = cbind(zh33ordered, list(rep(0, length(zh33[[1]]))))

#index configuration - manually update nrow and the row numbers in the order of T, B, Mono, Gr to reflect key file
zh33indexmatrix = matrix(nrow = length(zh33timepoints), ncol = 4, data = c(34:44,1:11,23:33,12:22))
zh33indexmatrixnoMono = matrix(nrow = length(zh33timepoints), ncol = 3, data = c(34:44,1:11,12:22))

#index configuration to generate Figure 2
# zh33indexmatrix = matrix(nrow = length(zh33timepoints), ncol = 3, data = c(31:45,1:15,16:30))

#create index matrix with 0 column instead of NA
zh33indexmatrixnoNA = zh33indexmatrix
zh33indexmatrixnoNA[is.na(zh33indexmatrix)] = nrow(zh33keyfile) + 1

zh33indexmatrixnoNAnoMono = zh33indexmatrixnoMono
zh33indexmatrixnoNAnoMono[is.na(zh33indexmatrixnoMono)] = nrow(zh33keyfile) + 1

#save data with no threshold applied
zh33alldatanothresh <- zh33ordered

#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zh33alldatanothresh = basethresholdinsamples(data = zh33alldatanothresh, samples = zh33indexmatrixnoNA, thresh = 100)

#apply threshold
zh33alldata = thresholdinsamples(data = zh33alldatanothresh, samples = zh33indexmatrixnoNA, thresh = 2000)

#create vector of sample names
zh33allnames = zh33keyfile$GIVENNAME

#get names
zh33allnames = append(as.vector(zh33allnames), "NULL SAMPLE")

#import GFP levels
zh33GFPraw = read.delim("ZH33_GFP.txt")
rownames(zh33GFPraw) = zh33GFPraw[[1]]
zh33GFP = (zh33GFPraw[c(1:30),c(4,5,2,3)])
colnames(zh33GFP) = c("T","B",'Mono',"Gr")
zh33GFPtimepoints = sort(unique(zh33GFPraw$Months))

setwd(wkdirsave)