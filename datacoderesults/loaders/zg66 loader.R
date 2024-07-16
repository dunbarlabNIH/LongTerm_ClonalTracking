#Rohan Hosuru and Jack Yang (based off Samson Koelle's Code)
#copywrite Cynthia Dunbar Lab
#07-12-23
#This script loads data from animal ZG66

wkdirsave = getwd()
#set this directory to data containing folder
setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/data")

#read main data
zg66 = read.delim("ZG66_combined_20230124_counts.txt")

#read key file containing proper names and sampleIDs
zg66keyfile = read.delim("zg66_key_031824.txt")
# zg66keyfile = read.delim("zg66_key_031824_fig2.txt")
# zg66keyfile = read.delim("zg66_key_edited_NAremoved.txt")

#create vector of timepoints
zg66timepoints = sort(unique(zg66keyfile$MONTH))
zg66timepointnames = as.character(zg66timepoints)

#order data
tomatch = match(zg66keyfile$FILENAME, colnames(zg66))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zg66ordered = zg66[,tomatch]

#add null column to data
zg66ordered = cbind(zg66ordered, list(rep(0, length(zg66[[1]]))))

#index configuration
zg66indexmatrix  = matrix(nrow = length(zg66timepoints), ncol = 4, data = c(37:48,1:12,25:36,13:24))
zg66indexmatrixnoMono = matrix(nrow = length(zg66timepoints), ncol = 3, data = c(37:48,1:12,13:24))
#index configuration to generate Figure 2
# zg66indexmatrix = matrix(nrow = length(zg66timepoints), ncol = 3, data = c(31:45,1:15,16:30))

#create index matrix with 0 column instead of NA
zg66indexmatrixnoNA = zg66indexmatrix
zg66indexmatrixnoNA[is.na(zg66indexmatrix)] = nrow(zg66keyfile) + 1

zg66indexmatrixnoNAnoMono = zg66indexmatrixnoMono
zg66indexmatrixnoNAnoMono[is.na(zg66indexmatrixnoMono)] = nrow(zg66keyfile) + 1

#save data with no threshold applied
zg66alldatanothresh = zg66ordered

#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zg66alldatanothresh = basethresholdinsamples(data=  zg66alldatanothresh, samples = zg66indexmatrixnoNA, thresh = 100)

#apply threshold
zg66alldata = thresholdinsamples(data = zg66alldatanothresh, samples = zg66indexmatrixnoNA, thresh = 2000)

#get names
zg66allnames = zg66keyfile$GIVENNAME
zg66allnames = append(as.vector(zg66allnames), "NULL SAMPLE")

#import GFP levels
zg66GFPraw = read.delim("ZG66_GFP.txt")
rownames(zg66GFPraw) = zg66GFPraw[[1]]
zg66GFP = (zg66GFPraw[c(1:30),c(4,5,2,3)])
colnames(zg66GFP) = c("T","B",'Mono',"Gr")
zg66GFPtimepoints = sort(unique(zg66GFPraw$Months))


setwd(wkdirsave)
