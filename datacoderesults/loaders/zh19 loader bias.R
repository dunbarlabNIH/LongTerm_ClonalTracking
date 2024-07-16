#Rohan Hosuru and Jack Yang (based off Samson Koelle's Code)
#copywrite Cynthia Dunbar Lab
#07-12-23
#This script loads data from animal ZH19

wkdirsave = getwd()
#set this directory to data containing folder
setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/data")

#read main data
zh19 = read.delim("ZH19_combined_20240620_counts.txt")

#read key file containing proper names and sampleIDs
zh19keyfile = read.delim("ZH19_key_06_28_24_NAremoved.txt")

#create vector of timepoints
zh19timepoints = sort(unique(zh19keyfile$MONTH))
zh19timepointnames = as.character(zh19timepoints)

#order data
tomatch = match(zh19keyfile$FILENAME, colnames(zh19))
fillervalue  = min(setdiff(seq(from = 1, to = 100, by = 1), tomatch[!is.na(tomatch)]))
tomatch[is.na(tomatch)] = fillervalue
zh19ordered = zh19[,tomatch]

#add null column to data
zh19ordered = cbind(zh19ordered, list(rep(0, length(zh19[[1]]))))

#index configuration - manually update nrow and the row numbers in the order of T, B, Mono, Gr to reflect key file
zh19indexmatrix = matrix(nrow = length(zh19timepoints), ncol = 4, data = c(34:44,
                                                                           1:11,
                                                                           23:33,
                                                                           12:22))

#index configuration to generate Figure 2
# zh19indexmatrix = matrix(nrow = length(zh19timepoints), ncol = 3, data = c(28:41,1:14,15:26,NA,27))

# this one is for no NA (for bias plots):
# zh19indexmatrix = matrix(nrow = length(zh19timepoints), ncol = 4, data = c(31:40,1:10,21:30,11:20))

#create index matrix with 0 column instead of NA
zh19indexmatrixnoNA = zh19indexmatrix
zh19indexmatrixnoNA[is.na(zh19indexmatrix)] = nrow(zh19keyfile) + 1

#save data with no threshold applied
zh19alldatanothresh <- zh19ordered

#apply base threshold which replicates preselection of samples in fastq pipeline (prepublication only)
zh19alldatanothresh = basethresholdinsamples(data = zh19alldatanothresh, samples = zh19indexmatrixnoNA, thresh = 100)

#apply threshold
zh19alldata = thresholdinsamples(data = zh19alldatanothresh, samples = zh19indexmatrixnoNA, thresh = 2000)

#create vector of sample names
zh19allnames = zh19keyfile$GIVENNAME

#get names
zh19allnames = append(as.vector(zh19allnames), "NULL SAMPLE")

#import GFP levels
zh19GFPraw = read.delim("ZH19_GFP.txt")
rownames(zh19GFPraw) = zh19GFPraw[[1]]
zh19GFP = (zh19GFPraw[c(1:22),c(4,5,2,3)])
colnames(zh19GFP) = c("T","B",'Mono',"Gr")
zh19GFPtimepoints = sort(unique(zh19GFPraw$Months))


setwd(wkdirsave)