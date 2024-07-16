library(tidyverse)
library(dplyr)
library(DelayedArray)
library(ggplot2)
library(gplots)

#test conducted with ZH33 Top 10/100 Clones Data
##How to generate for other monkeys
###Input appropriate counts and metadata into the barcodeTrackR App
###Use sample names from the key file for Gr
###Toggle number of top clones
###Download the heatmap data

##Set working directory
#read in data
nonMatrixHeatData <- read.delim("ZJ31_HeatmapData_Top100Clones.txt")

#You may notice the name of the previous variable is non matrix
#This is because the heatmap data is outputted in the long format
#For our analysis, we want to transform this into a matrix (wide format) of 
#fractional abundance values

#We want x rows for the number of barcodes and y columns for the number of samples
#Unique maintains order of barcode and samples
nRows <- length(unique(nonMatrixHeatData[,1]))
nCols <- length(unique(nonMatrixHeatData[,2]))

#The cool thing about making a matrix in R is one of the input values is the data
#we are interested. We specify the the nRows and nCols, and now we can specify
#how we want to fill in the matrix - by row or by column. For us, we will do
#by row (indicated by )
matrixHeatData <- matrix(data = nonMatrixHeatData[,3], nrow = nRows, ncol = nCols,
                         byrow = T)

#Assign column names
colnames(matrixHeatData) <- unique(nonMatrixHeatData[,2])

#calculate percent change - function for easy repeat use
#change = (final value - initial value) / initial value
#percentChange <- matrix of zeroes same row number as normalized, one less column
#for loop goes through all columns except last, goes through all rows
#if statement checks if the denominator will be zero (undefined), if so, then
#will give the percent change a value of NA and next jumps to the next loop iteration

calc_perChange <- function(normalized){
  num_Col <- ncol(normalized)
  num_Row <- nrow(normalized)
  percentChange <- matrix(rep(0, (num_Col-1) * num_Row), ncol = num_Col - 1,
                          nrow = num_Row)
  for (i in 1:(num_Col-1)){
    for (j in 1:(num_Row)){
      if(normalized[j,i] == 0){
        percentChange[j,i] <- NA
        next
      }
      percentChange[j,i] <- 100 * ((normalized[j,i+1] - normalized[j,i]) / normalized[j,i])
    }
  }
  return(percentChange)
}

percentChangeMatrix <- calc_perChange(matrixHeatData)

#sample names percent change is calculated between these samples
extract_compareNames <- function(normalized, percentChange){
  numSamples <- dim(normalized)[2]
  return(paste(colnames(normalized)[1:(numSamples-1)],colnames(normalized)[2:numSamples]))
}

Gr_compareNames <- extract_compareNames(matrixHeatData, percentChangeMatrix)

#Create txt file of combined sample name, we will use this to create
#xAxisNames as percentChange colnames are too chaotic
write.table(Gr_compareNames, file = "Gr_test_compareNames.txt")

#calculate the mean percent change per sample (colAverage)
#percentChange input, 2 for by column, mean, na.rm = T (drops NA in calculation)
#combine into one dataset
#Set xAxisNames, if using combined sample names, too large for figure
Gr_colAverage <- apply(percentChangeMatrix, 2, mean, na.rm = TRUE)

Gr_colSD <- apply(percentChangeMatrix, 2, sd, na.rm = TRUE)

#write in the refined names we want for each comparison
Gr_xAxisNames <- read.delim("ZJ31_Gr_GraphNames.txt", header = F)
Gr_xAxisNames <- dplyr::pull(Gr_xAxisNames, 3)
Gr_matrix <- cbind(Gr_xAxisNames, as.numeric(Gr_colAverage), as.numeric(Gr_colSD))

combined_df <- data.frame(Gr_matrix)
combined_df[,2:3] <- as.numeric(unlist(combined_df[,2:3]))
colnames(combined_df) <- c("Timepoint","Mean_Percent_Change", "SD")

write.table(combined_df,"Gr_matrix_ZJ31_Top100.txt",sep="\t",row.names=FALSE)

rightOrder <- c("6m vs.12m","12m vs.20m","20m vs.28m","28m vs.33m","33m vs.42m",
                "42m vs.44m","44m vs.52.5m","52.5m vs.63m","63m vs.75m")
combined_df$Timepoint <- factor(combined_df$Timepoint, levels = rightOrder)

#color by cell type
perPlot <- ggplot(combined_df, aes(x = factor(Timepoint, level = rightOrder), y = Mean_Percent_Change, group = 1)) +
                                  labs(x = "Timepoint Comparison", y = "Mean Percent Change") +
                                  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
                                  geom_line() +
                                  geom_point(size = 3, shape = 15) +
                                  theme_classic() +
                                  ggtitle("ZJ31 Top 100 Clones") +
                                  theme(axis.text.x = element_text(angle = 60, size = 20, hjust = 1),
                                        axis.text.y = element_text(size = 16),
                                        axis.title.x = element_text(size = 18),
                                        axis.title.y = element_text(size = 18),
                                        plot.title = element_text(size = 16, hjust = 0.5)) +
                                  scale_y_continuous(breaks = seq(-100,100,10))
  

perPlot <- perPlot + geom_errorbar(aes(ymin = Mean_Percent_Change - SD,
                                       ymax = Mean_Percent_Change + SD),
                                   width = 0.3, linewidth = 0.5)

perPlot
