### Boxplot Summary of Percent Change of Fractional Contribution within 
### Granulocyte Lineage

library(tidyverse)

setwd("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/SFigure2/Boxplot Data")

#Import all the data
ZG66_Top10 <- read.delim("Gr_matrix_ZG66_Top10.txt")
ZH19_Top10 <- read.delim("Gr_matrix_ZH19_Top10.txt")
ZH33_Top10 <- read.delim("Gr_matrix_ZH33_Top10.txt")
ZJ31_Top10 <- read.delim("Gr_matrix_ZJ31_Top10.txt")
ZJ38_Top10 <- read.delim("Gr_matrix_ZJ38_Top10.txt")

ZG66_Top100 <- read.delim("Gr_matrix_ZG66_Top100.txt")
ZH19_Top100 <- read.delim("Gr_matrix_ZH19_Top100.txt")
ZH33_Top100 <- read.delim("Gr_matrix_ZH33_Top100.txt")
ZJ31_Top100 <- read.delim("Gr_matrix_ZJ31_Top100.txt")
ZJ38_Top100 <- read.delim("Gr_matrix_ZJ38_Top100.txt")

#Formatting
## Only need one length value as Top 10 and 100 for each animal have same
## number of comparisons
ZG66 <- rep("ZG66", dim(ZG66_Top10)[1])
ZH19 <- rep("ZH19", dim(ZH19_Top10)[1])
ZH33 <- rep("ZH33", dim(ZH33_Top10)[1])
ZJ31 <- rep("ZJ31", dim(ZJ31_Top10)[1])
ZJ38 <- rep("ZJ38", dim(ZJ38_Top10)[1])

## Adding animal information as column
ZG66_Top10 <- data.frame(ZG66_Top10, ZG66)
ZH19_Top10 <- data.frame(ZH19_Top10, ZH19)
ZH33_Top10 <- data.frame(ZH33_Top10, ZH33)
ZJ31_Top10 <- data.frame(ZJ31_Top10, ZJ31)
ZJ38_Top10 <- data.frame(ZJ38_Top10, ZJ38)

colnames(ZG66_Top10)[4] <- "ID"
colnames(ZH19_Top10)[4] <- "ID"
colnames(ZH33_Top10)[4] <- "ID"
colnames(ZJ31_Top10)[4] <- "ID"
colnames(ZJ38_Top10)[4] <- "ID"

ZG66_Top100 <- data.frame(ZG66_Top100, ZG66)
ZH19_Top100 <- data.frame(ZH19_Top100, ZH19)
ZH33_Top100 <- data.frame(ZH33_Top100, ZH33)
ZJ31_Top100 <- data.frame(ZJ31_Top100, ZJ31)
ZJ38_Top100 <- data.frame(ZJ38_Top100, ZJ38)

colnames(ZG66_Top100)[4] <- "ID"
colnames(ZH19_Top100)[4] <- "ID"
colnames(ZH33_Top100)[4] <- "ID"
colnames(ZJ31_Top100)[4] <- "ID"
colnames(ZJ38_Top100)[4] <- "ID"

#Combine into Top 10 and Top 100 tables
All_Top10 <- bind_rows(ZG66_Top10, ZH19_Top10, ZH33_Top10, ZJ31_Top10, ZJ38_Top10)
All_Top10 <- All_Top10[,c(2,4)]
mean(All_Top10$Mean_Percent_Change)
sd(All_Top10$Mean_Percent_Change)

All_Top100 <- bind_rows(ZG66_Top100, ZH19_Top100, ZH33_Top100, ZJ31_Top100, ZJ38_Top100)
All_Top100 <- All_Top100[,c(2,4)]
mean(All_Top100$Mean_Percent_Change)
sd(All_Top100$Mean_Percent_Change)

#Boxplot
top10plot <- ggplot(All_Top10, aes(x = ID, y = Mean_Percent_Change, color = ID)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(5, "Set1"),rep("black",5))) +
  ggtitle("Mean Percent Change For Top 10 Clones Across All Timepoints Per Animal") +
  labs(x = "Animal", y = "Overall Mean Percent Change", color = "Animal") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  theme(axis.text.x = element_text(angle = 45, size = 20, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_y_continuous(breaks = seq(-10,10,2))

top10plot

top100plot <- ggplot(All_Top100, aes(x = ID, y = Mean_Percent_Change, color = ID)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(5, "Set1"),rep("black",5))) +
  ggtitle("Mean Percent Change For Top 100 Clones Across All Timepoints Per Animal") +
  labs(x = "Animal", y = "Overall Mean Percent Change", color = "Animal") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  theme(axis.text.x = element_text(angle = 45, size = 20, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_y_continuous(breaks = seq(-20,10,4))

top100plot
