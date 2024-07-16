#load up the necessary packages and libraries and functions
installed.packages("barcodetrackR")
installed.packages("magrittr")
installed.packages("SummarizedExperiment")
library(dplyr)
library(barcodetrackR)
library(rlang)
library(ggplot2)
library(RColorBrewer)
require("magrittr")
require("barcodetrackR")
require("SummarizedExperiment")
#we adjusted the barcodetrackR package and so here is the customized function that we used to calculate the cumulative counts
clonal_count_v3<-source("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/Figure3/Figure3F_to_H/Figure3F_functions/clonal_counts_v3.R") 
clonal_count_v3<-clonal_count_v3$value
#we adjusted the barcodetrackR package to accomodate for the low clonal contributions from later timepoints in the animal: 11021142 - so we read in a customized barcode ggheatmap function for that animal's analysis 
barcode_ggheatmap_first<-source("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/Figure3/Figure3F_to_H/Figure3F_functions/barcode_ggheatmap_42.R")
barcode_ggheatmap_first<-barcode_ggheatmap_first$value
autocorrelation.plot<-source("/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/Figure3/Figure3F_to_H/Figure3F_functions/autocorrelation_plot.R")
autocorrelation.plot<-autocorrelation.plot$value
########FIGURE 2B##############
#set up autocorrelation plot

#load data for comparison animals (TBI animals)
my_path<-"/Users/hosururv/Desktop/BJH_LongTermFollowUp_R_code_ByFigure/datacoderesults/Figure3/Data"

#ZG66
zg66.counts <- read.delim(file.path(my_path,"ZG66_combined_20230124_counts.txt"), header=T, row.names=1, sep="\t")
zg66.grans.samples <- c('ZG666hmCD33pGFPpGr.fastq','zg6612mGrDNA.fastq','zg6622mGFPpGrCD807.fastq','zg66_27m_GRANS_ficoll_CD85_R39_CTATAC_L007_R1_001.fastq','zg66_36m_Gr_CD33p_GFPp_20151214_pp_50trim_CD117_R6.fastq','pp_50trim_ZG66_47M__27102016__Gr_ficolled_R9_CD139.fastq','zg66_start_20170223_Gr_sampled_DE2400_YD23_N716_S85_L008_R1_001.fastq','ZG66_Mono_Gr_ficolled_sampled_LT2763_CD150_i511_S70_L008_R1_001.fastq','ZG66_20180724_Grans_CD160_sampled_CD160_4_S4_L001_R1_001.fastq','ZG66_82hm_20191024_Gr__PB_CD181_sampled_I522_S18_L008_R1_001.fastq','ZG66_106m_2021928_Gr_PB__CD211_sampled_8bp_i502_S35_L002_R1_001.fastq','ZG66_112m_202245_Gr_PB_NA_CD218_sampled_DA_8140_66_S66_L002_R1_001.fastq')
zg66.t.samples <- c('ZG666hmT.fastq','zg6612mT.fastq','zg6622mTCD803.fastq','zg66_27m_T_CD85_R34_CATGGC_L007_R1_001.fastq','zg66_36m_T_20151214_pp_50trim_CD117_R1.fastq','pp_50trim_ZG66_47M__27102016__T_R12_CD139.fastq','zg66_start_20170223_T_sampled_DE2400_YD23_N711_S81_L008_R1_001.fastq','ZG66_T_sampled_LT2763_CD150_i505_S66_L008_R1_001.fastq','ZG66_20180724_T_CD160_sampled_CD160_8_S8_L001_R1_001.fastq','ZG66_82hm_20191024_T__PB_CD181_sampled_I523_S19_L008_R1_001.fastq','ZG66_106m_2021928_T_PB__CD211_sampled_i512_32_S26_L002_R1_001.fastq','ZG66_112m_202245_T_PB_NA_CD218_sampled_DA_8140_68_S68_L002_R1_001.fastq')
zg66.b.samples <- c('ZG666hmB.fastq','zg6612mB.fastq','zg6622mBCD804.fastq','zg66_27m_B_CD85_R35_CATTTT_L007_R1_001.fastq','zg66_36m_B_resequenced_20151214_pp_50trim_CD122_R2.fastq','pp_50trim_ZG66_47M__27102016__B_R6_CD139.fastq','zg66_start_20170223_B_sampled_DE2400_YD23_N712_S82_L008_R1_001.fastq','ZG66_B_sampled_LT2763_CD150_i506_S67_L008_R1_001.fastq','ZG66_20180724_B_CD160_sampled_CD160_9_S9_L001_R1_001.fastq','ZG66_82hm_20191024_B__PB_CD181_sampled_I524_S20_L008_R1_001.fastq','ZG66_106m_2021928_B_PB__CD211_sampled_i512_33_S27_L002_R1_001.fastq','ZG66_112m_202245_B_PB_NA_CD218_sampled_DA_8140_69_S69_L002_R1_001.fastq')
zg66.mono.samples <- c('ZG666hmMono.fastq','zg6612mMono.fastq','zg6622mMonoCD805.fastq','zg66_27m_MONO_CD85_R36_CCAACA_L007_R1_001.fastq','zg66_36m_Mono_20151214_pp_50trim_CD117_R3.fastq','pp_50trim_ZG66_47M__27102016__Mono__R21_CD139.fastq','zg66_start_20170223_Mono_sampled_DE2400_YD23_N714_S83_L008_R1_001.fastq','ZG66_Mono_sampled_LT2763_CD150_i507_S68_L008_R1_001.fastq','ZG66_20180724_Mono_CD160_sampled_CD160_10_S10_L001_R1_001.fastq','ZG66_82hm_20191024_Mono__PB_CD181_sampled_I526_S21_L008_R1_001.fastq','ZG66_106m_2021928_Mono_PB__CD211_sampled_i512_35_S28_L002_R1_001.fastq','ZG66_112m_202245_Mono_PB_NA_CD218_sampled_DA_8140_70_S70_L002_R1_001.fastq')
zg66.time <- NULL
zg66.time$grans=c(6.5,12,22,27,36,47,50.5,55.5,67.5,82.5,106,112)
zg66.time$t=c(6.5,12,22,27,36,47,50.5,55.5,67.5,82.5,106,112)
zg66.time$b=c(6.5,12,22,27,36,47,50.5,55.5,67.5,82.5,106,112)
zg66.time$mono=c(6.5,12,22,27,36,47,50.5,55.5,67.5,82.5,106,112)


#ZH33
zh33.counts <- read.delim(file.path(my_path,"ZH33_combined_20231019_counts.txt"), header=T, row.names=1, sep="\t")
#Grans
zh33.grans.samples <- c('zh336hmGr.fastq','zh3312mGr.fastq','zh3321mGFPpGr.fastq','zh33093014GFPpGr.fastq','zh33_20150806_PB_Gr_Ficoll_CD_104_R12_CTTGTA_L005_R1_001.fastq','zh33_43m_Gr_Ficoll_reextract_20160108_pp_50trim_CD119_R17.fastq','ZH33_11012016_Gr_CD157_sampled_CD157_i506_S6_L007_R1_001.fastq','zh33_67m_Gr_sampled_CD156_i516_S13_L008_R1_001.fastq','ZH33_79m_2019_01_09_Grans__FX52_sampled_FX52_i505_S5_L001_R1_001.fastq','ZH33_88hm_20191017_Gr___CD181_sampled_I516_S13_L008_R1_001.fastq','ZH33_118m_20220204_GRAN_PB__CD221_sampled_new12bp_5_S5_R1_001.fastq')
zh33.t.samples <- c('zh336hmT.fastq','zh3312mT.fastq','zh3321mT.fastq','zh33T093014.fastq','zh33_20150806_PB_T_CD_104_R4_TGACCA_L005_R1_001.fastq','zh33_43m_T_20160108_pp_50trim_CD116_R10.fastq','ZH33_11012016_T_CD157_sampled_CD157_i515_S12_L007_R1_001.fastq','zh33_67m_T_sampled_CD156_i507_S7_L008_R1_001.fastq','ZH33_79m_2019_01_09_T__FX52_sampled_FX52_i506_S6_L001_R1_001.fastq','ZH33_88hm_20191017_CD3pos_T___CD181_sampled_I511_S9_L008_R1_001.fastq','ZH33_118m_2023922_T_1_Necropsy_NA_CD224_sampled_CD224_1_S87_L002_R1_001.fastq')
zh33.b.samples <- c('zh336hmB.fastq','zh3312mB.fastq','zh3321mB.fastq','zh33B093014.fastq','zh33_20150806_PB_B_CD_104_R1_ATCACG_L005_R1_001.fastq','zh33_43m_B_20160108_pp_50trim_CD118_R11.fastq','ZH33_11012016_B_CD157_sampled_CD157_i516_S13_L007_R1_001.fastq','zh33_67m_B_sampled_CD156_i510_S8_L008_R1_001.fastq','ZH33_79m_2019_01_09_B__FX52_sampled_FX52_i507_S7_L001_R1_001.fastq','ZH33_88hm_20191017_CD20pos_B___CD181_sampled_I512_S10_L008_R1_001.fastq','ZH33_118m_2023922_B_1_Necropsy_NA_CD224_sampled_CD224_2_S88_L002_R1_001.fastq')
zh33.mono.samples <- c('zh336hmMono.fastq','zh3312mMono.fastq','zh3321mMono.fastq','zh33Mono093014.fastq','zh33_20140806_PB_Mono_CD_104_R2_CGATGT_L005_R1_001.fastq','zh33_43m_Mono_20160108_pp_50trim_CD117_R12.fastq','ZH33_11012016_Mono_CD157_sampled_CD157_i518_S14_L007_R1_001.fastq','zh33__67m_Mono_sampled_CD156_i511_S9_L008_R1_001.fastq','ZH33_79m_2019_01_09_Mono__FX52_sampled_FX52_i510_S8_L001_R1_001.fastq','ZH33_88hm_20191017_CD14pos_bulk_mono___CD181_sampled_I514_S11_L008_R1_001.fastq','ZH33_118m_2022024_Mono_PB__CD221_sampled_new12bp_55_S54_R1_001.fastq')
zh33.time <- NULL
zh33.time$grans=c(6.5,12,21,28,38,43,53,67,79,88.5,118)
zh33.time$t=c(6.5,12,21,28,38,43,53,67,79,88.5,118)
zh33.time$b=c(6.5,12,21,28,38,43,53,67,79,88.5,118)
zh33.time$mono=c(6.5,12,21,28,38,43,53,67,79,88.5,118)

#ZH19
zh19.counts <- read.delim(file.path(my_path,"ZH19_combined_20240620_counts.txt"), header=T, row.names=1, sep="\t")
#Grans
zh19.grans.samples <- c('zh197mGr.fastq','zh1910mCD33pGFPp.fastq','zh19GFPpCD33p031615.fastq','pp_50trim_zh19_28m_Gr_DiegoDemultiplex.fastq','zh19_36m_Gr_Ficoll_20160509_pp_50trim_CD127_R1.fastq','zh19_44hm_PB_Gr_20170123_sampled_FX27_N706_S5_L008_R1_001.fastq','zh19_09272017_Gr_sampled_CD152_i511_S31_L008_R1_001.fastq','ZH19_62hm_20180725_Grans_CD160_sampled_CD160_1_S1_L001_R1_001.fastq','ZH19_78m_20191029_Gr_PB__CD183_sampled_CD183_I516_S37_L008_R1_001.fastq','ZH19_108m_2022531_Gr_LT__CD226_sampled_CD226_91_S91_L002_R1_001.fastq','ZH19_131m_202448_Gr_LT__CD226_sampled_CD226_92_S92_L002_R1_001.fastq')
zh19.t.samples <- c('zh197mTcd15.fastq','zh1910mT.fastq','zh19T031615.fastq.fastq','pp_50trim_zh19_28m_T_DiegoDemultiplex.fastq','zh19_36m_T_20160509_pp_50trim_CD127_R2.fastq','zh19_44hm_PB_T_20170123_sampled_FX27_N701_S1_L008_R1_001.fastq','zh19_09272017_CD3p_sampled_CD152_i515_S34_L008_R1_001.fastq','ZH19_62hm_20180725_T_CD161_sampled_CD161_12_S43_L002_R1_001.fastq','ZH19_78m_20191029_T_PB__CD183_sampled_CD183_I518_S38_L008_R1_001.fastq','ZH19_114m_20220912_T_PB__CD221_sampled_new12bp_44_S43_R1_001.fastq','ZH19_121m_20230531_T_PB__CD221_sampled_new12bp_49_S48_R1_001.fastq','ZH19_131m_202449_T_LT__CD226_sampled_CD226_81_S81_L002_R1_001.fastq')
zh19.b.samples <- c('zh197mB.fastq','zh1910mB.fastq','zh19B031615.fastq','pp_50trim_zh19_28m_B_CD110_DiegoDemultiplex.fastq','zh19_36m_B_20160509_pp_50trim_CD127_R3.fastq','zh19_44hm_PB_B_20170123_sampled_FX27_N702_S2_L008_R1_001.fastq','zh19_09272017_CD20p_sampled_CD152_i516_S35_L008_R1_001.fastq','ZH19_62hm_20180725_B_CD161_sampled_CD161_14_S44_L002_R1_001.fastq','ZH19_78m_20191029_B_repeat_PB__CD183_sampled_CD183_I514_S35_L008_R1_001.fastq','ZH19_114m_20220912_CD20_PB__CD221_sampled_new12bp_48_S47_R1_001.fastq','ZH19_121m_20230531_B_PB__CD221_sampled_new12bp_50_S49_R1_001.fastq','ZH19_131m_202449_B_LT__CD226_sampled_CD226_97_S97_L002_R1_001.fastq')
zh19.mono.samples <- c('zh197mMono.fastq','zh1910mMono.fastq','zh19Mono031615.fastq','pp_50trim_zh19_28m_Mono_DiegoDemultiplex.fastq','zh19_36m_Mono_CD14_20160509_pp_50trim_CD127_R5.fastq','zh19_44hm_PB_Mono_20170123_sampled_FX27_N703_S3_L008_R1_001.fastq','ZH19_20170927_CD14p_sampled_CD152_i521_S39_L008_R1_001.fastq','ZH19_62hm_20180725_Mono_CD160_sampled_CD160_23_S17_L001_R1_001.fastq','ZH19_78m_20191029_Mono_PB__CD183_sampled_CD183_I520_S40_L008_R1_001.fastq','ZH19_114m_20220912_Mono_PB__CD221_sampled_new12bp_45_S44_R1_001.fastq','ZH19_121m_20230531_Mono_PB__CD221_sampled_new12bp_53_S52_R1_001.fastq','ZH19_131m_202449_Mono_LT__CD226_sampled_CD226_82_S82_L002_R1_001.fastq')
zh19.time <- NULL
zh19.time$grans=c(7,10,22,28,36,44.5,52.5,62.5,78,108,131)
zh19.time$t=c(7,10,22,28,36,44.5,52.5,62.5,78,114,121,131)
zh19.time$b=c(7,10,22,28,36,44.5,52.5,62.5,78,114,121,131)
zh19.time$mono=c(7,10,22,28,36,44.5,52.5,62.5,78,114,121,131)


#ZJ31
zj31.counts <- read.delim(file.path(my_path,"ZJ31_combined_20220512_counts.txt"), header=T, row.names=1, sep="\t")
#Grans
zj31.grans.samples <- c('zj316mGr28cycles.fastq','zj31_12m_Gr_CD33p_CD_103_R2_CGATGT_L004_R1_001.fastq','zj31_20m_Gr_20160328_pp_50trim_CD124_R1.fastq','zj31_28m_Gr_Ficoll_20161130_sampled_DE2242_CD143_N706_S6_L008_R1_001.fastq','zj31_33m_Grans_F_twice_20170424_sampled_FX36_N701_S1_L008_R1_001.fastq','ZJ31_42m_2018_01_23_Grans__PB_FX48_sampled_FX48_i502_S26_L005_R1_001.fastq','ZJ31_44m_2018_03_19_Grans__PB_FX48_sampled_FX48_i515_S36_L005_R1_001.fastq','ZJ31__20181211_Grans___CD163_sampled_CD163_i529_S120_L006_R1_001.fastq','ZJ31_63m_20191021_Gr___CD181_sampled_I518_S14_L008_R1_001.fastq','ZJ31_75m_20201020_Gr__PB_CD193_sampled_i512_06_S6_L001_R1_001.fastq')
zj31.t.samples <- c('zj316mPBT28cyclesCD75R.fastq','zj31_12m_T_CD_103_R3_TTAGGC_L004_R1_001.fastq','zj31_20m_T_CD3p_20160328_pp_50trim_CD124_R2.fastq','zj31_28m_T_20161130_D701_sampled_DE2242_CD143_D701_S15_L008_R1_001.fastq','zj31_33m_T_20170424_sampled_FX36_N703_S3_L008_R1_001.fastq','ZJ31_42m_2018_01_23_T__PB_FX47_sampled_FX47_i526_S21_L004_R1_001.fastq','ZJ31_44m_2018_03_19_T__PB_FX48_sampled_FX48_i510_S32_L005_R1_001.fastq','ZJ31__20181211_T___CD164_sampled_CD164_i516_S61_L004_R1_001.fastq','ZJ31_63m_20191021_CD3pos_T___CD181_sampled_I505_S5_L008_R1_001.fastq','ZJ31_75m_20201020_T__PB_CD193_sampled_i512_01_S1_L001_R1_001.fastq')
zj31.b.samples <- c('zj316mB28cycles.fastq','zj31_12m_B_CD_103_R4_TGACCA_L004_R1_001.fastq','zj31_20m_B_20160328_pp_50trim_CD124_R3.fastq','zj31_28m_B_20161130_sampled_DE2242_CD143_N702_S2_L008_R1_001.fastq','zj31_33m_B_20170424_sampled_FX36_N702_S2_L008_R1_001.fastq','ZJ31_42m_2018_01_23_B__PB_FX47_sampled_FX47_i527_S22_L004_R1_001.fastq','ZJ31_44m_2018_03_19_B__PB_FX48_sampled_FX48_i511_S33_L005_R1_001.fastq','ZJ31__20181211_B___CD163_sampled_CD163_i527_S118_L006_R1_001.fastq','ZJ31_63m_20191021_CD20pos_B___CD181_sampled_I506_S6_L008_R1_001.fastq','ZJ31_75m_20201020_B__PB_CD193_sampled_i512_02_S2_L001_R1_001.fastq')
zj31.mono.samples <- c('zj316mMono28cyclesCD76R.fastq','zj31_12m_Mono_CD_103_R5_ACAGTG_L004_R1_001.fastq','zj31_20m_Mono_20160328_pp_50trim_CD124_R4.fastq','zj31_28m_Mono_20161130_sampled_DE2242_CD143_N703_S3_L008_R1_001.fastq','zj31_33m_20170424_Mono_trimmed_sampled_CD146_R24_S16_L002_R1_001.fastq','ZJ31_42m_2018_01_23_Mono__PB_FX47_sampled_FX47_i528_S23_L004_R1_001.fastq','ZJ31_44m_2018_03_19_Mono__PB_FX48_sampled_FX48_i512_S34_L005_R1_001.fastq','ZJ31__20181211_Mono___CD163_sampled_CD163_i528_S119_L006_R1_001.fastq','ZJ31_63m_20191021_CD14pos_bulk_mono___CD181_sampled_I507_S7_L008_R1_001.fastq','ZJ31_75m_20201020_Mono__PB_CD193_sampled_i512_03_S3_L001_R1_001.fastq')
zj31.time <- NULL
zj31.time$grans=c(6,12,20,28,33,42,44,52.5,63,75)
zj31.time$t=c(6,12,20,28,33,42,44,52.5,63,75)
zj31.time$b=c(6,12,20,28,33,42,44,52.5,63,75)
zj31.time$mono=c(6,12,20,28,33,42,44,52.5,63,75)


#ZJ38
zj38.counts <- read.delim(file.path(my_path,"ZJ38_combined_20240620_counts.txt"), header=T, row.names=1, sep="\t")
#Grans
zj38.grans.samples <- c('zj389hmCD33PGFPNGFPPCD93R9.fastq','zj38_12m_20150702sort_Grans_ficolled_CD102_R5_ACAGTG_L006_R1_001.fastq','zj38_21m_Gr_Ficoll_20160412_pp_50trim_CD126_R1.fastq','zj38_24m_Gr_Ficoll_20160706_pp_50trim_CD136_R28.fastq','ZJ38_44m_2018_03_14_Grans___FX48_sampled_FX49_i521_S16_L003_R1_001.fastq','ZJ38_56hm_20191031_Gr_PB__CD183_sampled_CD183_I523_S43_L008_R1_001.fastq','ZJ38_99m_2022912_Gr_LT__CD226_sampled_CD226_95_S95_L002_R1_001.fastq','ZJ38_117m_202448_Gr_LT__CD226_sampled_CD226_93_S93_L002_R1_001.fastq')
zj38.t.samples <- c('zj386mTCD83R5.fastq','zj38_12m_20150702sort_T_CD102_R10_TAGCTT_L006_R1_001.fastq','zj38_21m_T_20160412_pp_50trim_CD126_R2.fastq','zj38_24m_T_20160706_pp_50trim_CD136_R21.fastq','ZJ38_44m_2018_03_14_T___FX48_sampled_FX49_i518_S13_L003_R1_001.fastq','ZJ38_56hm_20191031_T_PB__CD183_sampled_CD183_I524_S44_L008_R1_001.fastq','ZJ38_107m_202367_T_PB_NA_CD223_sampled_CD223_29_S29_L001_R1_001.fastq','ZJ38_117m_202449_T_LT__CD226_sampled_CD226_84_S84_L002_R1_001.fastq')
zj38.b.samples <- c('zj386mBCD83R6.fastq','zj38_12m_20150702sort_B_CD102_R8_ACTTGA_L006_R1_001.fastq','zj38_21m_B_20160412_pp_50trim_CD126_R3.fastq','zj38_24m_B_20160706_pp_50trim_CD136_R23.fastq','ZJ38_56hm_20191031_B_PB__CD183_sampled_CD183_I526_S45_L008_R1_001.fastq','ZJ38_107m_202367_B_PB_NA_CD222_sampled_AG_8963_30_S30_L001_R1_001.fastq','ZJ38_117m_202449_B_LT__CD226_sampled_CD226_85_S85_L002_R1_001.fastq')
zj38.mono.samples <- c('zj389hmMonoCD93R13.fastq','zj38_12m_20150702sort_Mono_CD102_R7_CAGATC_L006_R1_001.fastq','zj38_21m_Mono_20160412_pp_50trim_CD126_R4.fastq','zj38_24m_Mono_20160706_pp_50trim_CD136_R12.fastq','ZJ38_44m_2018_03_14_Mono___FX48_sampled_FX49_i519_S14_L003_R1_001.fastq','ZJ38_56hm_20191031_Mono_PB__CD183_sampled_CD183_I527_S46_L008_R1_001.fastq','ZJ38_107m_202367_Mono_PB_NA_CD222_sampled_AG_8963_31_S31_L001_R1_001.fastq','ZJ38_117m_202449_Mono_LT__CD226_sampled_CD226_86_S86_L002_R1_001.fastq')
zj38.time <- NULL
zj38.time$grans=c(9.5,12,21,24,44,64,99,117)
zj38.time$t=c(6,12,21,24,44,64,107,117)
zj38.time$b=c(6,12,21,24,64,107,117)
zj38.time$mono=c(9.5,12,21,24,44,64,107,117)


auto_plot_grans<-(autocorrelation.plot(list(zg66.counts[,zg66.grans.samples],
                                      zh33.counts[,zh33.grans.samples],
                                      zh19.counts[,zh19.grans.samples],
                                      zj31.counts[,zj31.grans.samples],
                                      zj38.counts[,zj38.grans.samples]),
                                 list(zg66.time$grans,
                                      zh33.time$grans,
                                      zh19.time$grans,
                                      zj31.time$grans,
                                      zj38.time$grans),
                                 Monkey=c("ZG66","ZH33","ZH19","ZJ31","ZJ38"),"Grans"))

auto_plot_t<-(autocorrelation.plot(list(zg66.counts[,zg66.t.samples],
                                      zh33.counts[,zh33.t.samples],
                                      zh19.counts[,zh19.t.samples],
                                      zj31.counts[,zj31.t.samples],
                                      zj38.counts[,zj38.t.samples]),
                                 list(zg66.time$t,
                                      zh33.time$t,
                                      zh19.time$t,
                                      zj31.time$t,
                                      zj38.time$t),
                                 Monkey=c("ZG66","ZH33","ZH19","ZJ31","ZJ38"),"T Cell"))

auto_plot_b<-(autocorrelation.plot(list(zg66.counts[,zg66.b.samples],
                                      zh33.counts[,zh33.b.samples],
                                      zh19.counts[,zh19.b.samples],
                                      zj31.counts[,zj31.b.samples],
                                      zj38.counts[,zj38.b.samples]),
                                 list(zg66.time$b,
                                      zh33.time$b,
                                      zh19.time$b,
                                      zj31.time$b,
                                      zj38.time$b),
                                 Monkey=c("ZG66","ZH33","ZH19","ZJ31","ZJ38"),"B Cell"))

auto_plot_mono<-(autocorrelation.plot(list(zg66.counts[,zg66.mono.samples],
                                        zh33.counts[,zh33.mono.samples],
                                        zh19.counts[,zh19.mono.samples],
                                        zj31.counts[,zj31.mono.samples],
                                        zj38.counts[,zj38.mono.samples]),
                                   list(zg66.time$mono,
                                        zh33.time$mono,
                                        zh19.time$mono,
                                        zj31.time$mono,
                                        zj38.time$mono),
                                   Monkey=c("ZG66","ZH33","ZH19","ZJ31","ZJ38"),"Mono"))