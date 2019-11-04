library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)

#import reference bed files
hNSC_ac <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/hNSC_H3K27ac.bed")
POS_REF <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/wang_et_al_Pos_control_pCDS_gRNA_e29m1_sort_shuf_sort.bed")
POS_REF <- subset(POS_REF, select=-c(score))

SHUFFLE_REF <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/hNSC_Emb_brain_Shuffle_1000shuf_mapped500.bed")
SHUFFLE_FILTER <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/hNSC_Emb_brain_Shuffle_1000shuf_mapped500_filtered.bed")
SHUFFLE_EXCLUDE_ENHANCERS <- SHUFFLE_REF[!SHUFFLE_REF %over% SHUFFLE_FILTER]
SHUFFLE_EXCLUDE_TILES <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/Shuffle_gRNA_e29m1_sort_select_overDNAse.bed")

HACNS_REF <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/merge_HACNS_HARS.bed")
gain1 <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_7pcw_ac_Hu_gain.bed")
gain2 <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_8_5pcw_ac_Hu_gain.bed")
gain3 <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_12Fpcw_ac_Hu_gain.bed")
gain4 <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_12Opcw_ac_Hu_gain.bed")

vista <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/vista_enhancers.bed")

#import mageck biological effects from joint modeling across replicates
S01 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S01_NTC.gene_summary.txt", header= T)
S02 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S02_NTC.gene_summary.txt", header= T)
S03 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S03_NTC.gene_summary.txt", header= T)
S04 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S04_NTC.gene_summary.txt", header= T)
S05 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S05_NTC.gene_summary.txt", header= T)
S06 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S06_NTC.gene_summary.txt", header= T)
S07 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S07_NTC.gene_summary.txt", header= T)
S08 <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/beta/S08_NTC.gene_summary.txt", header= T)

colnames(S01)[1] <- "name"
colnames(S02)[1] <- "name"
colnames(S03)[1] <- "name"
colnames(S04)[1] <- "name"
colnames(S05)[1] <- "name"
colnames(S06)[1] <- "name"
colnames(S07)[1] <- "name"
colnames(S08)[1] <- "name"\

#import map of disrupted regions genomic coordinates
S01_map <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/mageck/lib/S01_lib.map")
S02_map <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/mageck/lib/S02_lib.map")
S03_map <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/mageck/lib/S03_lib.map")
S04_map <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/mageck/lib/S04_lib.map")
S05_map <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/mageck/lib/S05_lib.map")
S06_map <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/mageck/lib/S06_lib.map")
S07_map <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/pCDS_gRNA_e29m1_sort_hNSC_2gRNA_batch00.bed")
S08_map <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/pCDS_gRNA_e29m1_sort_hNSC_2gRNA_batch01.bed")

S07_map <- subset(S07_map, select=-c(score))
S07_map <- S07_map[!duplicated(S07_map$name),]

S08_map <- subset(S08_map, select=-c(score))
S08_map <- S08_map[!duplicated(S08_map$name),]

S01_map <- c(S01_map, POS_REF)
S02_map <- c(S02_map, POS_REF)
S03_map <- c(S03_map, POS_REF)
S04_map <- c(S04_map, POS_REF)
S05_map <- c(S05_map, POS_REF)
S06_map <- c(S06_map, POS_REF)
S07_map <- c(S07_map, POS_REF)
S08_map <- c(S08_map, POS_REF)


S01_mageck <- makeGRangesFromDataFrame(merge(S01_map, S01,by="name"), keep.extra.columns= TRUE)
S02_mageck <- makeGRangesFromDataFrame(merge(S02_map, S02,by="name"), keep.extra.columns= TRUE)
S03_mageck <- makeGRangesFromDataFrame(merge(S03_map, S03,by="name"), keep.extra.columns= TRUE)
S04_mageck <- makeGRangesFromDataFrame(merge(S04_map, S04,by="name"), keep.extra.columns= TRUE)
S05_mageck <- makeGRangesFromDataFrame(merge(S05_map, S05,by="name"), keep.extra.columns= TRUE)
S06_mageck <- makeGRangesFromDataFrame(merge(S06_map, S06,by="name"), keep.extra.columns= TRUE)
S07_mageck <- makeGRangesFromDataFrame(merge(S07_map, S07,by="name"), keep.extra.columns= TRUE)
S08_mageck <- makeGRangesFromDataFrame(merge(S08_map, S08,by="name"), keep.extra.columns= TRUE)

S0X_mageck<- sort(c(S01_mageck, S02_mageck, S03_mageck, S04_mageck, S05_mageck, S06_mageck, S07_mageck, S08_mageck))

S0X_spec <- S04[grep("_r", S04$name),]
S0X_spec$chr <- S0X_spec$name
S0X_spec$start <- 1
S0X_spec$end <- 20

S0X_mageck <- c(S0X_mageck, makeGRangesFromDataFrame(S0X_spec, keep.extra.columns=TRUE))

#annotate integrated sub-libraries based on reference features
S0X_mageck$neg_control <- 0
S0X_mageck[S0X_mageck %over% POS_REF]$neg_control <- 1

S0X_mageck$shuffle_control <- 0
S0X_mageck[S0X_mageck %over% SHUFFLE_FILTER]$shuffle_control <- 1

S0X_mageck$spec_control <- 0
S0X_mageck[grep("_r", S0X_mageck)]$spec_control <- 1


S0X_mageck[(S0X_mageck$t4.p.value == 0)]$t4.p.value <- 1 / (1+10*length(unique(S04_mageck$name)))
S0X_mageck[(S0X_mageck$t8.p.value == 0)]$t8.p.value <- 1 / (1+10*length(unique(S04_mageck$name)))
S0X_mageck[(S0X_mageck$t12.p.value == 0)]$t12.p.value <- 1 / (1+10*length(unique(S04_mageck$name)))

S0X_mageck$ac <- 0
S0X_mageck[S0X_mageck %over% hNSC_ac,]$ac <- 1

S0X_mageck$HACNS_HAR <- 0
S0X_mageck[S0X_mageck %over% HACNS_REF]$HACNS_HAR <- 1

S0X_mageck$GAIN_7pcw <- 0
S0X_mageck[S0X_mageck %over% gain1]$GAIN_7pcw <- 1

S0X_mageck$GAIN_8.5pcw <- 0
S0X_mageck[S0X_mageck %over% gain2]$GAIN_8.5pcw <- 1

S0X_mageck$GAIN_12_F_pcw <- 0
S0X_mageck[S0X_mageck %over% gain3]$GAIN_12_F_pcw <- 1

S0X_mageck$GAIN_12_O_pcw <- 0
S0X_mageck[S0X_mageck %over% gain4]$GAIN_12_O_pcw <- 1

S0X_mageck$GAIN <- 0
S0X_mageck[(S0X_mageck$GAIN_7pcw == 1) | (S0X_mageck$GAIN_8.5pcw == 1)| (S0X_mageck$GAIN_12_F_pcw == 1) | (S0X_mageck$GAIN_12_O_pcw == 1)]$GAIN <- 1

S0X_mageck$vista <- 0
S0X_mageck[S0X_mageck %over% vista]$vista <- 1

S0X_mageck$enhancer <- 0
S0X_mageck[grep("phastCons", S0X_mageck$name)]$enhancer <- 1

S0X_mageck$gene <- 0
S0X_mageck[(S0X_mageck %over% S07_map) | (S0X_mageck %over% S08_map)]$gene <- 1

S0X_mageck <- S0X_mageck[!S0X_mageck %over% SHUFFLE_EXCLUDE_ENHANCERS]
S0X_mageck <- S0X_mageck[!S0X_mageck %over% SHUFFLE_EXCLUDE_TILES]

S0X_mageck$t4.z.score <- (S0X_mageck$t4.beta - mean(S0X_mageck[S0X_mageck$shuffle_control == 1]$t4.beta)) / sd(S0X_mageck[S0X_mageck$shuffle_control == 1]$t4.beta)
S0X_mageck$t8.z.score <- (S0X_mageck$t8.beta - mean(S0X_mageck[S0X_mageck$shuffle_control == 1]$t8.beta)) / sd(S0X_mageck[S0X_mageck$shuffle_control == 1]$t8.beta)
S0X_mageck$t12.z.score <- (S0X_mageck$t12.beta - mean(S0X_mageck[S0X_mageck$shuffle_control == 1]$t12.beta)) / sd(S0X_mageck[S0X_mageck$shuffle_control == 1]$t12.beta)

output <- as.data.frame(S0X_mageck, row.names = NULL, row.names = NULL)

write.table(output, "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_spec.txt", row.names=F, sep=" ", quote= F)

