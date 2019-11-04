library(ggplot2)
library(rtracklayer)

#import reference files
enhancers <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/hNSC_Emb_brain_cluster_noPromoter.bed")
har_hacns <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/hNSC_Emb_brain_H3K27ac_noPromoter_merge_HACNS_HARS.bed")
shuffle <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/hNSC_Emb_brain_Shuffle_1000shuf_mapped500.sort.bed")
phastCons <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/phasCons_mammal.bed")
gain1 <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_7pcw_ac_Hu_gain.bed")
gain2 <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_8_5pcw_ac_Hu_gain.bed")
gain3 <- import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_12Fpcw_ac_Hu_gain.bed")
gain4 <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/GSE63648_12Opcw_ac_Hu_gain.bed")
vista <-import.bed("/Users/evangeller/Desktop/enhancer_analysis/bed/vista_enhancers.bed")

all_enhancers <- sort(c(enhancers, har_hacns))

#extend phastCon elements
start(phastCons) <- start(phastCons) - 15
end(phastCons) <- end(phastCons) + 15

#identify phastCon elements overlapping target set
phastCons_overlap <- phastCons[phastCons %over% all_enhancers]

all_enhancers_extend <- union(all_enhancers, phastCons_overlap, ignore.strand = TRUE)
all_enhancers_extend <- sort(c(all_enhancers_extend, shuffle))
start(all_enhancers_extend) <- start(all_enhancers_extend) - 15
end(all_enhancers_extend) <- end(all_enhancers_extend) + 15

#import jointly-modeled biological effects on proliferation
S0X_exp <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA.txt", header= TRUE), keep.extra.columns= TRUE)

#discard targeted regions that exhibit sign or timepoint inconsistently in proliferation
enhancer_tiles <- S0X_exp[(S0X_exp$enhancer == 1) | (S0X_exp$shuffle_control == 1),]
enhancer_tiles <- enhancer_tiles[enhancer_tiles$t4.phenotype != 'signIR']
enhancer_tiles <- enhancer_tiles[enhancer_tiles$t8.phenotype != 'signIR']
enhancer_tiles <- enhancer_tiles[enhancer_tiles$t12.phenotype != 'signIR']
enhancer_tiles <- enhancer_tiles[enhancer_tiles$t4.phenotype != 'dynamic']
enhancer_tiles <- enhancer_tiles[enhancer_tiles$t8.phenotype != 'dynamic']
enhancer_tiles <- enhancer_tiles[enhancer_tiles$t12.phenotype != 'dynamic']

enhancer_set <- all_enhancers_extend[all_enhancers_extend %over% enhancer_tiles]

#annotate enhancer intervals with features
enhancer_set$class <- 'NA'
enhancer_set[enhancer_set %over% all_enhancers]$class <- 'ac'
enhancer_set[enhancer_set %over%  shuffle]$class <- 'shuffle'
enhancer_set$HACNS_HAR <- 0
enhancer_set[enhancer_set %over% har_hacns]$HACNS_HAR <- 1
enhancer_set$GAIN_7pcw <- 0
enhancer_set[enhancer_set %over% gain1]$GAIN_7pcw <- 1
enhancer_set$GAIN_8.5pcw <- 0
enhancer_set[enhancer_set %over% gain2]$GAIN_8.5pcw <- 1
enhancer_set$GAIN_12_F_pcw <- 0
enhancer_set[enhancer_set %over% gain3]$GAIN_12_F_pcw <- 1
enhancer_set$GAIN_12_O_pcw <- 0
enhancer_set[enhancer_set %over% gain4]$GAIN_12_O_pcw <- 1
enhancer_set$GAIN <- 0
enhancer_set[(enhancer_set$GAIN_7pcw == 1) | (enhancer_set$GAIN_8.5pcw == 1)| (enhancer_set$GAIN_12_F_pcw == 1) | (enhancer_set$GAIN_12_O_pcw == 1)]$GAIN <- 1
enhancer_set$vista <- 0
enhancer_set[enhancer_set %over% vista]$vista <- 1
enhancer_categories <- quantile(enhancer_tiles[enhancer_tiles$t12.phenotype == 'negative']$PC1, probs = c(0.25, 0.5, 1))
names(enhancer_categories) <- c('severe', 'strong','all')


# for each enhancer summarize number of regions interrogated, number of proliferation phenotypes, total 'proliferation score' burden, and maximum magnitude proliferation phenotype

enhancer_count_tiles <- function(enhancer, timepoint) {
	sgRNAs <- enhancer_tiles[enhancer_tiles %over% enhancer]
	return(length(sgRNAs[(sgRNAs$t12.phenotype == 'neutral') | (sgRNAs$t12.phenotype == 'negative') | (sgRNAs$t12.phenotype == 'positive')]))
}

enhancer_count_phenotypes <- function(enhancer, timepoint) {
	sgRNAs <- enhancer_tiles[enhancer_tiles %over% enhancer]
	n_neg <- sum(mcols(sgRNAs)[,timepoint] == "negative")
	n_pos <- sum(mcols(sgRNAs)[,timepoint] == "positive")
	return(n_neg + n_pos)
}

enhancer_PC1_burden <- function(enhancer) {
	sgRNAs <- enhancer_tiles[enhancer_tiles %over% enhancer]
	return(sum(abs(mcols(sgRNAs)$PC1)))
}

enhancer_max_phenotype <- function(enhancer, timepoint) {
	sgRNAs <- enhancer_tiles[enhancer_tiles %over% enhancer]
	sgRNA_phenotypes <- sgRNAs[(sgRNAs$t12.phenotype == 'negative') | (sgRNAs$t12.phenotype == 'positive')]
	if(length(sgRNA_phenotypes) == 0) {
		return(0)
	}
	else {
		max_phenotype <- sgRNA_phenotypes[which.max(abs(sgRNA_phenotypes$PC1)),]
		return(max_phenotype$PC1)
	}
}

#compute summary for each H3K27ac enhancer included in the study
mcols(enhancer_set)$n.obs <- t(sapply(enhancer_set, enhancer_count_tiles, timepoint = 't12.phenotype'))
mcols(enhancer_set)$n.phenotypes.obs <- t(sapply(enhancer_set, enhancer_count_phenotypes, timepoint = 't12.phenotype'))
mcols(enhancer_set)$PC1.burden <- t(sapply(enhancer_set, enhancer_PC1_burden))
mcols(enhancer_set)$PC1.max.region <- t(sapply(enhancer_set, enhancer_max_phenotype, timepoint = 't12.phenotype'))


output <- as.data.frame(enhancer_set, row.names = NULL)

write.table(output, "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_enhancer_annotate_obs.txt", row.names=F, sep=" ", quote= F)