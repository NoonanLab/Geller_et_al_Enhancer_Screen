library(rtracklayer)
library(ggplot2)

#import biological effects across non-coding regions with TF binding assignments
S0X_enhancer <-  makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA_TF_assign.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_enhancer <- S0X_enhancer[S0X_enhancer$enhancer == 1]
S0X_enhancer_neg <- S0X_enhancer[S0X_enhancer$t12.phenotype == 'negative']

jaspar_2018 <- import.bed("~/Desktop/enhancer_analysis/JASPAR_2018/JASPAR2018_hg19_shortName.bed")

#assign categories of biological effect based on percentile of proliferation-decreasing phenotype
enhancer_neg_quantiles<- quantile(S0X_enhancer_neg$PC1, probs = c(0.25, 0.5, 1.0))
names(enhancer_neg_quantiles) <- c('severe', 'strong', 'all')

#calculate overrepresentation enrichment by hypergeometric testing 
enhancer_hyperGeometric <- function(category, phenotype) {
	overlap_features <- S0X_enhancer[S0X_enhancer %over% category]
	nonoverlap_features <- S0X_enhancer[!S0X_enhancer %over% category]
	if(phenotype %in% c('severe', 'strong', 'all')) {
		success_in_sample <- length(overlap_features[(overlap_features$PC1 <= enhancer_neg_quantiles[phenotype]) & (overlap_features$t12.phenotype == 'negative')])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		failure_in_bkgd <- length(S0X_enhancer$t12.phenotype == 'neutral') - success_in_bkgd
		sample_size <- length(S0X_enhancer[(S0X_enhancer$PC1 <= enhancer_neg_quantiles[phenotype]) & (S0X_enhancer$t12.phenotype == 'negative')])
		result <- phyper(success_in_sample -1, success_in_bkgd, failure_in_bkgd, sample_size, lower.tail= FALSE);

	}
	
	else {
		success_in_sample <- length(overlap_features[overlap_features$t12.phenotype == 'positive'])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		failure_in_bkgd <- length(S0X_enhancer$t12.phenotype == 'neutral') - success_in_bkgd
		sample_size <- length(S0X_enhancer[S0X_enhancer$t12.phenotype == 'positive'])
		result <- phyper(success_in_sample -1, success_in_bkgd, failure_in_bkgd, sample_size, lower.tail= FALSE);

	}
	
	return(result)


}

#calculate the fold-enrichment of TFBS for all motifs
TF_foldEnrich <- function(category, phenotype){

	overlap_features <- S0X_enhancer[S0X_enhancer %over% category]
	nonoverlap_features <- S0X_enhancer[!S0X_enhancer %over% category]

	if(phenotype %in% c('severe', 'strong', 'all')) {
		success_in_sample <- length(overlap_features[(overlap_features$PC1 <= enhancer_neg_quantiles[phenotype]) & (overlap_features$t12.phenotype == 'negative')])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		sample_bkgd <- length(S0X_enhancer$t12.phenotype == 'neutral')
		sample_size <- length(S0X_enhancer[(S0X_enhancer$PC1 <= enhancer_neg_quantiles[phenotype]) & (S0X_enhancer$t12.phenotype == 'negative')])
		result <- (success_in_sample / sample_size) / (success_in_bkgd / sample_bkgd)
	

	}
	
	else {
		success_in_sample <- length(overlap_features[overlap_features$t12.phenotype == 'positive'])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		sample_bkgd <- length(S0X_enhancer$t12.phenotype == 'neutral')
		sample_size <- length(S0X_enhancer[S0X_enhancer$t12.phenotype == 'positive'])
		result <- (success_in_sample / sample_size) / (success_in_bkgd /sample_bkgd)

	}
	
	return(result)

}


tf.name <- as.data.frame(colnames(mcols(S0X_enhancer)[55:626]))
colnames(tf.name) <- 'name'

#generate overrepresentation results for severe proliferation-decreasing phenotypes
tf.name$severe.enrichment <- 'NA'
tf.name$severe.p.value <- 'NA'
for(tf in tf.name$name) {
	tf.name[tf.name$name == tf,]$severe.enrichment <- TF_foldEnrich(jaspar_2018[jaspar_2018$name == tf], 'severe')
	tf.name[tf.name$name == tf,]$severe.p.value <- enhancer_hyperGeometric(jaspar_2018[jaspar_2018$name == tf], 'severe')
}
tf.name$severe.p.value.adjust <- p.adjust(tf.name$severe.p.value, method = 'BH')

#generate overrepresentation results for strong proliferation-decreasing phenotypes
tf.name$strong.enrichment <- 'NA'
tf.name$strong.p.value <- 'NA'
for(tf in tf.name$name) {
	tf.name[tf.name$name == tf,]$strong.enrichment <- TF_foldEnrich(jaspar_2018[jaspar_2018$name == tf], 'strong')
	tf.name[tf.name$name == tf,]$strong.p.value <- enhancer_hyperGeometric(jaspar_2018[jaspar_2018$name == tf], 'strong')
}
tf.name$strong.p.value.adjust <- p.adjust(tf.name$strong.p.value, method = 'BH')

#generate overrepresentation results for all proliferation-decreasing phenotypes
tf.name$all.enrichment <- 'NA'
tf.name$all.p.value <- 'NA'
for(tf in tf.name$name) {
	tf.name[tf.name$name == tf,]$all.enrichment <- TF_foldEnrich(jaspar_2018[jaspar_2018$name == tf], 'all')
	tf.name[tf.name$name == tf,]$all.p.value <- enhancer_hyperGeometric(jaspar_2018[jaspar_2018$name == tf], 'all')
}
tf.name$all.p.value.adjust <- p.adjust(tf.name$all.p.value, method = 'BH')

#generate overrepresentation results for all proliferation-increasing phenotypes
tf.name$positive.enrichment <- 'NA'
tf.name$positive.p.value <- 'NA'
for(tf in tf.name$name) {
	tf.name[tf.name$name == tf,]$positive.enrichment <- as.numeric(TF_foldEnrich(jaspar_2018[jaspar_2018$name == tf], 'positive'))
	tf.name[tf.name$name == tf,]$positive.p.value <- as.numeric(enhancer_hyperGeometric(jaspar_2018[jaspar_2018$name == tf], 'positive'))
}
tf.name$positive.p.value.adjust <- p.adjust(tf.name$positive.p.value, method = 'BH')

output <- as.data.frame(tf.name, row.names = NULL)

write.table(output, "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA_TF_assign_hypergeometric_BH.txt", row.names=F, sep=" ", quote= F)