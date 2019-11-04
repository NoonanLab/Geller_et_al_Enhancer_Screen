library(rtracklayer)
library(ggplot2)

#import H3K27ac enhancer annotation of proliferation phenotypes
enhancer_set <- read.table( "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_enhancer_annotate_obs.txt", header= TRUE)
enhancer_set <- enhancer_set[enhancer_set$n.obs > 0,]
enhancer_set <- enhancer_set[enhancer_set$class == 'ac',]

enhancer_ac <- enhancer_set[enhancer_set$class == 'ac',]
enhancer_shuffle <- enhancer_set[enhancer_set$class == 'shuffle',]
enhancer_shuffle_grange <- makeGRangesFromDataFrame(enhancer_shuffle)

#import conserved regions targeted for genetic disruption
S0X_exp <- read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA.txt", header= TRUE)

#assign background set to be all targeted regions excluding proliferation phenotypes filtered for sign and timepoint inconsistency
background_set <- S0X_exp[(S0X_exp$enhancer == 1) & (S0X_exp$t12.phenotype != 'signIR') & (S0X_exp$t12.phenotype != 'dynamic'),]
background_set$phenotype.perm <- 0
background_set[background_set$t12.phenotype != 'neutral',]$phenotype.perm <- 1


permutation_enhancer <- data.frame(background_set$name)

#generate 100,000 permutations of proliferation phenotypes across the full dataset on enhancer genetic disruptions
for(i in 1:100000) {
	permutation_enhancer <- cbind(permutation_enhancer, sample(background_set$phenotype.perm))
}

enhancer_granges <- makeGRangesFromDataFrame(enhancer_set, keep.extra.columns= TRUE)
background_granges <- makeGRangesFromDataFrame(background_set, keep.extra.columns = TRUE)

#calculate permutation-based p-value based on the number of occurences within the permutations that equal or exceed the number of proliferation phenotypes observed within the enhancer
enhancer_granges$p.val <- 'NA'

for(i in c(1:length(enhancer_granges))) {
	sgRNAs <- background_granges[background_granges %over% enhancer_granges[i]]
	sgRNAs_permutation <- permutation_enhancer[permutation_enhancer[,1] %in% sgRNAs$name,c(-1)]
	n.perm <- sum(colSums(sgRNAs_permutation) >= enhancer_granges[i]$n.phenotypes.obs)
	if(n.perm > 0) {
		enhancer_granges[i]$p.val <- n.perm / 100000
	}
	else {
		enhancer_granges[i]$p.val <- 1 / 100001
	}
}
enhancer_granges$p.val.adjust <- p.adjust(enhancer_granges$p.val, method = 'BH')

output <- as.data.frame(enhancer_granges, row.names = NULL)


write.table(output, "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_enhancer_annotate_obs_perm_shuf_pval.txt", row.names=F, sep=" ", quote= F)