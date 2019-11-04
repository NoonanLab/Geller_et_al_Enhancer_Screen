library(rtracklayer)
library(ggplot2)

#import enhancer-level proliferation phenotype summary
enhancer_set <- read.table( "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_enhancer_annotate_obs.txt", header= TRUE)

ac_set <- enhancer_set[enhancer_set$class == 'ac',]

ac_set$phenotype_density <- ac_set$n.phenotypes.obs / ac_set$n.obs

ac_set$feature.plot <- 'Enhancer'
ac_set[ac_set$GAIN == 1,]$feature.plot <- 'Gain'
ac_set[ac_set$HACNS_HAR == 1,]$feature.plot <- 'HAR'

ac_set$feature.plot <- factor(ac_set$feature.plot, c('Enhancer', 'Gain', 'HAR'))

#select H3K27ac enhancers harboring at least one proliferation phenotype
phenotype_set <- ac_set[ac_set$n.phenotypes.obs > 0,]

phenotype_set_har_gain <- phenotype_set[(phenotype_set$GAIN == 1) | (phenotype_set$HACNS_HAR == 1),]

#generate plots for enhancer-level summary characteristics

#Quantity of proliferation phenotypes residing within H3K27ac enhancers
histogram_counts <- ggplot(phenotype_set, aes(x = n.phenotypes.obs )) +  geom_histogram(binwidth = 1.0, alpha = 0.75, fill = "#bdbdbd") +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_histogram(data = phenotype_set_har_gain, binwidth = 1.0, fill = "#ef3b2c")

#Quantity of proliferation phenotypes residing within HGE/HAR enhancers
histogram_counts_hge_har <- ggplot(phenotype_set_har_gain, aes(x = n.phenotypes.obs )) +  geom_histogram(binwidth = 1.0, alpha = 0.75, fill = "#bdbdbd") +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_histogram(data = phenotype_set_har_gain, binwidth = 1.0, fill = "#ef3b2c")

#Density of proliferation phenotypes residing within H3K27ac enhancers
histogram_density <- ggplot(phenotype_set, aes(x = phenotype_density )) +  geom_histogram(binwidth = 0.05,  alpha = 0.75, fill = "#bdbdbd") +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_vline(aes(xintercept=mean(phenotype_density, na.rm=T)),color="#252525", linetype="dashed", size=0.45) +
	geom_histogram(data = phenotype_set_har_gain, binwidth = 0.05,  fill = "#ef3b2c")

#Density of proliferation phenotypes residing within HGE/HAR enhancers
histogram_density_hge_har <- ggplot(phenotype_set_har_gain	, aes(x = phenotype_density )) +  geom_histogram(binwidth = 0.05,  alpha = 0.75, fill = "#bdbdbd") +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_vline(aes(xintercept=mean(phenotype_density, na.rm=T)),color="#252525", linetype="dashed", size=0.45) +
	geom_histogram(data = phenotype_set_har_gain, binwidth = 0.05,  fill = "#ef3b2c")

#Cumulative proliferation burden versus number of conserved regions interrogated within each enhancer
proliferation_burden <- ggplot(phenotype_set, aes (x = n.obs , y = PC1.burden, color = feature.plot)) + geom_point(size= 0.8,alpha = 0.75, color = c("#525252"))  +
	geom_point(data = phenotype_set_har_gain, size = 0.8,color = c("#ef3b2c")) +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=19), axis.title.x=element_blank(), axis.title.y=element_blank())
