library(ggbio)
library(rtracklayer)
library(RColorBrewer)
library(tidyverse)


S0X_exp <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_filtered_spec.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_exp <- S0X_exp[!duplicated(S0X_exp$name)]


S0X_exp$t4.phenotype <- as.character(S0X_exp$t4.phenotype)
S0X_exp$t8.phenotype <- as.character(S0X_exp$t8.phenotype)
S0X_exp$t12.phenotype <- as.character(S0X_exp$t12.phenotype)

S0X_noSpec <- S0X_exp[S0X_exp$spec_control == 0]
S0X_phenotype <- S0X_noSpec

S0X.pca <- prcomp(mcols(S0X_phenotype[,c("t4.beta", "t8.beta", "t12.beta")]), center = FALSE, scale =FALSE) 

summary(S0X.pca)

S0X_phenotype$PC1 <- S0X.pca$x[,1]
S0X_phenotype$PC2 <- -1*S0X.pca$x[,2]
S0X_phenotype$PC3 <- S0X.pca$x[,3]


S0X_output <- as.data.frame(S0X_phenotype)
write.table(S0X_output, file = "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA.txt", sep = " ", quote = FALSE, row.names = FALSE)


S0X_pCDS <- S0X_phenotype[(S0X_phenotype$gene == 1)]
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t4.phenotype != 'signIR']
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t8.phenotype != 'signIR']
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t12.phenotype != 'signIR']



p <- ggplot(S0X_pCDS,aes(x=PC1,y=PC2, color = factor(t12.phenotype)))
p <- p + theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.4, color="#d9d9d9") + geom_vline(xintercept= 0, size = 0.4, color="#d9d9d9")  +lims(x= c(-3.5, 0.75), y = c(-0.6,0.6)) +
	geom_point(size = 0.4,  color = "#969696") + 
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_pCDS, grepl('neutral', t12.phenotype))), color = "#d9d9d9") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_pCDS, grepl('positive', t12.phenotype))), color = "#008d4b") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_pCDS, grepl('negative', t12.phenotype))), color = "#295ba9") +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_pCDS, grepl('^UBE2A$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_pCDS, grepl('^KIF20B$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_pCDS, grepl('^WDR12$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_pCDS, grepl('^FGFR1$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c')


ggsave("~/Desktop/enhancer_analysis/plots/PCA_pCDS_2d.pdf")
	
S0X_enhancer <- S0X_phenotype[S0X_phenotype$enhancer == 1]
S0X_enhancer <- S0X_enhancer[S0X_enhancer$t4.phenotype != 'signIR']
S0X_enhancer <- S0X_enhancer[S0X_enhancer$t8.phenotype != 'signIR']
S0X_enhancer <- S0X_enhancer[S0X_enhancer$t12.phenotype != 'signIR']

p <- ggplot(S0X_enhancer,aes(x=PC1,y=PC2, color = factor(t12.phenotype)))
p <- p + theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.4, color="#d9d9d9") + geom_vline(xintercept= 0, size = 0.4, color="#d9d9d9")  +
	geom_point(size = 0.4,  color = "#969696") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_enhancer, grepl('neutral', t12.phenotype))), color = "#d9d9d9") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_enhancer, grepl('positive', t12.phenotype))), color = "#008d4b") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_enhancer, grepl('negative', t12.phenotype))), color = "#295ba9") +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_enhancer, grepl('^S04_phastCons_311$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_enhancer, grepl('^S06_phastCons_3160$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_enhancer, grepl('^S02_phastCons_4483$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_enhancer, grepl('^S03_phastCons_4251$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_enhancer, grepl('^S03_phastCons_3417$', name))),  size = 0.6, stroke = 1.2, color = '#fd8d3c') 

ggsave("~/Desktop/enhancer_analysis/plots/PCA_enhancer_2d.pdf")


S0X_pCDS_neg <- S0X_pCDS[S0X_pCDS$t12.phenotype == 'negative']

pCDS_neg_quantiles<- quantile(S0X_pCDS_neg$PC1, probs = c(0.25, 0.5,1))

S0X_pCDS_neg <- as.data.frame(S0X_pCDS_neg)
hist_plot <- ggplot(S0X_pCDS_neg, aes(x=PC1, y = ..count.. + 0.1))
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_pCDS_neg, PC1 <= pCDS_neg_quantiles[3]), binwidth = 0.04, fill = '#023858')
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_pCDS_neg, PC1 <= pCDS_neg_quantiles[2]), binwidth = 0.04, fill = '#045a8d')
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_pCDS_neg, PC1 <= pCDS_neg_quantiles[1]), binwidth = 0.04, fill = '#0570b0')
hist_plot <- hist_plot +  scale_y_log10(expand = c(0, 0), limits = c(1, 1000), breaks = c(10,100,1000))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
				axis.title.y=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
hist_plot <- hist_plot + coord_fixed(ratio = 0.6) + annotation_logticks(sides = 'l') 

ggsave("~/Desktop/enhancer_analysis/plots/pCDS_neg_hist.pdf", width=3.25, height = 1.5)


S0X_pCDS_plot <- as.data.frame(S0X_pCDS)
hist_plot <- ggplot(S0X_pCDS_plot, aes(x=PC1, y = ..count.. + 0.1))
hist_plot <- hist_plot + geom_histogram(binwidth = 0.04, fill = '#525252')
hist_plot <- hist_plot + geom_histogram(data = subset(as.data.frame(S0X_phenotype), shuffle_control == 1), binwidth = 0.04, fill = '#bdbdbd')
hist_plot <- hist_plot +  scale_y_log10(expand = c(0, 0), limits = c(1, 1000), breaks = c(10,100,1000))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
				axis.title.y=element_blank(), axis.title.x = element_blank())
hist_plot <- hist_plot + coord_fixed(ratio = 0.6) + annotation_logticks(sides = 'l')

ggsave("~/Desktop/enhancer_analysis/plots/pCDS_background_hist.pdf", width=3.25)


S0X_enhancer <- S0X_phenotype[S0X_phenotype$enhancer == 1]
S0X_plot_enhancer <- S0X_enhancer[S0X_enhancer$t12.phenotype != 'IR']
S0X_enhancer_neg <- as.data.frame(S0X_enhancer[S0X_enhancer$t12.phenotype == 'negative'])
enhancer_neg_quantiles<- quantile(S0X_enhancer_neg$PC1, probs = c(0.25, 0.5,1))



hist_plot <- ggplot(S0X_enhancer_neg, aes(x=PC1, y = ..count.. + 0.1))
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_enhancer_neg, PC1 <= enhancer_neg_quantiles[3]), binwidth = 0.04, fill = '#023858')
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_enhancer_neg, PC1 <= enhancer_neg_quantiles[2]), binwidth = 0.04, fill = '#045a8d')
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_enhancer_neg, PC1 <= enhancer_neg_quantiles[1]), binwidth = 0.04, fill = '#0570b0')
hist_plot <- hist_plot +  scale_y_log10(expand = c(0, 0), limits = c(1, 1000), breaks = c(10,100,1000))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
				axis.title.y=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
hist_plot <- hist_plot + coord_fixed(ratio = 0.3) + annotation_logticks(sides = 'l') 

ggsave("~/Desktop/enhancer_analysis/plots/Enhancer_neg_hist.pdf", width=3.25, height = 1.5)



S0X_enhancer_plot <- as.data.frame(S0X_enhancer)

hist_plot <- ggplot(S0X_enhancer_plot, aes(x=PC1, y = ..count.. + 0.1))
hist_plot <- hist_plot + geom_histogram(data = subset(as.data.frame(S0X_pCDS_plot), gene == 1), binwidth = 0.04, fill = '#969696')
hist_plot <- hist_plot + geom_histogram(binwidth = 0.04, fill = '#525252')
hist_plot <- hist_plot + geom_histogram(data = subset(as.data.frame(S0X_phenotype), shuffle_control == 1), binwidth = 0.04, fill = '#bdbdbd')
hist_plot <- hist_plot +  scale_y_log10(expand = c(0, 0), limits = c(1, 1700), breaks = c(10,100,1000))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
				axis.title.y=element_blank(), axis.title.x = element_blank())
hist_plot <- hist_plot + coord_fixed(ratio = 0.5) + annotation_logticks(sides = 'l')

ggsave("~/Desktop/enhancer_analysis/plots/Enhancer_background_hist.pdf", width=3.25)


p <- ggplot(S0X_plot_enhancer,aes(x=PC1,y=PC2, color = factor(t12.phenotype)))
p <- p + theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.4, color="#d9d9d9") + geom_vline(xintercept= 0, size = 0.4, color="#d9d9d9")  +
	geom_point(size = 0.8,  color = "#bababa") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_plot_enhancer, grepl('neutral', t12.phenotype))), color = "#878787") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_plot_enhancer, grepl('positive', t12.phenotype))), color = "#008d4b") +
	geom_density_2d(size = 0.8, data = as.data.frame(subset(S0X_plot_enhancer, grepl('negative', t12.phenotype))), color = "#295ba9") +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_plot_enhancer, grepl('^S05_phastCons_318$', name))),  size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_plot_enhancer, grepl('^S03_phastCons_1483$', name))),  size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_plot_enhancer, grepl('^S03_phastCons_3644$', name))),  size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = PC1, y = PC2),data =  as.data.frame(subset(S0X_plot_enhancer, grepl('^S04_phastCons_1118$', name))),  size = 0.8, stroke = 1.2, color = '#fd8d3c') 



ggsave("~/Desktop/enhancer_analysis/plots/PCA_enhancer_2d.pdf")



