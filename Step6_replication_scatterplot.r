library(ggbio)
library(ggplot2)
library(rtracklayer)

#import rep1, rep2, and jointly-modeled biological effects
S0X_rep1 <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_rep1.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_rep2 <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_rep2.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_exp <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_filtered_spec.txt", header= TRUE), keep.extra.columns= TRUE)

S0X_rep1_filter <- S0X_rep1[S0X_rep1$name %in% S0X_exp$name]
S0X_rep2_filter <- S0X_rep2[S0X_rep2$name %in% S0X_exp$name]

S0X_exp <- S0X_exp[order(S0X_exp$name)]
S0X_rep1_filter <- S0X_rep1_filter[order(S0X_rep1_filter$name)]
S0X_rep2_filter <- S0X_rep2_filter[order(S0X_rep2_filter$name)]
S0X_full <- granges(S0X_exp)

#assign biological effects to each category
S0X_full$enhancer <- S0X_exp$enhancer
S0X_full$gene <- S0X_exp$gene
S0X_full$spec_control <- S0X_exp$spec_control
S0X_full$shuffle_control <- S0X_exp$shuffle_control

S0X_full$t4.beta.rep1 <- S0X_rep1_filter$t4.beta
S0X_full$t4.beta.rep2 <- S0X_rep2_filter$t4.beta
S0X_full$t8.beta.rep1 <- S0X_rep1_filter$t8.beta
S0X_full$t8.beta.rep2 <- S0X_rep2_filter$t8.beta
S0X_full$t12.beta.rep1 <- S0X_rep1_filter$t12.beta
S0X_full$t12.beta.rep2 <- S0X_rep2_filter$t12.beta
S0X_full$name <- S0X_exp$name


S0X_full$t4.phenotype <- as.character(S0X_exp$t4.phenotype)
S0X_full$t8.phenotype <- as.character(S0X_exp$t8.phenotype)
S0X_full$t12.phenotype <- as.character(S0X_exp$t12.phenotype)



S0X_pCDS <- S0X_full[((S0X_full$gene == 1 ) | (S0X_full$shuffle_control == 1)) & (S0X_full$t12.phenotype != 'dynamic') & (S0X_full$spec_control == 0)]

#assign phenotype categories
S0X_pCDS$plot_group <- 'Targeted Region'
S0X_pCDS[S0X_pCDS$t12.phenotype == 'negative']$plot_group <- 'Negative'
S0X_pCDS[S0X_pCDS$t12.phenotype == 'positive']$plot_group <- 'Positive'
S0X_pCDS[S0X_pCDS$shuffle_control == 1]$plot_group <- 'Genomic Background'
S0X_pCDS$plot_group <- factor(S0X_pCDS$plot_group, c('Targeted Region', 'Negative','Positive','Genomic Background'))
    
#plot rep1, rep2 gene genetic disruptions
scatter_t12 <- ggplot(S0X_pCDS[order(S0X_pCDS$plot_group)], aes(x = t12.beta.rep1, y = t12.beta.rep2, color = plot_group)) +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.5, color="#d9d9d9") + geom_vline(xintercept= 0, size = 0.5, color="#d9d9d9") + geom_abline(slope=1, intercept=0, size = 0.5, color="#d9d9d9") +
    geom_point(shape = 20, size = 1.6) + theme(legend.position="none") +
    scale_color_manual(values=c('#231F20','#295ba9','#008d4b','#D1D3D4')) +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group), data =  as.data.frame(subset(S0X_pCDS, grepl('^UBE2A$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^TCF7L1$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^CCND2$', name)))[1,], size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^SOX2$', name)))[1,], size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^DISC1$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^FGFR1$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^CEP135$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^MCPH1$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^CHD8$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^DIP2A$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^ASH2L$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^DYRK1A$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
    geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_pCDS, grepl('^ASPM$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	scale_fill_manual(values=c('#1B75BC','#009444')) + xlim(-3, 1.5 ) + ylim(-3, 1.5) + theme(aspect.ratio=0.9)

ggsave("~/Desktop/enhancer_analysis/figures/Scatter_t12_pCDS.pdf", plot = scatter_t12)


S0X_enhancer <- S0X_full[((S0X_full$enhancer == 1 ) | (S0X_full$shuffle_control == 1)) & (S0X_full$t12.phenotype != 'dynamic') & (S0X_full$spec_control == 0)]

#assign phenotype categories
S0X_enhancer$plot_group <- 'Targeted Region'
S0X_enhancer[S0X_enhancer$t12.phenotype == 'negative']$plot_group <- 'Negative'
S0X_enhancer[S0X_enhancer$t12.phenotype == 'positive']$plot_group <- 'Positive'
S0X_enhancer[S0X_enhancer$shuffle_control == 1]$plot_group <- 'Genomic Background'
S0X_enhancer$plot_group <- factor(S0X_enhancer$plot_group, c('Targeted Region', 'Negative', 'Positive','Genomic Background'))
    
#plot rep1, rep2 enhancer genetic disruptions
scatter_t12 <- ggplot(S0X_enhancer[order(S0X_enhancer$plot_group)], aes(x = t12.beta.rep1, y = t12.beta.rep2, color = plot_group)) +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.5, color="#d9d9d9") + geom_vline(xintercept= 0, size = 0.5, color="#d9d9d9") + geom_abline(slope=1, intercept=0, size = 0.5, color="#d9d9d9") +
    geom_point(shape = 20, size = 1.6) +
	scale_color_manual(values=c('#231F20','#295ba9','#008d4b','#D1D3D4')) + theme(legend.position="none")  +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S02_phastCons_429$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S03_phastCons_1339$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S03_phastCons_4251$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S05_phastCons_3882$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S05_phastCons_4245$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S05_phastCons_4089$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S03_phastCons_2954$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S02_phastCons_1572$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S02_phastCons_2676$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S01_phastCons_2069$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S06_phastCons_3160$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S01_phastCons_3008$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S03_phastCons_3417$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	geom_point(mapping = aes(x = t12.beta.rep1, y = t12.beta.rep2, fill = plot_group),data =  as.data.frame(subset(S0X_enhancer, grepl('^S04_phastCons_312$', name))), size = 0.8, stroke = 1.2, color = '#fd8d3c') +
	scale_fill_manual(values=c('#1B75BC','#009444'))  + xlim(-2, 1.5 ) + ylim(-2, 1.5) + theme(aspect.ratio=0.9)



ggsave("~/Desktop/enhancer_analysis/figures/Scatter_t12_enhancer.pdf")
