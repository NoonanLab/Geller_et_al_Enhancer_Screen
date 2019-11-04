library(ggbio)
library(ggplot2)
library(rtracklayer)

#import jointly-modeled biological effects for all genetic disruptions
S0X_exp <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_filtered_spec.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_exp <- S0X_exp[!duplicated(S0X_exp$name)]
S0X_exp <- S0X_exp[S0X_exp$spec_control != 1]

#remove sign irreproducible regions
S0X_exp <- S0X_exp[S0X_exp$t4.phenotype != 'signIR']
S0X_exp <- S0X_exp[S0X_exp$t8.phenotype != 'signIR']
S0X_exp <- S0X_exp[S0X_exp$t12.phenotype != 'signIR']

#assign categories to all biological effects
S0X_exp$feature.plot <- 'NA'
S0X_exp[(S0X_exp$enhancer == 1)]$feature.plot <- 'Enhancer'
S0X_exp[(S0X_exp$gene == 1)]$feature.plot <- 'Gene'
S0X_exp[(S0X_exp$shuffle_control == 1)]$feature.plot <- 'Background'

S0X_exp <- S0X_exp[S0X_exp$feature.plot != 'NA']
S0X_exp$feature.plot <- factor(S0X_exp$feature.plot, c('Background', 'Enhancer', 'Gene'))

#generate jitter plot at 4 cell divisions
jitter_t4 <- ggplot(S0X_exp, aes(x=feature.plot, y = t4.beta, color = feature.plot))  +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.4, color="#bdbdbd") +
	geom_hline(yintercept= -1, linetype = "dashed",size = 0.4, color="#bdbdbd") +
	geom_hline(yintercept= 1, linetype = "dashed",size = 0.4, color="#bdbdbd") +
	geom_jitter(width = 0.4, height = 0.25, size = 1.0) +
	geom_boxplot(color = "black", size = 0.5, outlier.shape = NA , width = 0.4) +
	stat_boxplot(geom = "errorbar",color = "black", size = 0.5, width = 0.2) +
    scale_color_manual(values=c('#bdbdbd','#DE1D40', '#737373')) + theme(legend.position="none")  + ylim(-2,2) +
    coord_fixed(ratio = 1.2)
    
ggsave("~/Desktop/enhancer_analysis/figures/Category_jitter_t4.pdf")

#generate jitter plot at 8 cell divisions
jitter_t8 <- ggplot(S0X_exp[S0X_exp$feature.plot != 'NA'], aes(x=feature.plot, y = t8.beta, color = feature.plot))  +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.4, color="#bdbdbd") +
	geom_hline(yintercept= -1, linetype = "dashed",size = 0.4, color="#bdbdbd") +
	geom_hline(yintercept= 1, linetype = "dashed",size = 0.4, color="#bdbdbd") +
	geom_jitter(width = 0.4, height = 0.25, size = 1.0) +
	geom_boxplot(color = "black", size = 0.5, outlier.shape = NA , width = 0.4) +
	stat_boxplot(geom = "errorbar",color = "black", size = 0.5, width = 0.2) +
    scale_color_manual(values=c('#bdbdbd','#DE1D40', '#737373')) + theme(legend.position="none")  + ylim(-2,2) +
    coord_fixed(ratio = 1.2)
    
ggsave("~/Desktop/enhancer_analysis/figures/Category_jitter_t8.pdf")

#generate jitter plot at 12 cell divisions
jitter_t12 <- ggplot(S0X_exp[S0X_exp$feature.plot != 'NA'], aes(x=feature.plot, y = t12.beta, color = feature.plot))  +
	theme(panel.background = element_rect(fill = "#ffffff"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1), text = element_text(size=20), axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank()) +
	geom_hline(yintercept= 0, size = 0.4, color="#bdbdbd") +
	geom_hline(yintercept= -1, linetype = "dashed",size = 0.4, color="#bdbdbd") +
	geom_hline(yintercept= 1, linetype = "dashed",size = 0.4, color="#bdbdbd") +
	geom_jitter(width = 0.4, height = 0.25, size = 1.0) +
	geom_boxplot(color = "black", size = 0.5, outlier.shape = NA , width = 0.4) +
	stat_boxplot(geom = "errorbar",color = "black", size = 0.5, width = 0.2) +
    scale_color_manual(values=c('#bdbdbd','#DE1D40', '#737373')) + theme(legend.position="none")  + ylim(-2,2) +
    coord_fixed(ratio = 1.2)


ggsave("~/Desktop/enhancer_analysis/figures/Category_jitter_t12.pdf")
