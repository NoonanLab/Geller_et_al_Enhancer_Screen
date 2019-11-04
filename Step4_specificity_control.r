library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(ggbio)

#import integrated library of biological effects following disruption and proliferation
S0X_batch <- makeGRangesFromDataFrame(read.table("~/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_filtered_spec.txt", header= TRUE), keep.extra.columns= TRUE)

#genes with predicted proliferation-decreasing biological effects
gene <- c('SOX2', 'SRSF1', 'CCND2', 'CTNNB1')

#extract all disruptions for each gene
SOX2 <- S0X_batch[grep('SOX2', S0X_batch$name)]
SOX2 <- SOX2[-c(58:59),]
SRSF1 <- S0X_batch[grep('SRSF1', S0X_batch$name)]
SRSF1 <- SRSF1[4:282,]
CCND2 <- S0X_batch[grep('CCND2', S0X_batch$name)]
CTNNB1 <- S0X_batch[grep('CTNNB1', S0X_batch$name)]

#compute full match sgRNA biological effects
SOX2.mean.full.match <- mean(SOX2[SOX2$gene == 1]$t12.beta)
SRSF1.mean.full.match <- mean(SRSF1[SRSF1$gene == 1]$t12.beta)
CCND2.mean.full.match <- mean(CCND2[CCND2$gene == 1]$t12.beta)
CTNNB1.mean.full.match <- mean(CTNNB1[CTNNB1$gene == 1]$t12.beta)

#compute proportion of full match sgRNA biolgoical effects
SOX2$t12.activity <- SOX2$t12.beta / SOX2.mean.full.match
SRSF1$t12.activity <- SRSF1$t12.beta / SRSF1.mean.full.match
CCND2$t12.activity <- CCND2$t12.beta / CCND2.mean.full.match
CTNNB1$t12.activity <- CTNNB1$t12.beta / CTNNB1.mean.full.match

#annotate sgRNA mismatches
SOX2$pattern <- 'Region'
SOX2[SOX2$gene == 1]$pattern <- 'Full Match'
SOX2[grep('r123_n1', SOX2$name)]$pattern <- 'R123, 1MM'
SOX2[grep('r1_n2', SOX2$name)]$pattern <- 'R1, 2MM'
SOX2[grep('r2_n2', SOX2$name)]$pattern <- 'R2, 2MM'
SOX2[grep('r3_n2', SOX2$name)]$pattern <- 'R3, 2MM'
SOX2[grep('r12_n2', SOX2$name)]$pattern <- 'R12, 2MM'
SOX2[grep('r1_n3', SOX2$name)]$pattern <- 'R1, 3MM'
SOX2[grep('r2_n3', SOX2$name)]$pattern <- 'R2, 3MM'
SOX2[grep('r3_n3', SOX2$name)]$pattern <- 'R3, 3MM'
SOX2[grep('r12_n3', SOX2$name)]$pattern <- 'R12, 3MM'
SOX2[grep('r123_n3', SOX2$name)]$pattern <- 'R123, 3MM'
SOX2[grep('r123_n4', SOX2$name)]$pattern <- 'R123, 4MM'
SOX2[grep('r3_n4', SOX2$name)]$pattern <- 'R3, 4MM'
SOX2$exp <- "yes"
SOX2[grep('r123_n1', SOX2$name)]$exp <- 'no'
SOX2[grep('r3_n2', SOX2$name)]$exp <- 'no'
SOX2[grep('r3_n3', SOX2$name)]$exp <- 'no'

#annotate sgRNA mismatches
SRSF1$pattern <- 'Region'
SRSF1[SRSF1$gene == 1]$pattern <- 'Full Match'
SRSF1[grep('r123_n1', SRSF1$name)]$pattern <- 'R123, 1MM'
SRSF1[grep('r1_n2', SRSF1$name)]$pattern <- 'R1, 2MM'
SRSF1[grep('r2_n2', SRSF1$name)]$pattern <- 'R2, 2MM'
SRSF1[grep('r3_n2', SRSF1$name)]$pattern <- 'R3, 2MM'
SRSF1[grep('r12_n2', SRSF1$name)]$pattern <- 'R12, 2MM'
SRSF1[grep('r1_n3', SRSF1$name)]$pattern <- 'R1, 3MM'
SRSF1[grep('r2_n3', SRSF1$name)]$pattern <- 'R2, 3MM'
SRSF1[grep('r3_n3', SRSF1$name)]$pattern <- 'R3, 3MM'
SRSF1[grep('r12_n3', SRSF1$name)]$pattern <- 'R12, 3MM'
SRSF1[grep('r123_n3', SRSF1$name)]$pattern <- 'R123, 3MM'
SRSF1[grep('r123_n4', SRSF1$name)]$pattern <- 'R123, 4MM'
SRSF1[grep('r3_n4', SRSF1$name)]$pattern <- 'R3, 4MM'
SRSF1$exp <- "yes"
SRSF1[grep('r123_n1', SRSF1$name)]$exp <- 'no'
SRSF1[grep('r3_n2', SRSF1$name)]$exp <- 'no'
SRSF1[grep('r3_n3', SRSF1$name)]$exp <- 'no'

#annotate sgRNA mismatches
CCND2$pattern <- 'Region'
CCND2[CCND2$gene == 1]$pattern <- 'Full Match'
CCND2[grep('r123_n1', CCND2$name)]$pattern <- 'R123, 1MM'
CCND2[grep('r1_n2', CCND2$name)]$pattern <- 'R1, 2MM'
CCND2[grep('r2_n2', CCND2$name)]$pattern <- 'R2, 2MM'
CCND2[grep('r3_n2', CCND2$name)]$pattern <- 'R3, 2MM'
CCND2[grep('r12_n2', CCND2$name)]$pattern <- 'R12, 2MM'
CCND2[grep('r1_n3', CCND2$name)]$pattern <- 'R1, 3MM'
CCND2[grep('r2_n3', CCND2$name)]$pattern <- 'R2, 3MM'
CCND2[grep('r3_n3', CCND2$name)]$pattern <- 'R3, 3MM'
CCND2[grep('r12_n3', CCND2$name)]$pattern <- 'R12, 3MM'
CCND2[grep('r123_n3', CCND2$name)]$pattern <- 'R123, 3MM'
CCND2[grep('r123_n4', CCND2$name)]$pattern <- 'R123, 4MM'
CCND2[grep('r3_n4', CCND2$name)]$pattern <- 'R3, 4MM'
CCND2$exp <- "yes"
CCND2[grep('r123_n1', CCND2$name)]$exp <- 'no'
CCND2[grep('r3_n2', CCND2$name)]$exp <- 'no'
CCND2[grep('r3_n3', CCND2$name)]$exp <- 'no'

#annotate sgRNA mismatches
CTNNB1$pattern <- 'Region'
CTNNB1[CTNNB1$gene == 1]$pattern <- 'Full Match'
CTNNB1[grep('r123_n1', CTNNB1$name)]$pattern <- 'R123, 1MM'
CTNNB1[grep('r1_n2', CTNNB1$name)]$pattern <- 'R1, 2MM'
CTNNB1[grep('r2_n2', CTNNB1$name)]$pattern <- 'R2, 2MM'
CTNNB1[grep('r3_n2', CTNNB1$name)]$pattern <- 'R3, 2MM'
CTNNB1[grep('r12_n2', CTNNB1$name)]$pattern <- 'R12, 2MM'
CTNNB1[grep('r1_n3', CTNNB1$name)]$pattern <- 'R1, 3MM'
CTNNB1[grep('r2_n3', CTNNB1$name)]$pattern <- 'R2, 3MM'
CTNNB1[grep('r3_n3', CTNNB1$name)]$pattern <- 'R3, 3MM'
CTNNB1[grep('r12_n3', CTNNB1$name)]$pattern <- 'R12, 3MM'
CTNNB1[grep('r123_n3', CTNNB1$name)]$pattern <- 'R123, 3MM'
CTNNB1[grep('r123_n4', CTNNB1$name)]$pattern <- 'R123, 4MM'
CTNNB1[grep('r3_n4', CTNNB1$name)]$pattern <- 'R3, 4MM'
CTNNB1$exp <- "yes"
CTNNB1[grep('r123_n1', CTNNB1$name)]$exp <- 'no'
CTNNB1[grep('r3_n2', CTNNB1$name)]$exp <- 'no'
CTNNB1[grep('r3_n3', CTNNB1$name)]$exp <- 'no'


#plot mismatch activity for each gene
SOX2_spec_control_plot <- ggplot(SOX2, aes(x=pattern, y = t12.activity, color = exp)) +
	geom_jitter( shape = 20, size = 0.6, width = 0.2) +
	theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
	ylab("") + xlab("") + scale_color_manual(values=c("red", "black")) + theme(legend.position="none")

SRSF1_spec_control_plot <- ggplot(SRSF1, aes(x=pattern, y = t12.activity, color = exp)) +
	geom_jitter( shape = 20, size = 0.6, width = 0.2) +
	theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
	ylab("") + xlab("") + scale_color_manual(values=c("red", "black")) + theme(legend.position="none")


CCND2_spec_control_plot <- ggplot(CCND2, aes(x=pattern, y = t12.activity, color = exp)) +
	geom_jitter( shape = 20, size = 0.6, width = 0.2) +
	theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
	ylab("") + xlab("") + scale_color_manual(values=c("red", "black")) + theme(legend.position="none")


CTNNB1_spec_control_plot <- ggplot(CTNNB1, aes(x=pattern, y = t12.activity, color = exp)) +
	geom_jitter( shape = 20, size = 0.6, width = 0.2) +
	theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
	ylab("") + xlab("") + scale_color_manual(values=c("red", "black")) + theme(legend.position="none")


spec_control <- as.data.frame(rbind(mcols(SOX2[,c('t12.activity', 'pattern')]),mcols(SRSF1[,c('t12.activity', 'pattern')]),mcols(CCND2[,c('t12.activity', 'pattern')])))

spec_control <- spec_control[spec_control$pattern != 'Region',]

levels(spec_control$pattern) <- unique(spec_control$pattern)

spec_control$pattern <- factor(spec_control$pattern, levels = c('Full Match', 'R123, 1MM', 'R1, 2MM', 'R2, 2MM', 'R3, 2MM', 'R12, 2MM', 'R1, 3 MM', 'R2, 3MM', 'R3, 3MM', 'R12, 3MM', 'R123, 3MM', 'R123, 4MM', 'R3, 4MM'))

spec_control <- spec_control[order(spec_control$pattern),]

spec_control <- spec_control[!is.na(spec_control$pattern),]

#plot cumulative mismatch activity across specificity control genes 
spec_control_plot <- ggplot(spec_control, aes(x=pattern, y = t12.activity)) +
	geom_boxplot(outlier.shape = NA) +
	theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
	ylab("Proportion of Full Match Phenotype") + xlab("")
	
