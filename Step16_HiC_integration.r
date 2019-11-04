library(GenomicInteractions)
library(rtracklayer)
library(GenomicFeatures)
library(ggbio)

#import high-resolution chromatin interactions generated from human neural precursor cells (NPC)
meta <- import("/Users/evangeller/Desktop/enhancer_analysis/hic_analysis/clean_master_hiccups_loops_nonsubsampled.bedpe", "bedpe")
npc_data <- makeGenomicInteractionsFromFile("/Users/evangeller/Desktop/enhancer_analysis/hic_analysis/clean_master_hiccups_loops_nonsubsampled.bedpe", type = "bedpe")
mcols(npc_data) <- mcols(meta)
#select NPC chromatin loop calls
npc_loops <- npc_data[npc_data$NA.15 == "NPC",]

#import proliferation phenotypes for genes and enhancers
S0X_exp <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_phenotype <- S0X_exp[(S0X_exp$t12.phenotype == "negative") | (S0X_exp$t12.phenotype == "positive"),]

#import gencode gene annotation
gencode <- import.gff("/Users/evangeller/Desktop/enhancer_analysis/bed/gencode.v19.annotation.gtf")
gencode_gene <- gencode[gencode$type == "gene"]
gencode_gene$id <- gencode_gene$gene_name

#select genes exhibiting a proliferation-altering phenotype following genetic disruption
gencode_gene_phenotype <- gencode_gene[gencode_gene$gene_name %in% S0X_phenotype$name,]

#define promoter regions for genes exhibiting proliferation phenotypes
gencode_gene_phenotype_promoter <- promoters(gencode_gene_phenotype , upstream = 10000, downstream = 10000)

expandRange = function(x, upstream=10000, downstream=10000) {
  strand_is_minus = strand(x) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(x)[on_plus] = start(x)[on_plus] - upstream
  start(x)[on_minus] = start(x)[on_minus] - downstream
  end(x)[on_plus] = end(x)[on_plus] + downstream
  end(x)[on_minus] = end(x)[on_minus] + upstream
  x
}
#extend promoter regions
gencode_gene_phenotype_promoter_gene <- expandRange(gencode_gene_phenotype)

#import enhancer-level summary of proliferation phenotypes
enhancer_set <- makeGRangesFromDataFrame(read.table( "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_enhancer_annotate_obs.txt", header= TRUE), keep.extra.columns= TRUE)
enhancer_set$id <- as.character(granges(enhancer_set))
ac_set <- enhancer_set[enhancer_set$class == 'ac',]
ac_set_phenotype <- ac_set[ac_set$n.phenotypes.obs > 0,]

#extend H3K27ac enhancer regions
start(ac_set_phenotype) <- start(ac_set_phenotype) - 10000
end(ac_set_phenotype) <- end(ac_set_phenotype) + 10000

#annotate high-resolution chromatin interactions based on overlap with enhancer or promoter regions that impact proliferation following genetic disruption
annotateInteractions(npc_data,list(enhancer = ac_set_phenotype, gene = gencode_gene_phenotype_promoter_gene ))
hcInt_enhancer_promoter <- npc_data[isInteractionType(npc_data, "enhancer", "gene") & (npc_data$NA.15 == "NPC")]
write.table(anchors(hcInt_enhancer_promoter), "/Users/evangeller/Desktop/enhancer_analysis/mageck/hNSC_hcInt_E-P.txt", row.names=F, sep=" ", quote= F)

#import manually processed file for downstream analysis of enhancer-promoter interactions
hcInt_processed <- read.table( "/Users/evangeller/Desktop/enhancer_analysis/mageck/hNSC_hcInt_E-P_processed.txt", header= F)

hcInt_enhancers <- ac_set_phenotype[ac_set_phenotype$id %in% hcInt_processed$V1,]
hcInt_results <- GRanges()

#make granges object from processed enhancer-promoter chromatin interactions
for(i in seq(1,nrow(hcInt_processed))) {
	enhancer <- hcInt_enhancers[hcInt_enhancers$id %in% hcInt_processed[i,]$V1,]
	enhancer$gene_name <- hcInt_processed[i,]$V2
	hcInt_results <- c(hcInt_results, enhancer)
	
}

#function association gene proliferation score with promoter identified in enhancer-promoter chromatin interactions
enhancer_gene_PC1 <- function(enhancer) {
	gene_id <- enhancer$gene_name
	gene_phenotype <- S0X_phenotype[S0X_phenotype$name %in% gene_id,]
	return(gene_phenotype$PC1)
}

#apply enhancer_gene_PC1 function to granges object of high-confidence enhancer-promoter chromatin interaction
mcols(hcInt_results)$gene.PC1 <- t(sapply(hcInt_results, enhancer_gene_PC1))

output <- as.data.frame(hcInt_results, row.names = NULL)

write.table(output, "/Users/evangeller/Desktop/enhancer_analysis/hic_analysis/Loop_call_E-P_int.txt", row.names=F, sep=" ", quote= F)

#import contact domains identified by arrowhead algorithm (Juicer software package from Aiden Lab) for regions containing a high-density of chromatin interactions
hNSC_contact_domain <- import.bed( "/Users/evangeller/Desktop/enhancer_analysis/hic_analysis/NPC_contact_domain.bed")

#subset proliferation-altering regions to gene targets
S0X_phenotype_gene <- S0X_phenotype[S0X_phenotype$gene == 1,]

#for each enhancer impact hNSC proliferation overlapping a contact domain, associate with proliferation-altering gene that also overlaps the contact domain
results <- GRanges()
for(i in seq_along(ac_set_phenotype)){
	local_domain <- hNSC_contact_domain[hNSC_contact_domain %over% ac_set_phenotype[i,]]
	local_gene_phenotypes <- gencode_gene_phenotype_promoter_gene[gencode_gene_phenotype_promoter_gene %over% local_domain]
	enhancer_gene_pair <- rep(ac_set_phenotype[i,], length(local_gene_phenotypes))
	enhancer_gene_pair$gene_name <- local_gene_phenotypes$id
	results <- c(results, enhancer_gene_pair)
}

#function association gene proliferation score with promoter identified by contact domains
enhancer_gene_PC1 <- function(enhancer) {
	gene_id <- enhancer$gene_name
	gene_phenotype <- S0X_phenotype[S0X_phenotype$name == gene_id,]
	return(gene_phenotype$PC1)
}

#annotate gene proliferation with enhancers associated to genes by contact domain interactions  
mcols(results)$gene.PC1 <- t(sapply(results, enhancer_gene_PC1))
output <- as.data.frame(results, row.names = NULL)
write.table(output, "/Users/evangeller/Desktop/enhancer_analysis/hic_analysis/Contact_domain_E-P_int.txt", row.names=F, sep=" ", quote= F)


#combine enhancer-promoter interactions identified by loop-calls and contact domains
full_results <- c(results, hcInt_results)
output <- as.data.frame(full_results, row.names = NULL)
write.table(output, "/Users/evangeller/Desktop/enhancer_analysis/hic_analysis/Loops_Contact_Domain_E-P_int.txt", row.names=F, sep=" ", quote= F)

S0X_enhancer <- S0X_exp[S0X_exp$enhancer == 1,]
S0X_enhancer_regions <- S0X_enhancer[S0X_enhancer %over% full_results,]
S0X_output <- as.data.frame(S0X_enhancer_regions)
write.table(S0X_output, file = "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA_HiC.txt", sep = " ", quote = FALSE, row.names = FALSE)


#plot gene proliferation score versus enhancer proliferation score for regulatory elements associated by chromatin interactions
full_results_neg <- full_results[(full_results$PC1.max.region < 0) & (full_results$gene.PC1 < 0), ]
full_results_evol <- full_results_neg[(full_results_neg$GAIN == 1) | (full_results_neg$HACNS_HAR == 1),]
ggplot(full_results_neg, aes(x=abs(gene.PC1), y= abs(PC1.max.region))) + geom_point(size = 1.2, alpha = 0.8, color = "#525252") + lims(x= c(0.5, 2.5), y = c(0.25,2.0)) +
	theme_classic() + geom_point(data = as.data.frame(full_results_evol), size = 1.2, color = '#ef3b2c') +
	theme(text = element_text(size=20), axis.title.x=element_blank(), axis.title.y=element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"))# get rid of legend panel bg
    
    
#plot beta scores for 4, 8, and 12 cell divisions for cortical enhancer interacting with multiple target genes
library(reshape2)
S0X_plot <- S0X_exp[(S0X_exp$name == "S05_phastCons_4084") | (S0X_exp$name == "PDS5B")  | (S0X_exp$name == "STARD13"),]
data_plot <- cbind(mcols(S0X_plot)$t4.beta,mcols(S0X_plot)$t8.beta,mcols(S0X_plot)$t12.beta)
colnames(data_plot) <- c("t4", "t8", "t12")
output <- melt(data_plot)
ggplot(output, aes(x = Var2, y = value, group = Var1)) + geom_point(size = 6.5, color = "#295ba9") + geom_line(linetype="dashed", color = "#525252", alpha = 0.6) + theme_classic() +
	theme(text = element_text(size=20),  axis.line = element_line(size = 1), axis.title.x=element_blank(), axis.title.y=element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    axis.text.x=element_blank(),
    axis.line.x=element_blank(),
    axis.ticks.x=element_blank(), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")) + # get rid of legend panel bg
    ylim(-1,0.25) + geom_hline(yintercept = 0)
    
#plot beta scores for 4, 8, and 12 cell divisions for human gain enhancer interacting with intellectual-disability associated gene
S0X_plot <- S0X_exp[(S0X_exp$name == "S05_phastCons_411") | (S0X_exp$name == "NSL1"),]
data_plot <- cbind(mcols(S0X_plot)$t4.beta,mcols(S0X_plot)$t8.beta,mcols(S0X_plot)$t12.beta)
colnames(data_plot) <- c("t4", "t8", "t12")
output <- melt(data_plot)
ggplot(output, aes(x = Var2, y = value, group = Var1)) + geom_point(size = 6.5, color = "#295ba9") + geom_line(linetype="dashed", color = "#525252", alpha = 0.6) + theme_classic() +
	theme(text = element_text(size=30),  axis.line = element_line(size = 3), axis.title.x=element_blank(), axis.title.y=element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    axis.text.x=element_blank(),
    axis.line.x=element_blank(),
    axis.ticks.x=element_blank(), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")) + # get rid of legend panel bg
    ylim(-1.4,0.25) + geom_hline(yintercept = 0, size = 3)
  

  
