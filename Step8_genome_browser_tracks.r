library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

#import full dataset of phenotyped biological effects
S0X_batch <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA.txt", header= TRUE), keep.extra.columns= TRUE)

strand(S0X_batch) <- "*"

#discard regions that are not enhancer targeting or genomic background controls
S0X_df <- S0X_batch[(S0X_batch$spec_control == 0)&(S0X_batch$gene == 0)]


strand(S0X_df) <- "*"

S0X_exp <- S0X_df[!duplicated(S0X_df$name)]


seqlevels(S0X_exp) <- seqlevels(Hsapiens)[1:23]
genome(S0X_exp) <- genome(Hsapiens)[1:23]
seqlengths(S0X_exp) <- seqlengths(Hsapiens)[1:23]


t4_export <- granges(S0X_exp)
t8_export <- granges(S0X_exp)
t12_export <- granges(S0X_exp)

t4_export$score <- S0X_exp$t4.z.score
t8_export$score <- S0X_exp$t8.z.score
t12_export$score <- S0X_exp$t12.z.score

#export bigWig files for NTC-normalized Z-scores at each timepoint
export(t4_export, '~/Desktop/enhancer_analysis/genome_browser/Beta_KO_Z_score_t4.bigWig', format='bigWig')
export(t8_export, '~/Desktop/enhancer_analysis/genome_browser/Beta_KO_Z_score_t8.bigWig', format='bigWig')
export(t12_export, '~/Desktop/enhancer_analysis/genome_browser/Beta_KO_Z_score_t12.bigWig', format='bigWig')


start(S0X_exp) <- start(S0X_exp) - 1
end(S0X_exp) <- end(S0X_exp) + 1


output <- as.data.frame(S0X_exp)

S0X_enhancer <- S0X_exp[S0X_exp$enhancer == 1]


#prepare colorized phenotype annotation of conserved regions and genome background controls
S0X_enhancer_neg <- S0X_enhancer[S0X_enhancer$t12.phenotype == 'negative']

S0X_enhancer_pos <- S0X_enhancer[S0X_enhancer$t12.phenotype == 'positive']

#partition proliferation-decreasing phenotypes by percentile of biological effect
enhancer_categories <- c(quantile(S0X_enhancer_neg$PC1, probs = c(0.25, 0.5, 1.0)),0)


names(enhancer_categories) <- c('severe', 'strong','all', 'positive')

output$strand <- "."
output$score <- "0"

#assign hex color code to proliferation-decreasing and proliferation-increasing phenotypes
output$t4.rgb <- '217,217,217'
output[output$t4.phenotype == 'negative' ,]$t4.rgb <- '41,91,169'
output[output$t4.phenotype == 'positive' ,]$t4.rgb <- '0,141,75'

output$t8.rgb <- '217,217,217'
output[output$t8.phenotype == 'negative' ,]$t8.rgb <- '41,91,169'
output[output$t8.phenotype == 'positive' ,]$t8.rgb <- '0,141,75'

output$t12.rgb <- '217,217,217'
output[output$t12.phenotype == 'negative' ,]$t12.rgb <- '41,91,169'
output[output$t12.phenotype == 'positive' ,]$t12.rgb <- '0,141,75'

#assign hex color code to proliferation-altering regions based on percentile of biological effect
output$PCA.rgb <- '217,217,217'
output[(output$t12.phenotype == 'negative')  & (output$PC1 <= enhancer_categories['all']),]$PCA.rgb <- '2,56,88'
output[(output$t12.phenotype == 'negative') & (output$PC1 <= enhancer_categories['strong']),]$PCA.rgb <- '4,90,141'
output[(output$t12.phenotype == 'negative') & (output$PC1 <= enhancer_categories['severe']),]$PCA.rgb <- '5,112,176'
output[output$t12.phenotype == 'positive' ,]$PCA.rgb <- '0,141,75'


write.table(output[,c('seqnames','start','end','name','score','strand','start','end','t4.rgb')], "~/Desktop/enhancer_analysis/genome_browser/Beta_classifier_phenotype_t4.bed", sep="\t", quote=F, col.names=F, row.names=F)

write.table(output[,c('seqnames','start','end','name','score','strand','start','end','t8.rgb')], "~/Desktop/enhancer_analysis/genome_browser/Beta_classifier_phenotype_t8.bed", sep="\t", quote=F, col.names=F, row.names=F)

write.table(output[,c('seqnames','start','end','name','score','strand','start','end','t12.rgb')], "~/Desktop/enhancer_analysis/genome_browser/Beta_classifier_phenotype_t12.bed", sep="\t", quote=F, col.names=F, row.names=F)

write.table(output[,c('seqnames','start','end','name','score','strand','start','end','PCA.rgb')], "~/Desktop/enhancer_analysis/genome_browser/Beta_classifier_phenotype_PCA.bed", sep="\t", quote=F, col.names=F, row.names=F)


