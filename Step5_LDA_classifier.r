library(MASS)
library(ggplot2)
library(rtracklayer)

#import jointly modeled biological effects following genetic disruption
S0X_exp <- as.data.frame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_spec.txt", header= TRUE))

#log transform permutation-based p-value
S0X_exp$log10p.t4 <- -log10(S0X_exp$t4.p.value) 
S0X_exp$log10p.t8 <- -log10(S0X_exp$t8.p.value)
S0X_exp$log10p.t12 <- -log10(S0X_exp$t12.p.value)

#compute 'crispr score'; composite of biological effect and log transformed p-value
S0X_exp$CS.t4 <- S0X_exp$log10p.t4 * S0X_exp$t4.beta
S0X_exp$CS.t8 <- S0X_exp$log10p.t8 * S0X_exp$t8.beta
S0X_exp$CS.t12 <- S0X_exp$log10p.t12 * S0X_exp$t12.beta

S0X_exp$t4.train.phenotype <- 'NA'
S0X_exp$t8.train.phenotype <- 'NA'
S0X_exp$t12.train.phenotype <- 'NA'

#assign proliferation-decreasing controls 'negative' biological effect
S0X_exp[S0X_exp$neg_control == 1,]$t4.train.phenotype <- 'negative'
S0X_exp[S0X_exp$neg_control == 1,]$t8.train.phenotype <- 'negative'
S0X_exp[S0X_exp$neg_control == 1,]$t12.train.phenotype <- 'negative'

#import rep 1, rep2 biological effects
S0X_rep1 <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_rep1.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_rep2 <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_rep2.txt", header= TRUE), keep.extra.columns= TRUE)

#annotate sign irreproducible regions across biological replicates
S0X_exp$signIR.t4 <- sign(S0X_rep1$t4.beta) != sign(S0X_rep2$t4.beta)
S0X_exp$signIR.t8 <- sign(S0X_rep1$t8.beta) != sign(S0X_rep2$t8.beta)
S0X_exp$signIR.t12 <- sign(S0X_rep1$t12.beta) != sign(S0X_rep2$t12.beta)

#assign top 100 regions from each timepoint as 'positive' biological effects
S0X_exp[order(-S0X_exp$t4.beta),][1:1000,][S0X_exp[order(-S0X_exp$t4.beta),][1:1000,]$signIR.t4 == TRUE,][1:100,]$t4.train.phenotype <-'positive'
S0X_exp[order(-S0X_exp$t8.beta),][1:1000,][S0X_exp[order(-S0X_exp$t8.beta),][1:1000,]$signIR.t8 == TRUE,][1:100,]$t8.train.phenotype <-'positive'
S0X_exp[order(-S0X_exp$t12.beta),][1:1000,][S0X_exp[order(-S0X_exp$t12.beta),][1:1000,]$signIR.t12 == TRUE,][1:100,]$t12.train.phenotype <-'positive'

#assign genomic backgroudn controls as 'neutral' biological effects
S0X_exp[S0X_exp$shuffle_control == 1,]$t4.train.phenotype <- 'neutral'
S0X_exp[S0X_exp$shuffle_control == 1,]$t8.train.phenotype <- 'neutral'
S0X_exp[S0X_exp$shuffle_control == 1,]$t12.train.phenotype <- 'neutral'

t4.train <- c(grep('negative', S0X_exp$t4.train.phenotype),grep('neutral', S0X_exp$t4.train.phenotype), grep('positive', S0X_exp$t4.train.phenotype))
t8.train <- c(grep('negative', S0X_exp$t8.train.phenotype),grep('neutral', S0X_exp$t8.train.phenotype), grep('positive', S0X_exp$t8.train.phenotype))
t12.train <- c(grep('negative', S0X_exp$t12.train.phenotype),grep('neutral', S0X_exp$t12.train.phenotype), grep('positive', S0X_exp$t12.train.phenotype))

tX.train.t4 <- S0X_exp[t4.train, c('CS.t4',  't4.train.phenotype')]
tX.train.t8 <- S0X_exp[t8.train, c('CS.t8', 't8.train.phenotype')]
tX.train.t12 <- S0X_exp[t12.train, c('CS.t12',  't12.train.phenotype')]

colnames(tX.train.t4) <- c('CS', 'train.phenotype')
colnames(tX.train.t8) <- c('CS','train.phenotype')
colnames(tX.train.t12) <- c('CS','train.phenotype')

tX.train <- rbind(tX.train.t4, tX.train.t8, tX.train.t12)

plot_t4 <- ggplot(tX.train.t4, aes(x=CS, fill=train.phenotype)) + geom_density(alpha = 0.8) + scale_fill_manual(values=c('#295ba9','#D1D3D4','#008d4b')) + theme_minimal()
plot_t8 <- ggplot(tX.train.t8, aes(x=CS, fill=train.phenotype)) + geom_density(alpha = 0.8) + scale_fill_manual(values=c('#295ba9','#D1D3D4','#008d4b')) + theme_minimal()
plot_t12 <- ggplot(tX.train.t12, aes(x=CS, fill=train.phenotype)) + geom_density(alpha = 0.8) + scale_fill_manual(values=c('#295ba9','#D1D3D4','#008d4b')) + theme_minimal()

plot_tX <- ggplot(tX.train, aes(x= CS, fill=train.phenotype)) + geom_density(alpha = 0.8) + scale_fill_manual(values=c('#295ba9','#D1D3D4','#008d4b')) + theme_minimal()

set.seed(31717)

#initialize linear-discriminant analysis (LDA) with uniform prior probability across classification groups
tX.lda <- lda(train.phenotype ~ CS, tX.train, prior = c(1,1,1)/3)

#classify all genetic disruptions at 4 cell divisions
S0X_exp$CS <- S0X_exp$CS.t4
tX.lda.p.t4 <- predict(object = tX.lda, newdata = S0X_exp)
S0X_exp$t4.predict <- tX.lda.p.t4$class
summary(S0X_exp$t4.predict )

#classify all genetic disruptions at 8 cell divisions
S0X_exp$CS <- S0X_exp$CS.t8
tX.lda.p.t8 <- predict(object = tX.lda, newdata = S0X_exp)
S0X_exp$t8.predict <- tX.lda.p.t8$class
summary(S0X_exp$t8.predict )

#classify all genetic disruptions at 12 cell divisions
S0X_exp$CS <- S0X_exp$CS.t12
tX.lda.p.t12 <- predict(object = tX.lda, newdata = S0X_exp)
S0X_exp$t12.predict <- tX.lda.p.t12$class
summary(S0X_exp$t12.predict )

#initialize proliferation phenotype assignment
S0X_exp$t4.phenotype <- 'neutral'
S0X_exp$t8.phenotype <- 'neutral'
S0X_exp$t12.phenotype <- 'neutral'

S0X_exp <- makeGRangesFromDataFrame(S0X_exp, keep.extra.columns=TRUE)

#assign biologically consistent proliferation phenotypes across multiple timepoints
S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict == 'negative') & (S0X_exp$t12.predict == 'negative') ]$t4.phenotype <- 'negative'
S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict == 'positive') & (S0X_exp$t12.predict == 'positive') ]$t4.phenotype <- 'positive'
S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict == 'negative') & (S0X_exp$t12.predict == 'negative') ]$t8.phenotype <- 'negative'
S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict == 'positive') & (S0X_exp$t12.predict == 'positive') ]$t8.phenotype <- 'positive'
S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict == 'negative') & (S0X_exp$t12.predict == 'negative') ]$t12.phenotype <- 'negative'
S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict == 'positive') & (S0X_exp$t12.predict == 'positive') ]$t12.phenotype <- 'positive'

S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict == 'negative') & (S0X_exp$t12.predict == 'negative') ]$t8.phenotype <- 'negative'
S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict == 'positive') & (S0X_exp$t12.predict == 'positive') ]$t8.phenotype <- 'positive'
S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict == 'negative') & (S0X_exp$t12.predict == 'negative') ]$t12.phenotype <- 'negative'
S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict == 'positive') & (S0X_exp$t12.predict == 'positive') ]$t12.phenotype <- 'positive'

S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict == 'neutral') & (S0X_exp$t12.predict == 'negative') ]$t12.phenotype <- 'negative'
S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict == 'neutral') & (S0X_exp$t12.predict == 'positive') ]$t12.phenotype <- 'positive'

#assign 'dynamic' for proliferation phenotypes not exhibiting biologicaly consistency
S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict != 'neutral') & (S0X_exp$t12.predict == 'neutral') ]$t4.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict != 'neutral') & (S0X_exp$t12.predict == 'neutral') ]$t8.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'neutral') & (S0X_exp$t8.predict != 'neutral') & (S0X_exp$t12.predict == 'neutral') ]$t12.phenotype <- 'dynamic'

S0X_exp[(S0X_exp$t4.predict != 'neutral') & (S0X_exp$t8.predict != 'neutral') & (S0X_exp$t12.predict == 'neutral') ]$t4.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict != 'neutral') & (S0X_exp$t8.predict != 'neutral') & (S0X_exp$t12.predict == 'neutral') ]$t8.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict != 'neutral') & (S0X_exp$t8.predict != 'neutral') & (S0X_exp$t12.predict == 'neutral') ]$t12.phenotype <- 'dynamic'

S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict != 'negative') & (S0X_exp$t12.predict == 'negative') ]$t4.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict != 'negative') & (S0X_exp$t12.predict == 'negative') ]$t8.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict != 'negative') & (S0X_exp$t12.predict == 'negative') ]$t12.phenotype <- 'dynamic'

S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict != 'negative') & (S0X_exp$t12.predict != 'negative') ]$t4.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict != 'negative') & (S0X_exp$t12.predict != 'negative') ]$t8.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'negative') & (S0X_exp$t8.predict != 'negative') & (S0X_exp$t12.predict != 'negative') ]$t12.phenotype <- 'dynamic'

S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict != 'positive') & (S0X_exp$t12.predict == 'positive') ]$t4.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict != 'positive') & (S0X_exp$t12.predict == 'positive') ]$t8.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict != 'positive') & (S0X_exp$t12.predict == 'positive') ]$t12.phenotype <- 'dynamic'

S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict != 'positive') & (S0X_exp$t12.predict != 'positive') ]$t4.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict != 'positive') & (S0X_exp$t12.predict != 'positive') ]$t8.phenotype <- 'dynamic'
S0X_exp[(S0X_exp$t4.predict == 'positive') & (S0X_exp$t8.predict != 'positive') & (S0X_exp$t12.predict != 'positive') ]$t12.phenotype <- 'dynamic'

#import rep 1 and rep2 biological effects
S0X_rep1 <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_rep1.txt", header= TRUE), keep.extra.columns= TRUE)
S0X_rep2 <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_rep2.txt", header= TRUE), keep.extra.columns= TRUE)


S0X_exp$signIR.t4 <- sign(S0X_rep1$t4.beta) != sign(S0X_rep2$t4.beta)
S0X_exp$signIR.t8 <- sign(S0X_rep1$t8.beta) != sign(S0X_rep2$t8.beta)
S0X_exp$signIR.t12 <- sign(S0X_rep1$t12.beta) != sign(S0X_rep2$t12.beta)

#annotate proliferation phenotypes exhibiting biological effect sign inconsistency 'signIR'
S0X_exp[(S0X_exp$t4.phenotype != "neutral") & (S0X_exp$signIR.t4 == TRUE)]$t4.phenotype <- 'signIR'
S0X_exp[(S0X_exp$t4.phenotype != "neutral") & (S0X_exp$signIR.t4 == TRUE)]$t8.phenotype <- 'signIR'
S0X_exp[(S0X_exp$t4.phenotype != "neutral") & (S0X_exp$signIR.t4 == TRUE)]$t12.phenotype <- 'signIR'

S0X_exp[(S0X_exp$t8.phenotype != "neutral") & (S0X_exp$signIR.t8 == TRUE)]$t4.phenotype <- 'signIR'
S0X_exp[(S0X_exp$t8.phenotype != "neutral") & (S0X_exp$signIR.t8 == TRUE)]$t8.phenotype <- 'signIR'
S0X_exp[(S0X_exp$t8.phenotype != "neutral") & (S0X_exp$signIR.t8 == TRUE)]$t12.phenotype <- 'signIR'

S0X_exp[(S0X_exp$t12.phenotype != "neutral") & (S0X_exp$signIR.t12 == TRUE)]$t4.phenotype <- 'signIR'
S0X_exp[(S0X_exp$t12.phenotype != "neutral") & (S0X_exp$signIR.t12 == TRUE)]$t8.phenotype <- 'signIR'
S0X_exp[(S0X_exp$t12.phenotype != "neutral") & (S0X_exp$signIR.t12 == TRUE)]$t12.phenotype <- 'signIR'

S0X_exp <- S0X_exp[S0X_exp$t12.phenotype != 'signIR']

#calculate the false positive rate and false discovery rate among the proliferation phenotype call set
train.positive <- S0X_exp[(S0X_exp$t4.train.phenotype == 'positive') | (S0X_exp$t8.train.phenotype == 'positive') | (S0X_exp$t12.train.phenotype == 'positive')]	
train.neutral <- S0X_exp[(S0X_exp$t4.train.phenotype == 'neutral') | (S0X_exp$t8.train.phenotype == 'neutral') | (S0X_exp$t12.train.phenotype == 'neutral')]	
train.negative <- S0X_exp[(S0X_exp$t4.train.phenotype == 'negative') | (S0X_exp$t8.train.phenotype == 'negative') | (S0X_exp$t12.train.phenotype == 'negative')]	

V <- length(train.neutral[(train.neutral$t12.phenotype == 'negative') | (train.neutral$t12.phenotype == 'positive')])
N <- length(train.neutral)
FPR <- V/N
S <- length(train.positive[train.positive$t12.phenotype == 'positive']) + length(train.negative[train.negative$t12.phenotype == 'negative'])
FDR <- V / (V+ S)

S0X_output <- subset(S0X_exp,select = -c(CS,signIR.t4, signIR.t8, signIR.t12))


write.table(S0X_output, file = "/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_filtered_spec.txt", sep = " ", quote = FALSE, row.names = FALSE)

