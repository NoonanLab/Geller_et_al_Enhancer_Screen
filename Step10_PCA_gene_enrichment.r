library(ggbio)
library(rtracklayer)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ReactomePA)
library(org.Hs.eg.db)

#import biological effects following genetic disruption of genes and enhancers
S0X_exp <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/S0X_mageck_annotate_classified_PCA.txt", header= TRUE), keep.extra.columns= TRUE)

#select gene disruptions and filter for sign and timepoint irreproducibility
S0X_pCDS <- S0X_exp[S0X_exp$gene == 1]
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t4.phenotype != 'signIR']
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t8.phenotype != 'signIR']
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t12.phenotype != 'signIR']
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t4.phenotype != 'dynamic']
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t8.phenotype != 'dynamic']
S0X_pCDS <- S0X_pCDS[S0X_pCDS$t12.phenotype != 'dynamic']


S0X_promoter <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/mageck/archive/S0X_mageck_promoter_fpkm_constraint.txt", header= TRUE), keep.extra.columns= TRUE)


#generate partitions based on percentile of proliferation score for proliferation-decreasing phenotypes
S0X_pCDS_neg <- S0X_pCDS[S0X_pCDS$t12.phenotype == 'negative']
S0X_pCDS_pos <- S0X_pCDS[S0X_pCDS$t12.phenotype == 'positive']
pCDS_neg_quantiles<- quantile(S0X_pCDS_neg$PC1, probs = c(0.25, 0.5,1))
pCDS_pos_quantiles<- quantile(S0X_pCDS_pos$PC1, probs = c(0))
pCDS_categories <- c(pCDS_neg_quantiles, pCDS_pos_quantiles)
names(pCDS_categories) <- c('severe', 'strong', 'all', 'positive')

#import estimate of intolerance to loss-of-function mutations within genes
gene_constraint <- read.csv("/Users/evangeller/Desktop/enhancer_analysis/functional_constraint/gene_constraint_pLI.csv")
S0X_constraint <- makeGRangesFromDataFrame(inner_join(as.data.frame(S0X_pCDS), gene_constraint[c("gene", "pLI")], by=c("name"="gene")), keep.extra.columns= TRUE)
constrained_90 <- S0X_constraint[S0X_constraint$pLI >= 0.9]

microceph <- unique(c("MCPH1", "WDR62", "CDK5RAP2", "CASC5", "ASPM", "CENPJ", "STIL", "CEP135", "CEP152", "ZNF335", "PHC1", "CDK6","CENPE", "SASS6", "MFSD2A", "ANKLE2", "CIT", "WDFY3", "NBS1", "ATR", "XLF", "LIG4", "XRCC4", "Ku70/80", "DNA-PK", "XRCC2", "MCPH1", "PNKP", "ATR", "RBBP8", "CENPJ", "CEP152", "CEP63", "NIN", "DNA2", "TRAIP", "NSMCE2"))
asd <- unique(c("CHD8", "SCN2A", "ARID1B", "NRXN1", "SYNGAP1", "DYRK1A", "CHD2", "ANK2", "KDM5B", "ADNP", "POGZ","SUV420H1","SHANK2","TBR1","GRIN2B","DSCAM","KMT2C","PTEN","SHANK3","TCF7L2","TRIP12","SETD5","TNRC6B","ASH1L","CUL3","KATNAL2","WAC","NCKAP1","RANBP17","KDM6B","ILF2","DNMT3A","GIGYF1","WDFY3","SPAST","KAT2B","MYT1L","SLC6A1","PHF2","ZNF559","BCL11A","MFRP","FOXP1","GABRB3","MIB1","P2RX5","TRIO","KMT2E","ETFB","NINL","CTTNBP2","INTS6","USP45","CAPN12","OR52M1","DIP2A","AKAP9","NAA15","APH1A","IRF2BPL","ACHE","NLGN3","MBD5","PTK7","ERBB2IP"))
ddd <- unique(c("CEP290","FLNA","SLX4","ARG1","POLR3A","FANCB","NTRK1","MYH9","BBS9","TRIM32","PAX2","FGFR3","NSD1","ADAR","MFRP","FLNB","NAA10","AGPS","OFD1","PAFAH1B1","COL1A1","COL1A1","CHST3","HRAS","RPGRIP1","RIPK4","BRIP1","OFD1","TBX3","LAMA1","PCCA","FKRP","COX15","COG7","SCN2A","MECP2","BRAF","FTL","CBL","POMGNT1","COL2A1","TGFB3","PIGO","GJA1","PLK4","ARSA","CLN3","LRP5","SLC6A3","PITX3","CHAMP1","XPC","FOXE3","CLN8","RTEL1","MKKS","NBAS","NPHS1","ACAN","GJA8","FLNA","NAGA","NT5C3A",
"ARX","SMARCA2","IQSEC2","VSX2","CEP290","NSDHL","DPAGT1","KBTBD13","SOX9","TRAPPC2","HOXC13","CDC6","IFT140","PIGV","MSX1","MAGEL2","SLC2A1","TWIST2","HADH","PRRT2","GDF5","PAFAH1B1","OTC","CRYBB2","SALL1","PDE6G","NALCN","FOXN1","ATP7A","BBS10","MTO1","XYLT1","CHRNG","CHSY1","POMT2","NOG","KIF1BP","IFT172","CC2D2A","INPP5E","SATB2","GATA6","DOCK8","SLC46A1","DMPK","AGA","SLC25A38","COL11A2","ABCC6","AFF2","FLNB","HOXA1","GPSM2","ACTB","EBP","PSMB8","NALCN","LAMC3","PITX2","XRCC4","PEX13",
"SCN11A","COL9A1","NKX2-5","CRB1","TCF4","GFM1","TUBA8","IGF2","COL1A1","FGFR2","BBS4","PCNT","BUB1B","EYA1","NAGS","TBX1","RPS6KA3","RPS6KA3","SLC25A15","PYGL","ATP7A","ALDH18A1","PYCR1","GDI1","COQ2","GLDC","L1CAM","MPDU1","COL2A1","FANCC","CNTNAP2","MRE11","DDB2","FAT4","MYH9","FGFR3","PIK3R1","VSX2","VPS13B","ERCC3","SMPD1","GLUD1","B3GALT6","HADHA","FLNA","DDHD1","SLC52A3","HUWE1","DYNC1H1","MEGF10","ADGRG1","MPLKIP","L1CAM","WT1","RECQL4","B3GALT6","CHRDL1","DVL1","DSTYK","ALG8","MASP1",
"ROR2","NDP","KAT6A","BCOR","HDAC8","HDAC8","FREM2","HSPG2","ALG1","AUTS2","MYH9","RARS2","GNPTAB","SLC35C1","SATB2","NDUFS8","SCN4A","PEX5","WDR60","IKBKG","DLL3","RASA1","FANCF","FH","POMGNT1","RUNX2","CRYGC","OCRL","RFX6","LRPPRC","NPR2","SRCAP","COASY","NPHP4","EVC","SCN1A","COL9A1","ZFYVE26","DEAF1","WDR35","ELN","RNASEH2A","NR2F2","SMARCB1","MOCS1","PAX6","DLAT","SH3PXD2B","FLNA","NR2F1","KAT6B","GLI3","RETREG1","PAH","TGFBR1","FKTN","GHR","PAX9","FHL1","CDH23","TRPM1","PEX13","AP1S2","STS",
"ALPL","GLI2","ASAH1","MMAA","PEX14","FLNB","PGAP2","STRA6","ZNF711","PRSS12","MYH3","FKRP","HCCS","STAMBP","NKX3-2","TGFBR2","PEX16","DDC","NFU1","DEPDC5","PIK3CA","INPP5E","DHODH","PRPS1","GDF5","GLB1","DMP1","RAX","BBS5","SLC6A1","GAMT","SLC2A10","FLNA","ARID1B","FGFR2","BRAF","AGL","TRPS1","BMP4","KLHL40","MGAT2","RELN","GALT","UGT1A1","LYST","ZEB2","TTC8","PRKAR1A","IFT43","FLVCR2","CYP2U1","ANTXR1","OCRL","NPHS2","HOXD13","EHMT1","CDT1","SLC13A5","ERCC6","MANBA","PTDSS1","RSPO4","SLC25A20",
"WNT10B","PEX1","TREX1","KCNA2","KCNA2","RAD21","GNPAT","HEXB","MESP2","CTNS","PCBD1","TEK","BCAP31","HSD17B4","ANKRD11","HOXD13","NF1","OFD1","IKBKG","PTF1A","ALDH18A1","WNT7A","SCN4A","SMARCA4","SLC5A5","FANCG","TUBB2A","SUMF1","PEX6","CIB2","CEP41","NF1","HCN1","TGFB1","PTEN","FGFR1","SGSH","SATB2","GRIK2","ORC1","FANCE","ARL6","PTEN","IDUA","STXBP1","DNMT3A","UVSSA","ROBO3","BBS1","IDUA","ALX1","NSD1","BMP4","WNT3","MT-TP","ALG12","OXCT1","LRP5","NFIX","MFSD8","ERCC2","PDSS2","GUSB","HSD17B10",
"RIT1","ZFP57","ZMPSTE24","RTEL1","GDF5","PEX26","MLYCD","IGF2","HGSNAT","GCDH","PIGL","DNAAF4","AMPD2","ABCC9","KDM6A","ASPM","SHOX","PRRT2","GDF6","TMEM67","ZIC3","HRAS","KCNC1","MITF","MMADHC","GALK1","NYX","CRB2","PC","TMEM237","PEX12","GATA4","SCN8A","RXYLT1","PITX2","GABRB3","SCN2A","FOXF1","SOX2","ARX","PEX26","SURF1","GRIN2B","KIF11","CTCF","EXT2","RPGRIP1L","COL1A1","OPHN1","GJB2","ACADM","MYO5A","SMAD3","LRRC6","CASK","ETFDH","TFAP2B","DDX3X","DDX3X","SALL4","GRHL3","IRF6","HYLS1","ZDHHC9",
"DYM","FOXC1","NDUFV1","DOLK","ATRX","IFT122","NKX2-5","SLC6A5","WDR19","CEP152","DSPP","DSPP","CLCN7","OTOGL","PTHLH","GJB2","TBC1D24","MAFB","C12orf65","PKHD1","DCX","FAM111A","MPI","POC1B","MCCC2","PDHA1","FGFR1","CRX","CRYAA","GJC2","JAGN1","PMM2","SOX10","UMPS","VSX2","PEX19","SLC4A4","FOXE1","ARL6","ABCB11","GNAO1","LEMD3","TMEM67","FLVCR1","COL2A1","THRA","BFSP2","TUBGCP6","FLNA","NECTIN4","FRMD7","COL1A1","GLUL","ORC4","COLEC11","CHD7","CEP290","ZIC2","GRIN2A","TBC1D24","WDR45","ROR2","NEU1",
"COL1A1","CDH3","COL1A1","PAX3","NKX2-1","SNRPB","CTSA","HOXA1","BMPER","GJA1","SLC19A3","SOX17","GRIN2B","AGK","BCL11A","GNAS","GNS","WDR60","FOXC1","ARX","CRYBB1","CRYBB1","SLC6A8","PHGDH","PEX26","ALDH5A1","PEX7","NDUFS4","POMT1","PTHLH","C8orf37","GLB1","ARX","SLC39A13","BICD2","B4GALT7","DHFR","CA2","CHUK","HSF4","ZMPSTE24","FGFR2","FGFR1","GALC","RNASEH2B","ACOX1","DHCR24","GFAP","MAP2K1","L1CAM","TBCE","JAG1","SOX3","SCARF2","MTRR","TP63","TTC8","LAMA2","PLP1","PTS","COL10A1","FBN2","AAAS",
"IMPAD1","DNAAF3","DNMT3B","FMR1","SMC3","GJC2","SMOC1","GDI1","SLC2A2","CTNS","LRP2","ESCO2","CCBE1","DYNC2H1","RNU4ATAC","RAPSN","POGZ","WNT7A","TUBA1A","ADSL","COG1","GNAS","ELOVL4","TUBB","HOXD13","ZC4H2","ZC4H2","PDHA1","CLDN19","ALS2","ERCC6","SHH","ARX","TMEM67","TCF12","TAT","TTC19","SOX10","PURA","GJC2","MNX1","ADA","GRIN2B","ATP7A","GATA6","TFAP2A","TBXAS1","LHX4","EDA","COL11A1","SLC22A5","SCN4A","KCNJ11","PSPH","GLB1","UBE2A","FZD6","FGF3","ISPD","CCNQ","WRAP53","MAP3K1","KRAS","TBX5",
"MITF","MKS1","CRYBA4","FLNB","COL9A2","PEX3","PAX6","CRYAA","EXOSC3","FANCD2","AUH","PEX2","EXT1","DBT","BCKDHA","QDPR","KMT2D","HINT1","UBR1","GLI3","ALG6","HDAC8","YY1","SLC17A5","ARID1A","MED12","SRY","CRYBA1","TP63","SCN1B","RAI1","SETD5","DDHD2","TP63","MMP13","RPS19","GAA","PTH1R","PIK3R2","KCTD1","TBX20","TMEM165","CASK","FLNA","FAM126A","HNF1B","MOCS2","SYNGAP1","HR","PTEN","ATRX","IKBKG","IRF6","SLC12A6","TP63","SLC35A2","FTCD","KCNQ2","ASPA","EP300","UPF3B","SAMHD1","ASL","HOXD13","PEX3",
"COL2A1","SLC26A2","SIX5","COQ9","PITX2","PPT1","SLC9A6","CENPJ","NAA10","MICU1","NDUFS1","EPG5","LRP4","SIX1","SETBP1","KRAS","WDR34","ETFA","CBS","PEX12","GDF5","PIGA","NIPBL","KCNQ1","EVC","BRAF","PCARE","WNT1","SIK1","CA8","ELAC2","FKRP","PTPN11","TWIST1","SOX10","EOGT","FOXE3","GUCY2C","OBSL1","SMARCA2","FLNA","RARB","RARB","MMAB","RAB39B","MTHFR","HMGCL","PIK3R1","ERCC6","ALDOA","NAGLU","PDHX","PITX3","FGFR3","NOTCH2","PAX6","PEX10","GDF5","TWIST1","ARSB","PIK3CA","NHS","ROR2","COX7B","FGFR1",
"PEPD","WDR35","SLC26A2","NGLY1","SC5D","COL11A2","ECEL1","ALMS1","GATAD2B","ARMC4","RYR1","PEX10","RECQL4","CIB2","FGF10","KCNT1","ABHD5","GABRB3","MATN3","MECP2","FAR1","NDUFS7","CTNS","GMPPB","CYC1","ANKH","PEX6","TRIM32","GNPTG","MECP2","TTC7A","PGK1","NAGA","FANCA","FAM20A","ACY1","EXT1","NOG","PROP1","GTPBP3","PTCHD1","MED12","ITGA3","DPAGT1","CC2D2A","DDR2","LARGE1","MYO5B","SLC16A2","BHLHA9","HSD17B4","HPRT1","GLI3","KCTD7","DAG1","COL11A2","KCNT1","ABCB7","PTH1R","CTC1","FGFR2","FAH","CEP290",
"SMC1A","CENPJ","COMP","PRPS1","CHST14","CRYBB3","IHH","CDON","DYNC2H1","FGFR3","LAMA1","FUCA1","CACNA1C","ELN","TRIP11","HOXD13","COL4A3","IDS","SPRED1","TGIF1","BCS1L","AIPL1","PEX5","PHF6","CCDC114","KIF1A","FAM20C","GRM6","RAF1","KIF1A","KIF1A","SRD5A3","IGF1R","IGF1R","MID1","PNPT1","PLCE1","VLDLR","MEF2C","MSX2","ALDH7A1","WT1","NMNAT1","TUBB4A","POMT2","AK2","HSD3B7","EGR2","PTPN11","CC2D2A","WDR19","IGSF1","SOX10","ERCC6L2","CEP57","HPSE2","IFT172","KIF7","LMBRD1","RMRP","LAMP2","HNF4A","LTBP2",
"COL1A1","TP63","SLC26A2","OTX2","FGFR2","NSDHL","ERCC6","KIF7","GATA2","SDHAF1","FBN1","FBN1","ATP8B1","MAN2B1","FOXG1","ADNP","PEX2","MECP2","RPGRIP1","GCH1","TRPV4","LIG4","MYH9","MYH9","POLR1C","GNAS","ARID1B","IFIH1","SIX3","CHRNA4","PCDH19","TBX4","MSX2","LTBP2","MKS1","PLP1","GNPTAB","HEXA","ATP6V1B1","RSPH1","SMAD4","CPLANE1","ALX4","SCO2","LARGE1","MAF","CCNO","ERCC2","RAPSN","NALCN","WDPCP","ALX3","FREM1","CKAP2L","PEX14","MITF","UBE3A","POLR3B","WNT5A","NF1","POMT2","UBE3B","PTEN","JAK3",
"COX10","PEX2","ALG3","SLC17A5","NDE1","KANSL1","COL11A2","WDR11","ALDH1A3","SOS1","HAX1","ITGA7","TSPAN7","CHD7","GCH1","RPE65","CUL4B","RASA1","TMCO1","FGFR2","NOG","PEX19","FOLR1","PAPSS2","ASAH1","SDHA","TSHR","COL9A3","FGFR3","GRIA3","PHIP","ERCC5","TSC2","PTH1R","NPHP1","HSD17B10","VDR","PRSS56","GRIN1","NF1","PIGT","SPG11","TGFBR1","HPRT1","TYRP1","PAK3","FBN1","CSPP1","LIG4","PAX6","FLNB","CCND2","TTC37","P3H1","SOX2","PAX6","GNAI3","STXBP1","ALDH4A1","TBCE","MTM1","CTDP1","LFNG","EVC2","MAN1B1",
"PIK3CA","PNPT1","ACAD9","SPR","POMT1","MCEE","ATIC","AP4E1","ESCO2","NDUFS4","CPS1","INPPL1","MC2R","GLB1","DYRK1A","POMGNT2","COL4A3BP","GNAS","AKT1","COL2A1","PITX2","PTF1A","BCKDHB","TAZ","ALX4","MTR","MYCN","ACP5","FTSJ1","GBA2","IL1RAPL1","RBM8A","KCNJ10","SHH","PLOD2","COL18A1","COG4","TYR","LRP5","IGHMBP2","TUBA1A","KLF1","MLC1","CRB1","NOG","TGFBR1","PTEN","FGFR3","SPEG","HIBCH","SCO1","LHX3","UROS","HSPG2","IDUA","NEB","L2HGDH","SETBP1","SMPD1","RAB3GAP2","TGDS","LBR","SEC23B","MTOR","SOX10",
"HOXA13","SOX3","SLC27A4","ALDH3A2","PAX3","TBX22","HPS1","GK","TUBB2B","RPGRIP1L","CCDC103","HR","CLPB","LRP5","DMD","EFTUD2","NOG","PITX3","MMP13","POLR3B","APTX","MCCC1","SLC35D1","SLC26A2","NPHP1","PMS2","TAB2","PEX16","LTBP3","PHGDH","COL11A2","GJA8","IFT80","FKBP14","FOXC2","DLD","FGFR3","CRYGD","PAX8","PSPH","ENPP1","PDE4D","EVC2","PTEN","TRAPPC9","CTNNB1","KDM5C","FLT4","PEX1","MPV17","TSC1","PNKP","FGFR2","L1CAM","MKKS","SKI","SHOX","CDH23","DCHS1","POLR1D","GDF6","SURF1","PEX1","PHOX2B","BRWD3","WDR62","PEX10","PPP2R5D","ERCC8","DKC1","NPC1","GJA3","FLNA","DPM1","DLG3","STAR","GALE","GLI3","MMACHC","FYCO1","BBS2","APOPT1","FBXL4","FBP1","TUSC3","PCYT1A","NPHP1","MAP2K2","PRPS1","GDF5","FGFR3","RAB23","TMPRSS6","GATA6","HMGCS2","ASS1","POLG","VPS33B","RNASET2","MGP","TCOF1","ERCC2","EDA","HDAC4","PEX26","MYH3","HOXD13","SYNGAP1","ROGDI","ERCC3","GJB2","NRAS","NBN","LMX1B","POMT1","COL2A1","DLD","ACAT1","TSHB","GJB2","COX15","NPHP3","DIS3L2","GTF2H5","COL9A2","HYAL1","HCFC1","IHH","NKX2-5","CEP83","SHH","WAC","CYP1B1","SF3B4","HNRNPU","BTD","SYP","UROC1","IVD","SIL1","PORCN","TGFBR2","FKTN","BBS7","ETHE1","CTSK","ATM","PRPS1","TMEM70","MECP2","PHOX2B","PGM1","NEK1","NEK1","KCTD7","FOXP1","SMARCAL1","CHD7","DHCR7","PNKP","RSPH3","LARP7","KMT2A","KCNJ11","KCNJ11","PHF8","ARSE","TP63","MAF","BIN1","DARS","NPHP3","FKTN","CDKN1C","CDH3","CHD2","CHM","IKBKG","TSHR","PTCH1","ERCC1","POC1A","BMPR1B","PEX5","MAF","MUT","TMEM67","GRIN2A","FANCI","SHH","GLIS3","PPP2R1A","ADAR","ADAR","MEGF8","SKIV2L","BBS12","ODAPH","SMARCA4","AFF4","IFITM5","NPC2","NDUFS1","EIF4A3","CRYGD","KCNQ2","MAB21L2","MAB21L2","CEP152","SLC33A1","FGD1","HPGD","COL2A1","PEX7","COQ8A","XPA","NODAL","BLM","SIX1","NKX2-1","AMT","CUL7","CTSD","RPGRIP1L","SCN1B","RAB3GAP1","EDNRB","SLC2A1","COMP","SDCCAG8","NPHP3","NFIX","PTCH1","ACAN","PQBP1","MCPH1","SPAG1","GLE1","CCDC39","NUBPL","NHS","MYH9","MYO5A","EZH2","PTH1R","GLMN","FOXC2","COL4A4","TK2","ETFB","GJA1","TXNL4A","ERCC4","COL2A1","ORC6","POU1F1","FOXC1","PDHA1","COL11A1","ERCC1","EFNB1","CC2D1A","FKTN","GM2A","MCOLN1","ERF","EDNRA","PEX7","DYM","TRIM37","AKR1D1","DMD","PGAP3","ACADVL","PCCB","FOXP3","HSF4","TP63","BSND","COX6B1","SMARCB1","PALB2","SOX9","TBC1D24","WDR34","PGM3","IKBKG","SUCLG1","LRP5","PAX6","PAK3","FANCI","CSTB","DSPP","TH","CCDC65","POMGNT1","FBN1","RNASEH2C","CCDC40","GUCY2C","COX10","COL4A3","GDF5","SBDS","RTTN","IGF1","TRPV4","FRAS1","SNX14","PDGFRB","NSD1","SALL4","FGFR2","GATM","AHDC1","DDOST","ENPP1","RAB18","HNF4A","NDUFS4","PAX6","MFRP","ANKH","CDKL5","KCNB1","PSAP","GMPPA","TPP1","GJA1","PIEZO2","ZBTB20","CASK","CEP290","VIPAS39","SCN8A","DDX11","TSEN54","KAT6B","FOXRED1","GJB2","GALNS","HCFC1","COG8","CRYBB2","PLOD1","POC1A","CLN8","DMD","HYDIN","TCTN3","KIF22","NDUFA1","ASXL1","CLN5","DYNC1H1","LEMD3","DCX","CREBBP","GPC3","SHOC2","SLC39A13","ZIC3","TSC2","EIF2AK3","ACTG1","GDF5","RECQL4","PAH","HLCS","USB1","HACE1","KIAA0586","MMP21","PRMT7","AHI1","TUBB","MAPRE2","SLC39A8","TAF1","DKC1","TERT","PARN","TINF2","TCN2","COL6A3","COL6A3","ADGRG6","ZIC1","SPATA5","PDGFRB","DLL4","ALDH18A1","GAS8","NUP107","SLC25A26","RTN4IP1","CEP104","MFSD2A","PRDM12","CNOT3","CSNK2A1","GNAI1","KCNQ3","MSL3","MYT1L","PPM1D","PUF60","QRICH1","KMT5B","TCF20","ZBTB18","WAC","SCN2A","UNC80","CCDC115","TANGO2","DVL3","FGFR1","TBCK","RERE","ITPR1","ITPR1","GNB1","FLAD1","SMO","COL6A1","CDC45","TMEM126B","PPP1CB","ATAD3A","GORAB","UBA5","TBX15","CHD4","CDK13","PRKD1","CREBBP","PPA2","NAA10","NAA10","PKD1L1","IL11RA","OTULIN","NANS","PIEZO2","BRPF1","HIVEP2","STAG1","SON","FGF12","EBF3","TBCD","HECW2","TRIP12","THOC6","CFAP410","BHLHA9","WDR26","ERF","ARMC9","CWC27","NACC1","PEX11B","FGFR1","TBC1D23","SLC25A24","CAD","SMC1A","RPL11","EDAR","CDKN1C","UFM1","UFC1","STAG2","SLC52A2"))


#select each category of genomic feature for overrepresentation testing
microceph_grange <- S0X_pCDS[S0X_pCDS$name %in% microceph]
macroceph_grange <- S0X_pCDS[S0X_pCDS$name %in% macroceph]
asd_grange <- S0X_pCDS[S0X_pCDS$name %in% asd]
ddd_grange <- S0X_pCDS[S0X_pCDS$name %in% ddd]
DDD_cnv <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/DDD_CNV.txt", header= TRUE), keep.extra.columns= TRUE)
DDD_cnv_rep <- DDD_cnv[DDD_cnv$grade == 1]
CNV_pCDS_DDD <- S0X_pCDS[S0X_pCDS %over% DDD_cnv_rep]
ASD_cnv <- makeGRangesFromDataFrame(read.table("/Users/evangeller/Desktop/enhancer_analysis/asd_denovo/Sanders_dnCNV.txt", header= TRUE), keep.extra.columns= TRUE)
CNV_pCDS_ASD <- S0X_pCDS[S0X_pCDS %over% ASD_cnv]

#function for testing overrepresentation using the hypergeometric distribution
pCDS_hyperGeometric <- function(category, quantiles){
	overlap_features <- S0X_pCDS[S0X_pCDS %over% category]
	nonoverlap_features <- S0X_pCDS[!S0X_pCDS %over% category]
	result <- rep(c('NA'), length(quantiles))
	names(result) <- names(quantiles)
	for(phenotype in c('severe', 'strong', 'all')) {
		success_in_sample <- length(overlap_features[(overlap_features$PC1 <= quantiles[phenotype]) & (overlap_features$t12.phenotype == 'negative')])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		failure_in_bkgd <- length(S0X_pCDS$t12.phenotype == 'neutral') - success_in_bkgd
		sample_size <- length(S0X_pCDS[(S0X_pCDS$PC1 <= quantiles[phenotype]) & (S0X_pCDS$t12.phenotype == 'negative')])
		result[phenotype] <- phyper(success_in_sample -1, success_in_bkgd, failure_in_bkgd, sample_size, lower.tail= FALSE);
	}
	for(phenotype in c('positive')) {
		success_in_sample <- length(overlap_features[overlap_features$t12.phenotype == 'positive'])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		failure_in_bkgd <- length(S0X_pCDS$t12.phenotype == 'neutral') - success_in_bkgd
		sample_size <- length(S0X_pCDS[S0X_pCDS$t12.phenotype == 'positive'])
		result[phenotype] <- phyper(success_in_sample -1, success_in_bkgd, failure_in_bkgd, sample_size, lower.tail= FALSE);
	}
	return(result)
}

#perform hyergeometric overrepresentation test for each category
Constraint_hyper <- pCDS_hyperGeometric(constrained_90,pCDS_categories)
Microceph_hyper <- pCDS_hyperGeometric(microceph_grange,pCDS_categories)
Macroceph_hyper <- pCDS_hyperGeometric(macroceph_grange,pCDS_categories)
ASD_hyper <- pCDS_hyperGeometric(asd_grange,pCDS_categories)
DDD_hyper <- pCDS_hyperGeometric(ddd_grange,pCDS_categories)
ASD_CNV_hyper <- pCDS_hyperGeometric(CNV_pCDS_ASD,pCDS_categories)
DDD_CNV_hyper <- pCDS_hyperGeometric(CNV_pCDS_DDD,pCDS_categories)


pCDS_foldEnrich <- function(category, quantiles){
	overlap_features <- S0X_pCDS[S0X_pCDS %over% category]
	nonoverlap_features <- S0X_pCDS[!S0X_pCDS %over% category]
	result <- rep(c('NA'), length(quantiles))
	names(result) <- names(quantiles)
	for(phenotype in c('severe', 'strong', 'all')) {
		success_in_sample <- length(overlap_features[(overlap_features$PC1 <= quantiles[phenotype]) & (overlap_features$t12.phenotype == 'negative')])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		sample_bkgd <- length(S0X_pCDS$t12.phenotype == 'neutral')
		sample_size <- length(S0X_pCDS[(S0X_pCDS$PC1 <= quantiles[phenotype]) & (S0X_pCDS$t12.phenotype == 'negative')])
		result[phenotype] <- (success_in_sample / sample_size) / (success_in_bkgd / sample_bkgd)
	}
	for(phenotype in c('positive')) {
		success_in_sample <- length(overlap_features[overlap_features$t12.phenotype == 'positive'])
		success_in_bkgd <- length(overlap_features[overlap_features$t12.phenotype == 'neutral'])
		sample_bkgd <- length(S0X_pCDS$t12.phenotype == 'neutral')
		sample_size <- length(S0X_pCDS[S0X_pCDS$t12.phenotype == 'positive'])
		result[phenotype] <- (success_in_sample / sample_size) / (success_in_bkgd /sample_bkgd)
	}
	return(result)
}

#calculate fold-enrichment relative to neutral gene genetic disruptions for each category
Constraint_foldEnrich <- as.numeric(pCDS_foldEnrich(constrained_90,pCDS_categories))
Microceph_foldEnrich <- as.numeric(pCDS_foldEnrich(microceph_grange,pCDS_categories))
ASD_foldEnrich <- as.numeric(pCDS_foldEnrich(asd_grange,pCDS_categories))
DDD_foldEnrich <- as.numeric(pCDS_foldEnrich(ddd_grange,pCDS_categories))
ASD_CNV_foldEnrich <- as.numeric(pCDS_foldEnrich(CNV_pCDS_ASD,pCDS_categories))
DDD_CNV_foldEnrich <- as.numeric(pCDS_foldEnrich(CNV_pCDS_DDD,pCDS_categories))

#prepare gene names for pathway enrichment analysis using Reactome database
background_genes <- S0X_pCDS$name
all_phenotype_genes <- S0X_pCDS[(S0X_pCDS$PC1 <= pCDS_categories['all']) & (S0X_pCDS$t12.phenotype == 'negative')]$name
strong_phenotype_genes <- S0X_pCDS[S0X_pCDS$PC1 <= pCDS_categories['strong']& (S0X_pCDS$t12.phenotype == 'negative')]$name
severe_phenotype_genes <- S0X_pCDS[S0X_pCDS$PC1 <= pCDS_categories['severe']& (S0X_pCDS$t12.phenotype == 'negative')]$name
positive_phenotype_genes <- S0X_pCDS[S0X_pCDS$t12.phenotype == 'positive']$name
hNSC_phenotype_genes <- hNSC_grange[hNSC_grange$t12.phenotype == 'negative']$name

#convert gene names to entrez id
background_id = bitr(background_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
all_id = bitr(all_phenotype_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
strong_id = bitr(strong_phenotype_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
severe_id = bitr(severe_phenotype_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
positive_id = bitr(positive_phenotype_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
hNSC_id = bitr(hNSC_phenotype_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#perform enrichment analyses of signaling pathways defined in the Reactome database
reactome_all <- enrichPathway(gene=all_id$ENTREZID, universe = background_id$ENTREZID,pvalueCutoff=0.05, qvalueCutoff = 0.05,readable=T)
write.csv(reactome_all , file = "/Users/evangeller/Desktop/enhancer_analysis/Reactome/all_reactome.csv", row.names=FALSE)
reactome_strong <- enrichPathway(gene=strong_id$ENTREZID, universe = background_id$ENTREZID,pvalueCutoff=0.05, readable=T)
write.csv(reactome_strong , file = "/Users/evangeller/Desktop/enhancer_analysis/Reactome/strong_reactome.csv", row.names=FALSE)
reactome_severe <- enrichPathway(gene=severe_id$ENTREZID, universe = background_id$ENTREZID,pvalueCutoff=0.05, qvalueCutoff = 0.05, readable=T)
write.csv(reactome_severe , file = "/Users/evangeller/Desktop/enhancer_analysis/Reactome/severe_reactome.csv", row.names=FALSE)
reactome_positive <- enrichPathway(gene=positive_id$ENTREZID, universe = background_id$ENTREZID,pvalueCutoff=0.05, readable=T)
write.csv(reactome_positive , file = "/Users/evangeller/Desktop/enhancer_analysis/Reactome/positive_reactome.csv", row.names=FALSE)
reactome_full <- enrichPathway(gene=c(all_id$ENTREZID,positive_id$ENTREZID), universe = background_id$ENTREZID,pvalueCutoff=0.05, qvalueCutoff = 0.05,readable=T)
write.csv(reactome_full , file = "/Users/evangeller/Desktop/enhancer_analysis/Reactome/full_reactome.csv", row.names=FALSE)

#generate histogram plot of proliferation-altering regions based on PC1 'proliferation score'
S0X_output <- as.data.frame(S0X_phenotype)
hist_plot <- ggplot(S0X_output, aes(x=PC1, y = ..count.. + 1)) 
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_output, gene == 1), binwidth = 0.05, fill = '#969696', color = 'black')
hist_plot <- hist_plot + geom_histogram(data = subset(S0X_output, shuffle_control == 1), binwidth = 0.05, fill = 'lightgrey', color = 'black')
hist_plot <- hist_plot + scale_y_log10() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
hist_plot <- hist_plot + coord_fixed(ratio = 0.5)
