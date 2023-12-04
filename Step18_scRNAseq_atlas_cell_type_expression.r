
# Load in required R packages
library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(DoubletFinder)
library(dplyr)
library(rliger)
library(SeuratWrappers)
library(phateR)
library(stringr)
library(clipr)
library(plotly)
library(tictoc)
library(purrr)

sample.names <- list.files("./Eze_et_al_2021_analysis/cortex_datasets")

# Read in 10X cortex datasets
tic()
Eze_dataset_list <- lapply(sample.names, function(x) {
  Read10X(data.dir = str_c("./Eze_et_al_2021_analysis/cortex_datasets/", x))
})
names(Eze_dataset_list) <- sample.names
toc()

# Create individual Seurat objects from each cortex dataset
tic()
Eze_dataset_list <- lapply(Eze_dataset_list, function(x) {
  CreateSeuratObject(x, min.cells = 50, min.features = 750)
})
toc()

# Remove low quality datasets based on minimum filtering steps
Eze_dataset_list <- Eze_dataset_list[!names(Eze_dataset_list) %in% c("CS20_frontal1", "CS20_frontal2", "GW16_V1", "GW16_motor", "GW16_Somato",
                                                                     "GW15_temporal", "GW17_parietal", "GW15_parietal")]

# Calculate total percentage of mitochondrial gene expression per cell
tic()
Eze_dataset_list <- lapply(Eze_dataset_list, function(x) {
  AddMetaData(x, metadata = PercentageFeatureSet(x, pattern = "^MT-"), col.name = "percent.mt")
})
toc()

# Setup function for filtering cells by quality metrics
quality_metric_subset_function <- function(x) {
  
  mean_nCount <- mean(x$nCount_RNA)
  upper_nCount <- mean(x$nCount_RNA) + 2*sd(x$nCount_RNA);
  lower_nCount <- mean(x$nCount_RNA) - 2*sd(x$nCount_RNA);
  mean_nFeature <- mean(x$nFeature_RNA)
  upper_nFeature <- mean(x$nFeature_RNA) + 2*sd(x$nFeature_RNA);
  lower_nFeature <- mean(x$nFeature_RNA) - 2*sd(x$nFeature_RNA);
  mean_MT <- mean(x$percent.mt)
  upper_MT <- mean(x$percent.mt) + 2*sd(x$percent.mt);
  
  x$mean_nCount <- mean_nCount
  x$upper_nCount <- upper_nCount;
  x$lower_nCount <- lower_nCount;
  x$mean_nFeature <- mean_nFeature
  x$upper_nFeature <- upper_nFeature;
  x$lower_nFeature <- lower_nFeature;
  x$mean_MT <- mean_MT;
  x$upper_MT <- upper_MT;
  
  subset(x, subset = nCount_RNA < upper_nCount & 
           nCount_RNA > lower_nCount & 
           nFeature_RNA < upper_nFeature & 
           nFeature_RNA > lower_nFeature & 
           percent.mt < 10);
}

# Filter cells by established quality metrics
tic() 
Eze_dataset_list <- lapply(Eze_dataset_list, quality_metric_subset_function)
toc()
# Merge all cortex datasets
tic() # 684.424 sec
human_cortex_seurat_merged <- merge(Eze_dataset_list[[1]], y = Eze_dataset_list[-1])
toc()
# Normalize merged Seurat object
tic() # 231.397 sec
human_cortex_seurat <- NormalizeData(human_cortex_seurat_merged, normalization.method = "LogNormalize")
toc()
# Define highly-variable feature set
tic() # 358.366 sec
human_cortex_seurat <- FindVariableFeatures(human_cortex_seurat)
toc()

# Integrate human cortex datasets by their sample identities and define clusters
human_cortex_seurat <- RunFastMNN(object.list = SplitObject(human_cortex_seurat, split.by = "orig.ident"), k = 20, d = 50)
human_cortex_seurat <- FindNeighbors(human_cortex_seurat, reduction = "mnn", dims = 1:30)
human_cortex_seurat <- FindClusters(human_cortex_seurat, resolution = 0.6)

# Identify cell type markers for each cluster in integrated dataset
human_cortex_seurat.markers <- FindAllMarkers(human_cortex_seurat, min.pct = 0.25, only.pos=TRUE)
write.table(human_cortex_seurat.markers, file = "./Eze_et_al_2021_analysis/human_cortex_markers")
human_cortex_seurat.markers <- read.table("./Eze_et_al_2021_analysis/human_cortex_markers")

# Label clusters with cell type identification based on cell type markers
human_cortex_cell_types <- list('0' = 'IMCN', '1' = 'CIN', '2' = 'ULCN1', '3' = 'ULCN2', '4' = 'AS+ORG', '5' = 'ULCN3',
                                '6' = 'DLCN', '7' = 'IPC', '8' = 'IPC (G2/M)', '9' = 'RG (G1/S)', '10' = 'MES',
                                '11' = 'OL', '12' = 'RG (G2/M)', '13' = 'ORG', '14' = 'MG', '15' = 'AS', '16' = 'PGG')
human_cortex_seurat$old.ident <- Idents(human_cortex_seurat)
human_cortex_seurat <- RenameIdents(object = human_cortex_seurat, human_cortex_cell_types)
Idents(human_cortex_seurat) <- fct_relevel(Idents(human_cortex_seurat), sort)

# Subcluster AS+ORG cluster and rename IPC identity to add cell cycle information
human_cortex_seurat <- RenameIdents(human_cortex_seurat, list("ORG" = "INPG", "IPC" = "IPC (G1/S)"))
human_cortex_seurat <- FindSubCluster(human_cortex_seurat, cluster = "AS+ORG", graph.name = "RNA_snn",
                                      subcluster.name = "AS.ORG_sc")
Idents(human_cortex_seurat) <- human_cortex_seurat$AS.ORG_sc
human_cortex_seurat <- RenameIdents(human_cortex_seurat, list("AS+ORG_1" = "ORG (G1/S)", "AS+ORG_4" = "ORG (G1/S)",
                                                              "AS+ORG_5" = "ORG (G2/M)", "AS+ORG_0" = "AS",
                                                              "AS+ORG_2" = "AS", "AS+ORG_3" = "AS",
                                                              "AS+ORG_6" = "AS"))
human_cortex_progenitor_seurat <- subset(human_cortex_seurat, idents = c("RG (G2/M)", "IPC (G2/M)", "ORG (G1/S)",
                                                                         "RG (G1/S)", "IPC (G1/S)"))
human_cortex_progenitor_seurat <- ScaleData(human_cortex_progenitor_seurat)

# Creating matrices of average expression from proliferation assay hits
Geller_gene_hits <- read.table(file = "./../Evan_revisions/Geller_gene_hits.csv", sep = ",", header = TRUE)
rownames(Geller_gene_hits) <- Geller_gene_hits$Gene.Symbol

human_cortex_progenitor_avg_exp <- AverageExpression(human_cortex_progenitor_seurat,
                                                     features = Geller_gene_hits$Gene.Symbol, slot = "scale.data")
human_cortex_progenitor_avg_exp <- human_cortex_progenitor_avg_exp$RNA
human_cortex_progenitor_avg_exp <- human_cortex_progenitor_avg_exp[,c("IPC (G1/S)", "IPC (G2/M)", "RG (G1/S)", "RG (G2/M)", "ORG (G1/S)")]
num_centers = 10
human_cortex_progenitor_kmeans <- kmeans(human_cortex_progenitor_avg_exp, centers = num_centers)

avg_exp_by_cluster <- lapply(1:num_centers,
                             function(x) human_cortex_progenitor_avg_exp[human_cortex_progenitor_kmeans$cluster == x,])
avg_exp_by_cell_type_ann <- lapply(avg_exp_by_cluster, function(x) colMeans(x))
t4_pheno_by_gene_ann <- lapply(avg_exp_by_cluster, function(x) {Geller_gene_hits[rownames(x),]$t4.phenotype})
t8_pheno_by_gene_ann <- lapply(avg_exp_by_cluster, function(x) {Geller_gene_hits[rownames(x),]$t8.phenotype})
t12_pheno_by_gene_ann <- lapply(avg_exp_by_cluster, function(x) {Geller_gene_hits[rownames(x),]$t12.phenotype})

# Create heatmap annotations for cell cycle, cell type, average expression
cell_type_colors <- c("#F79256", "#FBD1A2", "#7DCFB6", "#00B2CA", "#1D4E89")
names(cell_type_colors) <- c("IPC (G1/S)", "IPC (G2/M)", "RG (G1/S)", "RG (G2/M)", "ORG (G1/S)")
cell_cycle_colors <- c("light grey", "#386C0B", "light grey", "#386C0B", "light grey")
names(cell_cycle_colors) <- c("G1/S", "G2/M", "G1/S", "G2/M", "G1/S")
col_ann_vector <- lapply(avg_exp_by_cell_type_ann, function(x) {HeatmapAnnotation(avg_exp = x,
                                                                                  cell_cycle = names(cell_cycle_colors),
                                                                                  cell_type = names(cell_type_colors),
                                                                                  col = list("avg_exp" = colorRamp2(c(min(x), max(x)), c("white", "orange")),
                                                                                             "cell_type" = cell_type_colors,
                                                                                             "cell_cycle" = cell_cycle_colors),
                                                                                  annotation_label = c("Avg Exp", "Cell Cycle", "Cell Type"),
                                                                                  annotation_name_gp = gpar(fontsize = 8, fontfamily = "Helvetica"),
                                                                                  annotation_name_rot = 30,
                                                                                  show_legend = FALSE)})

# Create heatmap annotation for proliferation phenotype
row_ann_vector <- lapply(1:num_centers, function(x) {rowAnnotation(t4_pheno = t4_pheno_by_gene_ann[[x]],
                                                          t8_pheno = t8_pheno_by_gene_ann[[x]],
                                                          t12_pheno = t12_pheno_by_gene_ann[[x]],
                                                          col = list(t4_pheno = c("negative" = "#BA110C",
                                                                     "neutral" = "light grey",
                                                                     "positive" = "light green"),
                                                                     t8_pheno = c("negative" = "#BA110C",
                                                                                  "neutral" = "light grey",
                                                                                  "positive" = "light green"),
                                                                     t12_pheno = c("negative" = "#BA110C",
                                                                                  "neutral" = "light grey",
                                                                                  "positive" = "light green")),
                                                          annotation_label = c("T4 Pheno", "T8 Pheno", "T12 Pheno"),
                                                          annotation_name_gp = gpar(fontsize = 8, fontfamily = "Helvetica"),
                                                          show_legend = FALSE)})

# Create heatmap for each cluster
ann_heatmaps_by_cluster <- lapply(1:10, function(x) {
                                  draw(Heatmap(avg_exp_by_cluster[[x]],
                                          name = glue::glue("Cluster { x }"),
                                          top_annotation = col_ann_vector[[x]],
                                          left_annotation = row_ann_vector[[x]],
                                          row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
                                          column_names_gp = gpar(fontfamily = "Helvetica"),
                                          show_column_dend = FALSE,
                                          show_row_dend = FALSE,
                                          row_names_rot = 30,
                                          cluster_columns = TRUE,
                                          show_heatmap_legend = FALSE), padding = unit(c(2, 2, 5, 2), "mm"))})
names(ann_heatmaps_by_cluster) <- c("IPC-biased", "ORG-biased", "Low Enrichment", "G2/M-biased", "Ribosomal Genes",
                                    "G2/M-biased", "ORG-biased", "ORG-biased", "IPC-biased", "RG(G1/S)-biased")
names(avg_exp_by_cluster) <- names(ann_heatmaps_by_cluster)

saveRDS(avg_exp_by_cluster, file = "~/Dropbox/Noonan Lab/workspace/Evan_revisions/avg_exp_by_cluster.rds")
avg_exp_by_cluster_Evan <- readRDS(file = "~/Dropbox/Noonan Lab/workspace/Evan_revisions/avg_exp_by_cluster.rds")

saveRDS(ann_heatmaps_by_cluster, file = "~/Dropbox/Noonan Lab/workspace/Evan_revisions/annotated_heatmaps_by_cluster.rds")
ann_heatmaps_by_cluster <- readRDS(file = "~/Dropbox/Noonan Lab/workspace/Evan_revisions/annotated_heatmaps_by_cluster.rds")

# Create supplementary table of cluster assignments for genes

avg_exp_by_cluster <- lapply(avg_exp_by_cluster, function(x) {as.data.frame(x) %>% mutate(gene = rownames(x))})
avg_exp_by_cluster <- lapply(1:10, function(x) {avg_exp_by_cluster[[x]] %>% 
    mutate(cluster = rep.int(x, dim(avg_exp_by_cluster[[x]])[1]))})
cluster_gene_assignment <- do.call("rbind", avg_exp_by_cluster) %>% select(gene, cluster)
cluster_name_list <- names(avg_exp_by_cluster)
names(cluster_name_list) <- 1:10
cluster_gene_assignment$cluster_name <- cluster_name_list[cluster_gene_assignment$cluster]

write.table(cluster_gene_assignment, 
            file = "~/Dropbox/Noonan Lab/workspace/Evan_revisions/cluster_assignment_supp_table.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE)
cluster_gene_assignment <- read.table(file = "~/Dropbox/Noonan Lab/workspace/Evan_revisions/cluster_assignment_supp_table.tsv",
                                      sep = "\t", header = TRUE)

# Density plots to show gene expression on UMAP

png(filename = "./raw_plots/Geller_density_POU3F3.png", width = 750, height = 650)
plot_density(human_cortex_progenitor_seurat, features = c("POU3F3")) + 
  scale_color_gradient(low = "#D3D3D3", high = "#CC5801") + 
  scale_x_continuous(limits = c(1, 13)) + 
  scale_y_continuous(limits = c(-3, 7)) + 
  th_cell_type_dimplot_theme + 
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
dev.off()

# UMAP representation of just fetal human cortex progenitors

human_cortex_cell_type_barcodes <- lapply(c("IPC (G1/S)", "IPC (G2/M)",
                                            "RG (G1/S)", "RG (G2/M)",
                                            "ORG (G1/S)"),
                                          function(x) {
                                            names(Idents(human_cortex_progenitor_seurat))[Idents(human_cortex_progenitor_seurat) == x]
                                          })
png(filename = "./raw_plots/human_cortex_progenitor_UMAP.png", width = 750, height = 650)
DimPlot(human_cortex_progenitor_seurat, pt.size = 1, raster = FALSE, order = rev(c("IPC (G1/S)", "IPC (G2/M)", "RG (G1/S)", "RG (G2/M)", "ORG (G1/S)"))) + 
  scale_x_continuous(limits = c(1, 13)) + 
  scale_y_continuous(limits = c(-3, 7)) + 
  scale_color_manual(values = list("IPC (G1/S)" = "#F79256", "IPC (G2/M)" = "#FBD1A2", "RG (G1/S)" = "#7DCFB6", "RG (G2/M)" = "#00B2CA", "ORG (G1/S)" = "#1D4E89")) +
  th_cell_type_dimplot_theme + 
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
dev.off()

# UMAP representation of fetal human cortex with highlighted progenitors

DimPlot(human_cortex_seurat, cells.highlight = human_cortex_cell_type_barcodes,
        cols.highlight = c("#F79256", "#FBD1A2", "#7DCFB6", "#00B2CA", "#1D4E89"), pt.size = 1, raster = FALSE) + 
  th_cell_type_dimplot_theme +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

# Create background gene set for GO analysis from supp table 2

Geller_supp2 <- read.csv("./Table_S2.csv")
Geller_gene_targets <- Geller_supp2$name[Geller_supp2$name %in% rownames(human_cortex_progenitor_seurat)]


