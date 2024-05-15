### Workflow for more in-depth analysis of bulk RNA-seq analysis ###
### Emphasis on pathways ###
### See developmentcookbookheatmaps.R for heatmap emphasis ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("vsn")
BiocManager::install('EnhancedVolcano')
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("org.Mm.eg.db")
install.packages("pathfindR")
library(DESeq2)
library(dplyr)
library(apeglm)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(pathfindR)
library(clusterProfiler)
library(pathview)
library(org.Mm.eg.db)

### Insert the cell types of interest from the raw data ###
### Edit the Raw data and copy the columns of ENSEMBL and cell types of interest to a new excel sheet for df ###
### Save the excel file as csv with separator field of comma (,) ###
### Create one more csv file as described above but include gene names as a column as well. This will become useful for Pathview ###

### Define data and metadata ###
df <- read.csv("EC_p9_p21_Ctrl.csv")
rownames(df) <- df$gene_id
df <- df[, -which(names(df) == "gene_id")]

metadata <- read.table("metadata_EC_p9_p21_wt.txt", header=TRUE)
rownames(metadata) <- metadata$Sample
metadata <- dplyr::select(metadata, -Sample)

### Start DeSeq Workflow ###
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = metadata,
                              design= ~ Condition)

dds <- DESeq(dds) # Create the object.

res <- results(dds) # Retrieve the results.
res

resultsNames(dds) # Retrieve information on the comparison condition string.
mcols(res)$description # Retrieve information on the column descriptions of res object.

resLFC <- lfcShrink(dds, coef="Condition_EC_Control_p9_vs_EC_Control_p21", type="apeglm") # Calculate the log fold
                                                                                          # changes in shrunk dataset.
resLFC

resOrdered <- res[order(res$padj),] # Order res object for increasing p-values.

summary(res)
sum(res$padj < 0.05, na.rm=TRUE) # How many adjusted p-values were less than 0.05?

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, "ENSMUSG00000027875", intgroup="Condition") # Plot for Hmgcs2 absolute counts.

res_padj005 <- res[which(res$"padj" < 0.05), ] # Filter for padj lower than 0.05.

# Filter the res_padj01 object for LFC of lower than -0.5 or higher than 0.5.
res_filtered_adjp_LFC <- res_padj005[which(res_padj005$log2FoldChange < -0.5 | res_padj005$log2FoldChange > 0.5), ]
res_p9_up <- res_filtered_adjp_LFC[which(res_filtered_adjp_LFC$log2FoldChange > 0) ,] # Filter for genes upregulated in WT.
res_p21_up <- res_filtered_adjp_LFC[which(res_filtered_adjp_LFC$log2FoldChange < 0) ,] # Filter for genes upregulated in cKO.
write.table(res_p9_up, "res005_p9up.txt") # Write the WT upregulated genes.
write.table(res_p21_up, "res005_p21up.txt") # Write the cKO upregulated genes.

### Data transformations and visualization ###

# Observe that we use absolute counts for the visualizations.
vsd <- vst(dds, blind=FALSE) # Variance stabilizing transformation of the dataset with non-blind dispersion.
rld <- rlog(dds, blind=FALSE) # Regularized logarithm transformation of the dataset with non-blind dispersion.
head(assay(vsd), 3)
head(assay(rld), 3)

ntd <- normTransform(dds) # Simple log2(n + 1).

# Assessing the effect of transformation on the data variance.
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

### Data quality assessment by sample clustering and visualization ###

# Heatmap clustering of the sample count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
heatdf <- as.data.frame(colData(dds)[,"Condition"])
colnames(heatdf) <- c("Condition")
rownames(heatdf) <- rownames(metadata)

pheatmap(
  assay(rld)[select, ],
  cluster_rows = TRUE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  annotation_col = heatdf,
  fill = colorRampPalette(c("blue", "white", "red"))(100)
)  

# Heatmap clustering of the sample-to-samples distances
sampleDistsntd <- dist(t(assay(ntd)))
sampleDistsvsd <- dist(t(assay(vsd)))
sampleDistsrld <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDistsrld)
rownames(sampleDistMatrix) <- paste(rld$Condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDistsrld,
         clustering_distance_cols=sampleDistsrld,
         col=colors)

# PCA of the samples
plotPCA(rld, intgroup="Condition")

EnhancedVolcano(res,
                lab = rownames(res),
                pCutoff = 10e-25,
                FCcutoff = 0.5,
                #drawConnectors = TRUE,
                x = 'log2FoldChange',
                y = 'padj')

### PathfindR for investigating protein interaction networks ###
# Create the format of the required input for PathfindR.
gene_names <- read.csv('EC_p9_p21_Ctrl_genenames.csv')
gene_names <- gene_names$gene_name
gene_names

res_index <- data.frame(Gene_symbol = rownames(res))
data <- res[, c('log2FoldChange', 'padj')]

data <- cbind(data, gene_names)
data

names(data)[1] <- "logFC"  # Renaming the second column to "New_X"
names(data)[2] <- "FDR_p"
names(data)[3] <- "Gene_symbol"
data <- data[, c("Gene_symbol", "logFC", "FDR_p")]
data <- as.data.frame(data)
data <- data[complete.cases(data$FDR_p), ] # Remove NAs.

data005 <- data[which(data$"FDR_p" < 0.05), ] # Filter for padj lower than 0.05.

# Filter the data005 dataframe for LFC of lower than -0.5 or higher than 0.5.
data_filtered_adjp_LFC <- data005[which(data005$logFC < -0.5 | data005$logFC > 0.5), ]

mmu_output <- run_pathfindR(input = data_filtered_adjp_LFC,
                            p_val_threshold = 0.05,
                            convert2alias = FALSE,
                            gene_sets = "mmu_KEGG",
                            pin_name_path = "mmu_STRING")


# display the heatmap of hierarchical clustering
clustered_mmuoutput <- cluster_enriched_terms(mmu_output, plot_hmap=TRUE)
visualize_terms(clustered_mmuoutput, hsa_KEGG=FALSE) # Optional if you want to generate PNG files with clusters.
clustered_mmuoutput_fuzzy <- cluster_enriched_terms(mmu_output, method="fuzzy")

### Pathway mapping from the fold changes with PathView ###

# Assuming your genes are in SYMBOL format
gene_symbols <- data_filtered_adjp_LFC$Gene_symbol

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids_filtered <- na.omit(entrez_ids)

# Perform KEGG enrichment analysis using Entrez IDs
kegg_result <- enrichKEGG(gene = entrez_ids_filtered, organism = 'mmu', pvalueCutoff = 0.05)
kegg_result <- as.data.frame(kegg_result)
# Check the results
kegg_result

fold_changes_aligned_with_entrez <- setNames(data_filtered_adjp_LFC$logFC, entrez_ids) # Align the foldchanges and gene ids.

datakegg <- data.frame(gene=names(entrez_ids), foldChange=data_filtered_adjp_LFC$logFC) # Create a new df for pathview

# Create heatmaps for your maps of interest from KEGG pathways.
pathways = c("mmu00010", "mmu00020", "mmu00061", "mmu00062", "mmu00071", "mmu00100", "mmu00120",
             "mmu00140", "mmu00190", "mmu00230", "mmu00561", "mmu00564", "mmu00565", "mmu00590", 
             "mmu00600", "mmu00650", "mmu02010", "mmu03320", "mmu04020", "mmu04066", "mmu04146",
             "mmu04216", "mmu04370", "mmu04510", "mmu04514", "mmu04530", "mmu04540", "mmu04550",
             "mmu04930", "mmu05010", "mmu05022", "mmu05200")
pathview(gene.data=fold_changes_aligned_with_entrez, pathway.id=pathways, 
         specie ="mmu", gene.idtype="ENTREZ", limit=list(gene=2), low=list(gene = "blue", cpd = "blue"))
