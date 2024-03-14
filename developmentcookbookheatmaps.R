BiocManager::install("DESeq2")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("tidyr")
library(ggplot2)
library(DESeq2)
library(dplyr)
library(tidyr)

### Read the data and calculate size factors with the help of deseq2 ###

df <- read.csv("Olig_p21_cko_wt.csv")
rownames(df) <- df$gene_id
df <- df[, -which(names(df) == "gene_id")]

metadata <- read.table("metadata_OL_p21_cko_wt.txt", header=TRUE)
rownames(metadata) <- metadata$Sample
metadata <- dplyr::select(metadata, -Sample)
  
# Read the data as a deseq2 object
dds <- DESeqDataSetFromMatrix(countData = df,
                                colData = metadata,
                                design= ~ Condition)

dds <- DESeq(dds) # Create the object.

res <- results(dds) # Important for retrieving p-values later on.

# Retrieve the size factor calculations from DESeq2
sizefactors <- dds$sizeFactor

### Heatmap as in Düking et al. ###
matrixgenenames <- read.csv('Olig_p21_cko_wt_genenames.csv')
rownames(matrixgenenames) <- matrixgenenames$gene_id
matrixgenenames <- matrixgenenames %>%
  mutate %>%
  dplyr::select(-gene_id)


# Normalize your gene counts to size factors.
columns_to_normalize <- c("SI.cert.ctrl.oligo.p21.1", "SI.cert.ctrl.oligo.p21.2",  
                          "SI.cert.ctrl.oligo.p21.3", "SI.cert.ctrl.oligo.p21.4",
                          "SI.cert.ctrl.oligo.p21.5", "SI.cert.ctrl.oligo.p21.6",
                          "SI.hmgcs2.cKO.oligos.p21.1", "SI.hmgcs2.cKO.oligos.p21.2",
                          "SI.hmgcs2.cKO.oligos.p21.3", "SI.hmgcs2.cKO.oligos.p21.5")

sizefactored_matrixgenenames <- sweep(matrixgenenames[, columns_to_normalize, drop = FALSE], 2, sizefactors, "/")
sizefactored_matrixgenenames <- cbind(matrixgenenames["gene_name"], sizefactored_matrixgenenames)

wt_samples <- c("SI.cert.ctrl.oligo.p21.1", "SI.cert.ctrl.oligo.p21.2",  
                "SI.cert.ctrl.oligo.p21.3", "SI.cert.ctrl.oligo.p21.4",
                "SI.cert.ctrl.oligo.p21.5", "SI.cert.ctrl.oligo.p21.6")

ko_samples <- c("SI.hmgcs2.cKO.oligos.p21.1", "SI.hmgcs2.cKO.oligos.p21.2",
                "SI.hmgcs2.cKO.oligos.p21.3", "SI.hmgcs2.cKO.oligos.p21.5")

# Calculate the mean values for WT.
wt_means <- rowMeans(sizefactored_matrixgenenames[, wt_samples]) 

# Avoid infinity
epsilon <- 0.0000001
wt_means_epsilon <- wt_means + epsilon

# Calculate the mean values for KO.
ko_means <- rowMeans(sizefactored_matrixgenenames[, ko_samples])

# Normalize KO values by WT means
normalized_ko <- sweep(sizefactored_matrixgenenames[, ko_samples], 1, wt_means_epsilon, FUN="/")

# Avoid infinity
normalized_ko_epsilon <- normalized_ko + epsilon

# Calculate log2 fold change for normalized KO values
log2_fold_change_ko <- log2(normalized_ko_epsilon)

gene_names <- sizefactored_matrixgenenames$gene_name

results_df <- data.frame(Gene = gene_names,
                         WT_Mean = wt_means,
                         KO_Mean = ko_means,
                         log2_fold_change_ko,
                         p_adj = res$padj)

colnames(results_df) <- c("Gene", "WT_Mean", "KO_Mean",                   
                          "log2(cKOp21.1/WT)", "log2(cKOp21.2/WT)",
                          "log2(cKOp21.3/WT)", "log2(cKOp21.4/WT)", "p_adj")

# Glycolysis

glycolysis_genes <- c('Slc2a1', 'Hk1', 'Hk2', 'Gpi1', 'Pfkm', 'Pfkp', 'Aldoa',
                      'Aldob', 'Tpi1', 'Gapdh', 'Pgk1', 'Pgam1', 'Eno1', 'Eno3',
                      'Pkm', 'Ldha', 'Ldhb')

glycolysis <- results_df[results_df$Gene %in% glycolysis_genes,]

# Remove low count Aldoa duplicate and Aldob.
glycolysis <- glycolysis[-which(rownames(glycolysis) %in% c("ENSMUSG00000114515", "ENSMUSG00000028307")), ]

# Convert the glycolysis dataframe into long format for ggplot2.
glycolysis_long <- glycolysis %>%
  pivot_longer(cols = starts_with("log2"), 
               names_to = "Condition", 
               values_to = "Log2FoldChange") %>%
  mutate(Condition = gsub("log2\\(|/WT\\)", "", Condition)) # Clean up the Condition names

# Plot the heatmap.
ggplot(glycolysis_long, aes(x = Condition, y = Gene, fill = Log2FoldChange)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-2, 2), 
                       oob = scales::squish) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Log2 Fold Change")
