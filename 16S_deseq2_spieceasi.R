### Importing Qiime2 Output and Creating Phyloseq Object ###
devtools::install_github("jbisanz/qiime2R")
install.packages("tidyverse")
install.packages("stringr")

library(phyloseq)
library(qiime2R)
library(dplyr)
library(tidyverse)
library(stringr)
library(DESeq2)
library(ggplot2)
library(SpiecEasi)
library(igraph)
library(Matrix)

### Phyloseq object creation ###

otu_table <- read_qza("trial_table.qza")$data
otu_table <- otu_table[, c('sample-3', 'sample-4', 'sample-5', 'sample-6', 'sample-7',
                           'sample-8', 'sample-9', 'sample-10', 'sample-11', 'sample-12')]
# OTU table count filtering
row_sums <- rowSums(otu_table, na.rm = TRUE)
# Identify rows where the sum is >= 100
rows_to_keep <- row_sums >= 100

# Filter rows based on the identified rows to keep
otu_table <- otu_table[rows_to_keep, ]

taxonomy <- read_qza("trial_taxonomy.qza")$data

# Function to create columns for taxonomy ranks.
transform_taxa_df <- function(taxa_vector) {
  # Initialize an empty dataframe
  taxa_df <- data.frame(phylum = character(), class = character(), order = character(), 
                        family = character(), genus = character(), species = character(), stringsAsFactors = FALSE)
  
  # Iterate over the taxa_vector
  for (taxa_str in taxa_vector) {
    taxa_parts <- str_split(taxa_str, pattern = ";", simplify = FALSE)[[1]]
    
    # Extract parts using regex and return NA if not found
    phylum_part <- ifelse(any(grepl("p__", taxa_parts)), taxa_parts[grep("p__", taxa_parts)], NA)
    class_part <- ifelse(any(grepl("c__", taxa_parts)), taxa_parts[grep("c__", taxa_parts)], NA)
    order_part <- ifelse(any(grepl("o__", taxa_parts)), taxa_parts[grep("o__", taxa_parts)], NA)
    family_part <- ifelse(any(grepl("f__", taxa_parts)), taxa_parts[grep("f__", taxa_parts)], NA)
    genus_part <- ifelse(any(grepl("g__", taxa_parts)), taxa_parts[grep("g__", taxa_parts)], NA)
    species_part <- ifelse(any(grepl("s__", taxa_parts)), taxa_parts[grep("s__", taxa_parts)], NA)
    
    # Combine into a dataframe
    taxa_row <- data.frame(phylum = phylum_part, class = class_part, order = order_part, 
                           family = family_part, genus = genus_part, species = species_part, stringsAsFactors = FALSE)
    
    # Bind the row to the main dataframe
    taxa_df <- rbind(taxa_df, taxa_row)
  }
  
  return(taxa_df)
}

# Apply the function to the `taxa` column and bind the results to the original dataframe
taxa_df <- transform_taxa_df(taxonomy$Taxon)
taxonomy <- bind_cols(taxonomy, taxa_df)
rownames(taxonomy) <- taxonomy$Feature.ID
taxonomy <- taxonomy %>%
  mutate %>%
  select(-Feature.ID, -Taxon)



metadata <- read.csv("manifest_wtko_sd.txt", sep="\t")
rownames(metadata) <- metadata$sample.id
metadata <- metadata %>%
  mutate %>%
  select(-sample.id)
tree <- read_qza("rooted-tree.qza")$data

# Assuming you have already loaded or created the otu_table, taxonomy, and tree objects

# Check OTU IDs in the OTU table and taxonomy table
otu_ids <- rownames(otu_table)
taxa_ids <- rownames(taxonomy)

# Find mismatches
mismatched_ids <- setdiff(otu_ids, taxa_ids)
mismatched_ids_reverse <- setdiff(taxa_ids, otu_ids)

# Print mismatched IDs
print(mismatched_ids)
print(mismatched_ids_reverse)

# If using a phylogenetic tree, ensure its tip labels match OTU IDs as well
tree_tip_labels <- tree$tip.label
mismatched_ids_tree <- setdiff(otu_ids, tree_tip_labels)
mismatched_ids_tree_reverse <- setdiff(tree_tip_labels, otu_ids)

# Print mismatches with the tree
print(mismatched_ids_tree)
print(mismatched_ids_tree_reverse)


ps <- phyloseq(otu_table(otu_table, taxa_are_rows = TRUE),
               tax_table(as.matrix(taxonomy)),
               sample_data(metadata),
               phy_tree(tree))

### DESeq2 differential abundance calculation ###

dds <- phyloseq_to_deseq2(ps, ~ Genotype) # Create deseq2 object.

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds, fitType="local")

resultsNames(dds)

res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]


alpha <- 0.01
sigtab <- res[(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab <- sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab <- posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "phylum", "class", "family", "genus", "species")]
head(posigtab)

negsigtab <- sigtab[sigtab[, "log2FoldChange"] < 0, ]
negsigtab <- negsigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "phylum", "class", "family", "genus", "species")]
head(negsigtab)


theme_set(theme_bw())
sigtabgen <- subset(sigtab, !is.na(species))
# Phylum order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$species, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$species = factor(as.character(sigtabgen$species), levels=names(x))
ggplot(sigtabgen, aes(y=species, x=log2FoldChange, color=species)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

### Spieceasi network construction ###

se.wtkosd <- spiec.easi(ps, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=200))

mb <- adj2igraph(getRefit(se.wtkosd),  vertex.attr=list(name=taxa_names(ps)))
plot_network(mb, ps, type='taxa', color="genus")
