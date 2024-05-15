### Repository contains workflow for network analysis of 16S rRNA sequencing data ###
### Repository assumes *relative abundance*, not absolute abundances of taxonomical ranks ###

install.packages("devtools")
install_github("zdk123/SpiecEasi")
install_github("hallucigenia-sparsa/seqgroup", build_vignettes = TRUE)
BiocManager::install("RCy3")
install_github("frankkramer-lab/ndexr")

library(devtools)
library(SpiecEasi)
library(seqgroup)
library(phyloseq)
library(RCy3)
library(dplyr)
library(ndexr)
library(igraph)
library(stringr)

df <- read.csv("~/Saher/16S_analysis/conet/16s_relative_abundance_transposed.csv")

taxa_info <- rownames(df)

transform_taxa <- function(taxa_str) {
  # Split the string into components based on the semicolon separator
  taxa_parts <- str_split(taxa_str, pattern = ";", simplify = TRUE)
  
  # Extract the genus and species parts
  genus_part <- taxa_parts[grep("g__", taxa_parts)]
  species_part <- taxa_parts[grep("s__", taxa_parts)]
  
  # Combine the genus and species with a semicolon if both are present, or return an NA if either is missing
  if (length(genus_part) > 0 || length(species_part) > 0) {
    return(paste(genus_part, species_part, sep = ";"))
  } else {
    return(NA)
  }
}

df$transformed_taxa <- sapply(taxa_info, transform_taxa)
df_clean <- df[!is.na(df$transformed_taxa), ]
unique_transformed_taxa <- df_clean$transformed_taxa[duplicated(df_clean$transformed_taxa) | duplicated(df_clean$transformed_taxa, fromLast = TRUE)]
df_clean <- df_clean[!df_clean$transformed_taxa %in% unique_transformed_taxa, ]
rownames(df_clean) <- df_clean$transformed_taxa
df_clean <- df_clean %>%
  mutate %>%
  select(-transformed_taxa)


taxalineages <- read.table("taxalineages.csv", sep=',', header=1)
rownames(taxalineages) <- taxalineages$Index

mcyc_pathway_abundances <- read.delim("~/Saher/16S_analysis/conet/mcyc_pathway_abundances.tsv", comment.char="#")
metadata <- read.delim("~/Saher/16S_analysis/conet/metadata.txt")
metadata_wtko_sd <- read.delim("~/Saher/16S_analysis/conet/metadata_wtko_sd.txt")
rownames(metadata) <- as.character(metadata$Samples)
metadata <- metadata %>%
  mutate %>%
  select(-Samples)
rownames(metadata_wtko_sd) <- as.character(metadata_wtko_sd$Samples)
metadata_wtko_sd <- metadata_wtko_sd %>%
  mutate %>%
  select(-Samples)
groups=as.vector(metadata$Group)
groups_wtko_sd =as.vector(metadata_wtko_sd$Group)

compareGroups(df_clean,groups=groups_wtko_sd,property="alpha",pvalViz = TRUE)
compareGroups(df_clean,groups=groups_wtko_sd,property="beta",pvalViz = TRUE)
compareGroups(mcyc_pathway_abundances,groups=groups_wtko_sd,property="beta",pvalViz = TRUE)
compareGroups(df,groups=groups_wtko_sd,property="beta",method="DM")

controlsd.indices=(c(2, 3, 4, 22, 23))
kosd.indices=(c(16, 17, 19, 20, 21))
par(mfrow=c(1,1))
par(mar=c(4, 4, 2, 2))
taxon.color.map=groupBarplot(df[,control.indices],topTaxa = 20, randSampleNum = 20, extendTaxonColorMap = TRUE, main="WT_SD")
taxon.color.map=groupBarplot(df[,ko.indices],topTaxa = 20, randSampleNum = 20, taxon.color.map = taxon.color.map, extendTaxonColorMap = TRUE, main="KO_SD")

groupcolor <- c("WT_SD", "KO_SD", "WT_KD", "KO_KD")
# Define the new contrasting colors
colors <- c('#377eb8', '#ff7f00', '#e41a1c', '#984ea3')
# Create the map of predefined colors for the groups
groupColors <- setNames(colors, groupcolor)
# Display the color map
print(groupColors)

png("seqPCoA20.png", width = 1200, height = 1200, res = 100)
# Create your plot here
seqPCoA(df,groups=groups, topTaxa=20, groupColors=groupColors)
dev.off()  # Close the PNG device and save the plot

taxon="d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;__;__"
png("wtsd_kosd_lachnospiraceae.png", width = 1600, height = 900, res = 100)
compareDistribs(df,taxon=taxon,groups=groups_wtko_sd,group1="WT_SD",group2="KO_SD")
dev.off()

clusters=findClusters(df,k=NA, method="pam", qualityIndex = "CH")

png("seqPCoAClusters.png", width = 1200, height = 900, res = 100)
seqPCoA(df,groups=groups,clusters=clusters, topTaxa=30)
dev.off()

clus.table = table(groups, clusters)
chisq.test(clus.table)

# def.par=par(no.readonly = TRUE) # save previous par settings
layout(matrix(c(1,3,2,4),2,2,byrow = TRUE), c(2,3,2,3), TRUE) # define space for 4 figures
kosd.network=buildNetwork(df_clean[c('sample.3', 'sample.4', 'sample.5', 'sample.6', 'sample.7')], method=c("spearman","bray"), repNum=0, nameLevel="species") # run CoNet on IBD samples
wtsd.network=buildNetwork(df_clean[c('sample.8', 'sample.9', 'sample.10', 'sample.11', 'sample.12')], method=c("spearman","bray"), repNum=0, nameLevel="species") # run CoNet on control samples
kosd.aracne=buildNetwork(df_clean[c('sample.3', 'sample.4', 'sample.5', 'sample.6', 'sample.7')], method="aracne",repNum=100)
wtsd.aracne=buildNetwork(df_clean[c('sample.8', 'sample.9', 'sample.10', 'sample.11', 'sample.12')], method="aracne",repNum=100)

png("network_wtsd2.png", width = 2000, height = 1600, res = 100)
plot(wtsd.network,main="WT_SD Network")
dev.off()

png("network_kosd2.png", width = 2000, height = 1600, res = 100)
plot(kosd.network, main="KO_SD Network")
dev.off()

edges <- as_edgelist(kosd.aracne)

# Assuming no explicit 'interaction' type is available, just use a placeholder or specific relation if known
# For simplicity, this example uses "pp" (protein-protein interaction) as a placeholder

prev=10 # minimum occurrence in 10 samples
df_incidence=df_clean[c(controlsd.indices, kosd.indices)]
df_incidence[df_incidence>0]=1
rowSums=rowSums(df_incidence)
indices.prev=which(rowSums>=prev)
hmgcs2status <- groups_wtko_sd
hmgcs2status[groups_wtko_sd=='WT_SD'] = 0
hmgcs2status[groups_wtko_sd=='KO_SD'] = 1
df_extended=rbind(df_clean[c(controlsd.indices, kosd.indices)][indices.prev,],as.numeric(hmgcs2status))
rownames(df_extended)[33] <- "Hmgcs2"

wtko_sd_extended=buildNetwork(df_extended, method=c("spearman", "bray"), repNum=0,  nameLevel="species") # run CoNet on all samples
plot(wtko_sd_extended,main="Global network")

png("wtko_sd_hmgcs2.png", width = 2000, height = 1600, res = 100)
plot(wtko_sd_extended,main="WT vs KO Network")
dev.off()

createNetworkFromIgraph(wtko_sd_extended)

edges <- as_edgelist(wtko_sd_extended)


prev=20
mcyc_pathway_incidence=mcyc_pathway_abundances
mcyc_pathway_incidence[mcyc_pathway_incidence>0]=1
rowSums=rowSums(mcyc_pathway_incidence)
indices.functions.prev=which(rowSums>=prev)
bip.network=barebonesCoNet(abundances=df_clean[c(controlsd.indices, kosd.indices)][indices.prev,],metadata=t(mcyc_pathway_abundances[c(controlsd.indices, kosd.indices)][indices.functions.prev,]),min.occ=0,methods="pearson",T.up=0.7,T.down=-0.7)
plot(bip.network,main="Global species-function network")


prev=150 # minimum occurrence for functions
ibd_functions_incidence=ibd_functions
ibd_functions_incidence[ibd_functions_incidence>0]=1
rowSums=rowSums(ibd_functions_incidence)
indices.functions.prev=which(rowSums>=prev)
# CoNet expects metadata to have samples as rows, therefore functions have to be transposed
bip.network=barebonesCoNet(abundances=ibd_taxa[indices.prev,],metadata=t(ibd_functions[indices.functions.prev,]),min.occ=0,methods="pearson",T.up=0.7,T.down=-0.7)
