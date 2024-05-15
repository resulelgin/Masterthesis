### Repository contains the analysis workflow of sctenifoldKNK output ###
### Outcome of sctenifoldKNK is expected to be obtained from sctenifoldknk.sh, calling sctenifoldknk.R, saving the output as .rds ###

endotheliaeae14hmgcs2ko <- readRDS("EAE14EndotheliaSC_hmgcs2ko_output.rds")
endothelianaivehmgcs2ko <- readRDS("EndotheliaNaiveSC_hmgcs2ko_output.rds")
microglianaiveoxct1ko <- readRDS("microglia_naive_oxct1ko_output.rds")
microgliaeae14oxct1ko <- readRDS("microglia_eae14_oxct1ko_output.rds")
youngcolon <- readRDS("youngcolon_hmgcs2ko.rds")

# Check the differentially regulated genes.
head(endotheliaeae14hmgcs2ko$diffRegulation, n = 50)
head(endothelianaivehmgcs2ko$diffRegulation, n = 50)
head(microglianaiveoxct1ko$diffRegulation, n= 50)
head(microgliaeae14oxct1ko$diffRegulation, n= 50)
head(youngcolon$diffRegulation, n=70)

ec_naive_hmgcs2_sig <- subset(endothelianaivehmgcs2ko$diffRegulation, endothelianaivehmgcs2ko$diffRegulation$p.adj < 0.05)
ec_eae14_hmgcs2_sig <- subset(endotheliaeae14hmgcs2ko$diffRegulation, endotheliaeae14hmgcs2ko$diffRegulation$p.adj < 0.05)
microglia_naive_oxct1_sig <- subset(microglianaiveoxct1ko$diffRegulation, microglianaiveoxct1ko$diffRegulation$p.adj < 0.05)
microglia_eae14_oxct1_sig <- subset(microgliaeae14oxct1ko$diffRegulation, microgliaeae14oxct1ko$diffRegulation$p.adj < 0.05)
youngcolon_hmgcs2_sig <- subset(youngcolon$diffRegulation, youngcolon$diffRegulation$p.adj < 0.10)

which(endothelianaivehmgcs2ko[["tensorNetworks"]][["WT"]]@Dimnames[[1]] == "Hmgcs2") # KO gene tensor position 1682
which(endotheliaeae14hmgcs2ko[["tensorNetworks"]][["WT"]]@Dimnames[[1]] == "Hmgcs2") # KO gene tensor position 1793
which(microglianaiveoxct1ko[["tensorNetworks"]][["WT"]]@Dimnames[[1]] == "Oxct1") # KO gene tensor position 7718
which(microgliaeae14oxct1ko[["tensorNetworks"]][["WT"]]@Dimnames[[1]] == "Oxct1") # KO gene tensor position 11669

# Manual Check #
kogene_position = 7718
gene_position = 7718
matrix_dimension = 9471
gene_correlation = (gene_position - 1) * matrix_dimension + kogene_position
microglianaiveoxct1ko[["tensorNetworks"]][["WT"]]@x[gene_correlation]

# Inference of the directionalities using the tensor correlations. Caution is advised, experimental.
correlation_retrieval <- function(matrix, ko, significance_list) {
  
  kogene_position <- which(matrix[["tensorNetworks"]][["WT"]]@Dimnames[[1]] == ko)
  matrix_dimension <- matrix[["tensorNetworks"]][["WT"]]@Dim[1]
  gene_positions <- as.numeric(rownames(significance_list))
  
  # Initialize an empty array to store results
  correlations_list <- c()
  
  # Loop over each gene in gene_positions
  for (gene_position in (gene_positions)) {
    
    # Calculate the resulting value for each gene
    resulting_value <- (gene_position - 1) * matrix_dimension + kogene_position
    
    # Retrieve the correlation values
    correlations <- matrix[["tensorNetworks"]][["WT"]]@x[resulting_value]
    
    # Append the correlations to the vector
    correlations_list <- append(correlations_list, correlations)
  }
  
  # Boolean logic to infer the directionality of the gene expression responses 
  significance_list$Correlations <- correlations_list
  significance_list$AdjustedLogfoldchange <- ifelse(
    significance_list$Correlations >= 0,
    -significance_list$FC,  # Apply minus if correlation is positive
    significance_list$FC   # Keep as is if correlation is negative
  )
  return (significance_list)
}
  

ec_naive_hmgcs2_adjusted <- correlation_retrieval(endothelianaivehmgcs2ko, "Hmgcs2", ec_naive_hmgcs2_sig)
ec_eae14_hmgcs2_adjusted <- correlation_retrieval(endotheliaeae14hmgcs2ko, "Hmgcs2", ec_eae14_hmgcs2_sig)
microglia_naive_oxct1_adjusted <- correlation_retrieval(microglianaiveoxct1ko, "Oxct1", microglia_naive_oxct1_sig)
microglia_eae14_oxct1_adjusted <- correlation_retrieval(microgliaeae14oxct1ko, "Oxct1", microglia_eae14_oxct1_sig)
youngcolon_adjusted <- correlation_retrieval(youngcolon, "Hmgcs2", youngcolon_hmgcs2_sig)
