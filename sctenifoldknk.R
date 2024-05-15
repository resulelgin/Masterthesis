### Workflow for efficient scRNA-seq in-silico gene knockout using scTenifoldKNK ###
### Workflow is expected to run with HPC ###
### Workflow optimally works with ~300 GB RAM, ~20 cores in less than 24 hours ###
### See sctenifoldknk.sh for sbatch submission ###

options(repos = c(CRAN = "https://cran.rstudio.com/"))
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")  # Set local library path
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)
# Install packages if not already installed
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  install_github('cailab-tamu/scTenifoldKnk')
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages('dplyr')
}

if (!requireNamespace("scTenifoldNet", quietly = TRUE)) {
  install.packages('scTenifoldNet')
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages('igraph')
}
if (!requireNamespace("rslurm", quietly = TRUE)) {
  install.packages('rslurm')
}

library(dplyr)
library(igraph)
library(scTenifoldKnk)
library(scTenifoldNet)
library(rslurm)

# Read data
X <- read.csv("/scratch/users/resul.elgin/Saher/sctenifoldknk/youngcolon/youngcolon.csv")
rownames(X) <- X$X
X <- X %>%
select(-X)

X_filtered <- X[rowSums(X) != 0, ]


X_hmgcs2ko <- scTenifoldKnk(X_filtered, gKO='Hmgcs2', qc = FALSE)
saveRDS(X_hmgcs2ko, file="youngcolon_hmgcs2ko.rds")
