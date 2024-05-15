#!/bin/bash
#SBATCH --job-name=sctenifoldyoungcolonhmgcs2
#SBATCH --account=all
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=resulelgin@gmail.com
#SBATCH --mem=150G
#SBATCH --cpus-per-task=20

# Load necessary modules
module load r  # Adjust the version as necessary

# Run your R script
Rscript sctenifoldknk.R
