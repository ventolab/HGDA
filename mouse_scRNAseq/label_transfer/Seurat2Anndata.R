#!/usr/bin/env R

###################
# Import packages #
###################

# Required packages
packages <- c("Seurat", "argparse")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# SeuratDisk (not currently available on CRAN)
if ("SeuratDisk" %in% rownames(installed.packages())) {
	library("SeuratDisk")
} else {
	if (!requireNamespace("remotes", quietly = TRUE)) {
  	install.packages("remotes")
	}
	remotes::install_github("mojaveazure/seurat-disk")
	library("SeuratDisk")
}

##################################
# Define command line parameters #
##################################

parser <- ArgumentParser()

parser$add_argument("--seurat_object", type="character",
                    help = "Path to Seurat object you wish to convert to Anndata")

parser$add_argument("--object_name",
                    type="character",
                    help = "Name of the saved object")

parser$add_argument("--outdir", type = "character", 
		    help = "Path to where you want to store the resulting Anndata object")

args <- parser$parse_args()

seurat_object <- args$seurat_object

object_name <- args$object_name

outdir <- args$outdir

print("Selected arguments: ")
print(paste0("Seurat Object: ", seurat_object))
print(paste0("Object name: ", object_name))
print(paste0("Outdir: ", outdir))

####################################
# Convert Seurat object to Anndata #
####################################

# Load Seurat object
SeuratObject <- readRDS(file = seurat_object)
print("Seurat object loaded") 

# Save intermediate format (h5Seurat)
SaveH5Seurat(SeuratObject, filename = paste0(outdir, object_name, ".h5Seurat"))
print("Seurat object saved as .h5seurat")

# Convert and save to h5ad file 
Convert(paste0(outdir, object_name, ".h5Seurat"), dest = "h5ad")

print("Seurat object was successfully converted to Anndata")
