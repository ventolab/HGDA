#!/usr/bin/env R

###################
# Import packages #
###################

# Required packages
packages <- c("argparse", "biomaRt", "dplyr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

##################################
# Define command line parameters #
##################################

parser <- ArgumentParser()

parser$add_argument("--genes", type="character",
                    help = "Path to txt file with mouse ENSEMBL gene IDs you wish to convert to human.")

parser$add_argument("--outdir", type = "character", 
                    help = "Path to where you want to store the resulting converted human ENSEMBL IDs and gene names.")

args <- parser$parse_args()

mouse_genes <- args$genes
outdir <- args$outdir

print("Selected arguments: ")
print(paste0("Mouse genes list: ", mouse_genes))
print(paste0("Outdir: ", outdir))

################################################
# Read in gene list & convert gene ENSEMBL IDs #
################################################

# Connect to "genes" database on biomart and choose Mus Musculus dataset 
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")

mouse_ids <- scan(file = mouse_genes, what = "character", sep = "\n")

# Query for cross species comparison 
df_genes <- getBM(mart = ensembl,
      filters = c("ensembl_gene_id"),
      values = mouse_ids, 
      attributes = c("ensembl_gene_id", 
                     "external_gene_name", 
                     "hsapiens_homolog_ensembl_gene", 
                     "hsapiens_homolog_associated_gene_name", 
                     "hsapiens_homolog_orthology_type")
)

# Filter for one2one orthologs across species 
df_genes <- df_genes %>% filter(hsapiens_homolog_orthology_type == "ortholog_one2one") 

# Save datafram to .csv file 
write.csv(df_genes, paste0(outdir, "human_converted_ids.csv"), row.names = FALSE)

print(paste0('Mouse gene IDs have been successfully converted to human and saved :', outdir, 'human_converted_ids.csv'))
