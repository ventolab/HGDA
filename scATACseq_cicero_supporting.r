####### CICERO ON SUPPORTING CELLS ####### 

# Load libraries
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
library(SeuratDisk)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ensembldb)
library(JASPAR2018)
library(TFBSTools)
library(patchwork)
library(universalmotif)
library(MotifDb)
library(TFBSTools)
library(patchwork)
library(chromVAR)
library(motifmatchr)
library(cicero)

# Set directories
outdir <- "/nfs/team292/vl6/my_MULTIOME_dir/supporting_apr2021/"
experiment_prefix <- 'supporting_'

# Load Seurat object
ATAC_Seurat_withChromVar = readRDS(file = paste0(outdir, experiment_prefix, "_chromVar_binary.rds"))
DefaultAssay(ATAC_Seurat_withChromVar) <- "peaks"
ATAC_Seurat_withChromVar
Idents(ATAC_Seurat_withChromVar) <- ATAC_Seurat_withChromVar@meta.data$cell_type
Idents(ATAC_Seurat_withChromVar) <- factor(x = Idents(ATAC_Seurat_withChromVar),
                                    levels = c('coelEpi', 'sKITLG', 'sLGR5', 'sPAX8b', 'sPAX8m', 'preGC_I_OSR1',
                                               'ovarianSurf', 'preGC_II', 'preGC_II_hypoxia', 'preGC_III_Notch',
                                               'Sertoli','FetalLeydig-like'))

gonads_colors = c('#7b9e99','#91bd80','#bdb380','#d4db81','#70ccbe','#cc8fdb','#edb7b7','#d9abb7','#e08b8b', 
                  '#e64e74', '#aad3f2', '#60bddb')

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_supporting/celltype_seurat.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(object = ATAC_Seurat_withChromVar,  label = TRUE, repel = TRUE,
        label.size = 7, cols = gonads_colors, pt.size = 0.7)
dev.off()


pdf(file = "/home/jovyan/MULTIOME/figures/supporting/sex_seurat.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(object = ATAC_Seurat_withChromVar, group.by = "sex",  repel = TRUE,
        label.size = 7.5, cols = c("plum", "lightskyblue"), pt.size = 0.7)
dev.off()


# convert to CellDataSet format 
supporting_cds <- as.cell_data_set(x = ATAC_Seurat_withChromVar)

supporting_input_cds <- monocle3::detect_genes(supporting_cds)

# Ensure there are no peaks included with zero reads
supporting_input_cds <- supporting_input_cds[Matrix::rowSums(exprs(supporting_input_cds)) != 0,]

# make cicero object
supporting_cicero <- make_cicero_cds(supporting_input_cds, reduced_coordinates = reducedDims(supporting_input_cds)$UMAP)
supporting_cicero

# get the chromosome sizes from the Seurat object
genome <- seqlengths(ATAC_Seurat_withChromVar)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
#genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns_supporting <- run_cicero(supporting_cicero, genomic_coords = genome.df, sample_num = 100)

temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

ccans_supporting <- generate_ccans(conns_supporting)

write.csv(x = conns_supporting, file = paste0(outdir, experiment_prefix, "_conns.csv"))
write.csv(x = ccans_supporting, file = paste0(outdir, experiment_prefix, "_ccans.csv"))

# How many CCANs do I have with a co-accessibility threshold of 0.2? 
length(unique(ccans_supporting$CCAN)) # 3552 --> CCANs are not numbered consecutively!!!!! 

# Add gene info to CCANS peaks 
adata_var = read.csv(paste0(outdir, experiment_prefix, "adata_var_for_cicero.csv"))
adata_var_reduced = adata_var %>% dplyr::select(c("peaks_formatted", "peak_width", "exon", 
                                                  "gene_id", "gene", "gene_name", "annotation", 
                                                  "promoter", "tss_distance", "ENCODE_blacklist"))
adata_var_reduced$Peak <- adata_var_reduced$peaks_formatted
ccans_supporting_anno <- merge(ccans_supporting, adata_var_reduced, by = "Peak")

# Function to retrieve CCANs containing genes of interest and plot it 
retrieve_plot_CCAN <- function(ccons, ccans, gene_of_interest, cutoff){
  # Retrieve CCAN of interest
  ccan_gene <- ccans %>% dplyr::filter(gene_name == gene_of_interest)
  ccan_number <- ccan_gene$CCAN[1]
  ccan_full <- ccans %>% dplyr::filter(CCAN == ccan_number)
  print(ccan_number)
  return(ccan_full)
  # Get parameters for plotting
  chromosome <- str_split(ccan_full$Peak[1], pattern = "-")[[1]][1]
  start <- as.numeric(str_split(ccan_full$Peak[1], pattern = "-")[[1]][2])
  end <- as.numeric(str_split(ccan_full$Peak[nrow(ccan_full)], pattern = "-")[[1]][3])
  # Plot CCAN 
  plot_connections(ccons, chromosome, start, end,
                   gene_model = gene_anno, 
                   coaccess_cutoff = cutoff, 
                   connection_width = .5, 
                   connection_color = "darkgreen",
                   peak_color = "blue",
                   alpha_by_coaccess = TRUE, 
                   gene_model_color = "orange",
                   collapseTranscripts = "longest")
  
}

# Compute average number of peaks per CCAN
avg_peaks_CCAN <- mean(table(ccans_gc$CCAN))

############# Find CCANs enriched in individual cell states ###########
# 1. calculated the fraction of cells of each cell type that have signal at a peak and assumed
# that the distribution of reads per cell across cell types is close to uniform

write.csv(x = ccans_supporting_anno, file = paste0(outdir, experiment_prefix, "_ccans_annotated.csv"))

