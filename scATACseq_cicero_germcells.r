####### CICERO ON GERM CELLS ####### 

# Load libraries
library(Signac)
library(Seurat)
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
outdir <- "/nfs/team292/vl6/my_MULTIOME_dir/germcells_apr2021/"
experiment_prefix <- 'germcells_'

# Load Seurat object
ATAC_Seurat_withChromVar = readRDS(file = paste0(outdir, experiment_prefix, "_ATAC_Seurat.rds"))
DefaultAssay(ATAC_Seurat_withChromVar) <- "peaks"
ATAC_Seurat_withChromVar

Idents(ATAC_Seurat_withChromVar) <- ATAC_Seurat_withChromVar@meta.data$clusters
Idents(ATAC_Seurat_withChromVar) <- factor(x = Idents(ATAC_Seurat_withChromVar),
                                    levels = c('PGC' ,'GC_mitotic','oogonia_STRA8','oogonia_meiotic','pre-spermatogonia' ))

gonads_colors = c('#7b9e99', '#bdb380', '#d9abb7', '#edb7b7', '#779eed')

# Add fragment files 
fragments <- CreateFragmentObject(
  path = paste0(outdir, "merged/fragments_germcells.tsv.gz"),
  cells = colnames(ATAC_Seurat_withChromVar), 
  validate.fragments = TRUE
)

Fragments(ATAC_Seurat_withChromVar) <- fragments
ATAC_Seurat_withChromVar


# convert to CellDataSet format 
gc.cds <- as.cell_data_set(x = ATAC_Seurat_withChromVar)

input_cds <- monocle3::detect_genes(gc.cds)

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

# make cicero object
gc.cicero <- make_cicero_cds(input_cds, reduced_coordinates = reducedDims(input_cds)$UMAP)
gc.cicero

# get the chromosome sizes from the Seurat object
genome <- seqlengths(ATAC_Seurat_withChromVar)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
#genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns_gc <- run_cicero(gc.cicero, genomic_coords = genome.df, sample_num = 100)

temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

ccans_gc <- generate_ccans(conns_gc) 

write.csv(x = conns_gc, file = paste0(outdir, experiment_prefix, "_conns.csv"))
write.csv(x = ccans_gc, file = paste0(outdir, experiment_prefix, "_ccans.csv"))

# How many CCANs do I have with a co-accessibility threshold of 0.08? 
length(unique(ccans_gc$CCAN)) # 1920 --> CCANs are not numbered consecutively!!!!! 

# Add gene info to CCANS peaks 
ccans_gc = read.csv(paste0(outdir, experiment_prefix, "_ccans.csv"))
conns_gc = read.csv(paste0(outdir, experiment_prefix, "_conns.csv"))
adata_var = read.csv(paste0(outdir, experiment_prefix, "adata_var_for_cicero.csv"))
adata_var_reduced = adata_var %>% dplyr::select(c("peaks_formatted", "peak_width", "exon", 
                                                  "gene_id", "gene", "gene_name", "annotation", 
                                                  "promoter", "tss_distance", "ENCODE_blacklist"))
adata_var_reduced$Peak <- adata_var_reduced$peaks_formatted
ccans_gc_anno <- merge(ccans_gc, adata_var_reduced, by = "Peak")

# Add CICERO connections to Seurat object
links <- ConnectionsToLinks(conns = conns_gc, ccans = ccans_gc, threshold = 0.1)
Links(ATAC_Seurat_withChromVar) <- links
#Links(ATAC_Seurat_withChromVar_noXYGC) <- links

# Add Annoation 
gene.ranges <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# convert to UCSC style
seqlevelsStyle(gene.ranges) <- "UCSC"
genome(gene.ranges) <- "hg38"

# set gene annotations
Annotation(ATAC_Seurat_withChromVar) <- gene.ranges

# get gene annotation information
Annotation(ATAC_Seurat_withChromVar)

# Try coverage plot 
Idents(ATAC_Seurat_withChromVar) = ATAC_Seurat_withChromVar@meta.data$clusters_sex_final
ATAC_Seurat_withChromVar_noXYGC = subset(ATAC_Seurat_withChromVar, idents = c("PGC_female",
                                                "PGC_male", 'GC_mitotic_female', 'GC_mitotic_male',
                                          'Oogonia_STRA8', 'Oogonia_meiotic', 'Gonocyte'))
Idents(ATAC_Seurat_withChromVar_noXYGC) <- factor(x = Idents(ATAC_Seurat_withChromVar_noXYGC), 
                                                  levels = c("PGC_female",
                                                             "PGC_male", 'GC_mitotic_female', 'GC_mitotic_male',
                                                             'Oogonia_STRA8', 'Oogonia_meiotic', 'Gonocyte'))
Idents(ATAC_Seurat_withChromVar) = ATAC_Seurat_withChromVar@meta.data$clusters
table(Idents(ATAC_Seurat_withChromVar))
pdf(file = "/home/jovyan/MULTIOME/figures/germcells/STRA8_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5) # The height of the plot in inches
CoveragePlot(ATAC_Seurat_withChromVar, region = "chr7-134986138-135268491", window = 2000)
dev.off()

pdf(file = "/home/jovyan/MULTIOME/figures/germcells/ZGLP1_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
CoveragePlot(ATAC_Seurat_withChromVar, region = "chr19-10300000-10320000", window = 2000)
dev.off()

pdf(file = "/home/jovyan/MULTIOME/figures/germcells/TP63_coverage.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
CoveragePlot(ATAC_Seurat_withChromVar, region = "chr3-189629389-189899276", window = 2000)
dev.off()

pdf(file = "/home/jovyan/MULTIOME/figures/germcells/MXD4_coverage.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
CoveragePlot(ATAC_Seurat_withChromVar, region = "chr4-2243432-2268109", window = 2000)
dev.off()

pdf(file = "/home/jovyan/MULTIOME/figures/germcells/DMRTB1_coverage.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
CoveragePlot(ATAC_Seurat_withChromVar, region = "chr1-53454399-53469488", window = 2000)
dev.off()

pdf(file = "/home/jovyan/MULTIOME/figures/germcells/HOXA5_coverage.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
CoveragePlot(ATAC_Seurat_withChromVar, region = "chr7-27139052-27147681", window = 2000)
dev.off()

pdf(file = "/home/jovyan/MULTIOME/figures/germcells/HOXB2_coverage.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
CoveragePlot(ATAC_Seurat_withChromVar, region = "chr17-48539894-48549109 ", window = 2000)
dev.off()





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
  print(start)
  end <- as.numeric(str_split(ccan_full$Peak[nrow(ccan_full)], pattern = "-")[[1]][3])
  print(end)
  # Plot CCAN
  plot_connections(ccons, chromosome, 9003504, 12059913,
                   gene_model = gene_anno, 
                   coaccess_cutoff = cutoff, 
                   connection_width = 1.2, 
                   connection_color = "gray",
                   peak_color = "purple",
                   alpha_by_coaccess = TRUE, 
                   gene_model_color = "blue",
                   collapseTranscripts = "longest")
  
}


###### Genes of interest for germ cells ######
ZGLP1 = retrieve_plot_CCAN(conns_gc, ccans_gc_anno, "CTD-2369P2.10", 0)
pdf(file = "/home/jovyan/MULTIOME/figures/germcells/ZGLP1_CCAN.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
retrieve_plot_CCAN(conns_gc, ccans_gc_anno, "CTD-2369P2.10", 0)
dev.off()

pdf(file = "/home/jovyan/MULTIOME/figures/germcells/STRA8_CCAN.pdf",   # The directory you want to save the file in
    width = 3.5, # The width of the plot in inches
    height = 3.2) # The height of the plot in inches
retrieve_plot_CCAN(conns_gc, ccans_gc_anno, "STRA8", 0)
dev.off()


# Compute average number of peaks per CCAN
avg_peaks_CCAN <- mean(table(ccans_gc$CCAN))

############# Find CCANs enriched in individual cell states ###########
# 1. calculated the fraction of cells of each cell type that have signal at a peak and assumed
# that the distribution of reads per cell across cell types is close to uniform

write.csv(x = ccans_gc_anno, file = paste0(outdir, experiment_prefix, "_ccans_annotated.csv"))

