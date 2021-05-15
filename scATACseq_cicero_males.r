####### CICERO ON MALE GONADAL CELLS ####### 

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
outdir <- "/nfs/team292/vl6/my_MULTIOME_dir/males_apr2021/"
experiment_prefix <- 'males_'

# Load Seurat object for females intra-gonadal cells 
ATAC_Seurat_males = readRDS(file = paste0(outdir, experiment_prefix, "_full.rds"))
DefaultAssay(ATAC_Seurat_males) <- "peaks"
ATAC_Seurat_males
Idents(ATAC_Seurat_males) <- ATAC_Seurat_males@meta.data$cell_type
table(Idents(ATAC_Seurat_males))
Idents(ATAC_Seurat_males) <- factor(x = Idents(ATAC_Seurat_males),
                                      levels = c('Germ cells',
                                                 'coelEpi',  'sPAX8',  'Sertoli',
                                                 'FetalLeydig', 'Ti',  'Gi', 
                                                 'M_prog_ISL1', 'M_MGP',
                                                 'PV','Epithelial', 'Endothelial', 'Immune',  'Erythroid', 'Neural' ))


gonads_colors = c('#ff0000', # germs
                  '#366b36', '#4b944a', # sup
                  '#779eed', '#71a2c7', '#60bddb', '#94714e',
                  '#e6dec8', '#f5eecb', 
                  '#ffc266', '#d98200', '#e36a1e',  '#ffb485', '#ffd919',  # other
                   '#b5b5b5')

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_males/male_umap_celltype_seurat.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(object = ATAC_Seurat_males, label = TRUE, repel = TRUE,
        label.size = 7, cols = gonads_colors)
dev.off()

# convert to Seurat object to CellDataSet format 
cds_males <- as.cell_data_set(x = ATAC_Seurat_males)

input_cds_males <- monocle3::detect_genes(cds_males)

# Ensure there are no peaks included with zero reads
input_cds_males <- input_cds_males[Matrix::rowSums(exprs(input_cds_males)) != 0,]

# Make cicero object
cicero_males <- make_cicero_cds(input_cds_males, reduced_coordinates = reducedDims(input_cds_males)$UMAP)
cicero_males

# Get the chromosome sizes from the Seurat object
genome_males <- seqlengths(ATAC_Seurat_males)

# Use chromosome 1 to save some time
# Omit this step to run on the whole genome
#genome <- genome[1]

# Convert chromosome sizes to a dataframe
genome.df_males <- data.frame("chr" = names(genome_males), "length" = genome_males)

# Run cicero
conns_males <- run_cicero(cicero_males, genomic_coords = genome.df_males, sample_num = 100)

temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

ccans_males <- generate_ccans(conns_males, coaccess_cutoff_override = 0.2) # --> COACCESSIBILITY CUTOFF USED: 0.2

write.csv(x = conns_males, file = paste0(outdir, experiment_prefix, "_conns.csv"))
write.csv(x = ccans_males, file = paste0(outdir, experiment_prefix, "_ccans.csv"))

# How many CCANs do I have with a co-accessibility threshold of ? 
length(unique(ccans_males$CCAN)) # 4588 --> CCANs are not numbered consecutively!!!!! 

# Add gene info to CCANS peaks 
conns_males = read.csv(paste0(outdir, experiment_prefix, "_conns.csv"))
ccans_males = read.csv(paste0(outdir, experiment_prefix, "_ccans.csv"))

adata_var_males = read.csv(paste0(outdir, experiment_prefix, "adata_var_for_cicero.csv"))
library(dplyr)
adata_var_males_reduced = adata_var_males %>% dplyr::select(c("peaks_formatted", "peak_width", "exon", 
                                                  "gene_id", "gene", "gene_name", "annotation", 
                                                  "promoter", "tss_distance", "ENCODE_blacklist"))
adata_var_males_reduced$Peak <- adata_var_males_reduced$peaks_formatted
ccans_males_anno <- merge(ccans_males, adata_var_males_reduced, by = "Peak")

# Compute average number of peaks per CCAN
avg_peaks_CCAN <- mean(table(ccans_males_anno$CCAN)) ### 14 peaks per CCAN

# Plot networks with coverage plots
# Add fragment files 
fragments <- CreateFragmentObject(
  path = paste0(outdir, "merged/fragments.tsv.gz"),
  cells = colnames(ATAC_Seurat_males), 
  validate.fragments = TRUE
)

Fragments(ATAC_Seurat_males) <- fragments
ATAC_Seurat_males

# Add CICERO connections to Seurat object
links_males <- ConnectionsToLinks(conns = conns_males, ccans = ccans_males, threshold = 0)
Idents(ATAC_Seurat_males) <- ATAC_Seurat_males@meta.data$cell_type

# Subset to mesenchymal cells only to look at differential accessibility of TFs distinguishing intragonadal from extragonadal mesenchymal cells
ATAC_Seurat_males_mesenchymal = subset(ATAC_Seurat_males, idents = c('FetalLeydig', 'Ti', 'Gi',  'M_MGP',
                                                                     'M_prog_ISL1'))
Idents(ATAC_Seurat_males_mesenchymal) = factor(x = Idents(ATAC_Seurat_males_mesenchymal), 
                                               levels = c('FetalLeydig', 'Ti', 'Gi',  
                                                          'M_prog_ISL1', 'M_MGP'))
table(Idents(ATAC_Seurat_males_mesenchymal))
Links(ATAC_Seurat_males_mesenchymal) <- links_males


# Add Annoation 
gene.ranges <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# convert to UCSC style
seqlevelsStyle(gene.ranges) <- "UCSC"
genome(gene.ranges) <- "hg38"

# set gene annotations
Annotation(ATAC_Seurat_males) <- gene.ranges
Annotation(ATAC_Seurat_males_mesenchymal) <- gene.ranges

# get gene annotation information
Annotation(ATAC_Seurat_males_mesenchymal)

############# Find CCANs enriched in individual cell states ###########
# 1. calculated the fraction of cells of each cell type that have signal at a peak and assumed
# that the distribution of reads per cell across cell types is close to uniform

write.csv(x = ccans_males_anno, file = paste0(outdir, experiment_prefix, "_ccans_annotated.csv"))

############ Plot CCANs with SNPs ################
library(stringr)
retrieve_CCAN <- function(ccons, ccans, gene_of_interest, cutoff){
  # Retrieve CCAN of interest
  ccan_gene <- ccans %>% dplyr::filter(gene_name == gene_of_interest)
  ccan_number <- ccan_gene$CCAN[1]
  ccan_full <- ccans %>% dplyr::filter(CCAN == ccan_number)
  print(ccan_number)
  # Get parameters for plotting
  chromosome <- str_split(ccan_full$Peak[1], pattern = "-")[[1]][1]
  print(chromosome)
  start <- as.numeric(str_split(ccan_full$Peak[1], pattern = "-")[[1]][2])
  print(start)
  end <- as.numeric(str_split(ccan_full$Peak[nrow(ccan_full)], pattern = "-")[[1]][3])
  print(end)
  return(ccan_full)
  
}


male_cols <- c('#71a2c7', '#60bddb', '#94714e',
               '#e6dec8', '#f5eecb')
#  GATA2: 128470000-128500000
ccan_GATA2 <- retrieve_CCAN(conns_males, ccans_males_anno, "GATA2", 0.25)
pdf(file = "/home/jovyan/MULTIOME_april2021/figures_males/males_GATA2_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_males_mesenchymal, region = "chr3-128476000-128496000", window = 2000)
p & scale_fill_manual(values = male_cols)
dev.off()

ccan_GATA4 <- retrieve_CCAN(conns_males, ccans_males_anno, "GATA4", 0.25)
pdf(file = "/home/jovyan/MULTIOME_april2021/figures_males/males_GATA4_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_males_mesenchymal, region = "chr8-11660781-11790330", window = 2000)
p & scale_fill_manual(values = male_cols)
dev.off()

ccan_ARX <- retrieve_CCAN(conns_males, ccans_males_anno, "ARX", 0.25)
pdf(file = "/home/jovyan/MULTIOME_april2021/figures_males/males_ARX2_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_males_mesenchymal, region = "chrX-25004418-25021216", window = 2000)
p & scale_fill_manual(values = male_cols)
dev.off()

ccan_NR2F1 <- retrieve_CCAN(conns_males, ccans_males_anno, "NR2F1", 0.25)
pdf(file = "/home/jovyan/MULTIOME_april2021/figures_males/males_NR2F1_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_males_mesenchymal, region = "chr5-93584943-93598705", window = 2000)
p & scale_fill_manual(values = male_cols)
dev.off()

# LHX9: 197910000-197940000
ccan_LHX9 <- retrieve_CCAN(conns_males, ccans_males_anno, "LHX9", 0.25)
pdf(file = "/home/jovyan/MULTIOME_april2021/figures_males/males_LHX9_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_males_mesenchymal, region = "chr1-197910000-197930000", window = 2000)
p & scale_fill_manual(values = male_cols)
dev.off()


