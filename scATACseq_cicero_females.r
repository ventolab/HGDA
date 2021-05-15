####### CICERO ON FEMALE GONADAL CELLS and INTRAGONADAL CELLS only for PCOS GWAS INTERSECTION ####### 

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
outdir <- "/nfs/team292/vl6/my_MULTIOME_dir/females_apr2021/"
experiment_prefix <- 'females_'

# Load Seurat object for females intra-gonadal cells 
ATAC_Seurat_females = readRDS(file = paste0(outdir, experiment_prefix, "_full.rds"))
DefaultAssay(ATAC_Seurat_females) <- "peaks"
ATAC_Seurat_females
Idents(ATAC_Seurat_females) <- ATAC_Seurat_females@meta.data$cell_type
table(Idents(ATAC_Seurat_females))
Idents(ATAC_Seurat_females) <- factor(x = Idents(ATAC_Seurat_females),
                                      levels = c('Germ cells',
                                                 'coelEpi',  'sPAX8', 'preGC_I_OSR1',
                                                 'ovarianSurf','preGC_II_KITLG', 'preGC_III_GJA1',
                                                  'Oi',  'Gi', 
                                                 'M_prog_ISL1', 'M_MGP', 'M_MullDuct_LGR5', 
                                                 'PV','Epithelial', 'Endothelial', 'Immune', 'Neural' ))

gonads_colors = c('#ff0000', # germs
                  '#366b36', '#4b944a', # sup
                  '#d9439a', '#ffb5ca', '#ff6390',  '#f582c5', # preGC
                  '#fcbdc4', '#94714e', # Gi and Oi
                  '#e6dec8', '#f5eecb',  '#edcba8',
                  '#ffc266', '#d98200', '#e36a1e',  '#ffb485',  # other
                  '#b5b5b5')

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/umap_celltype_seurat.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(object = ATAC_Seurat_females,  label = TRUE, repel = TRUE,
        label.size = 7, cols = gonads_colors)
dev.off()

##### Select intragonadal cells only for PCOS GWAS study 
ATAC_Seurat_females = subset(ATAC_Seurat_females, idents = c('Germ cells', 'coelEpi', 'ovarianSurf', 
                                                             'sPAX8', 'preGC_I_OSR1', 'preGC_II_KITLG', 
                                                             'preGC_III_GJA1', 'Oi', 'Gi', 'PV', 'Immune',
                                                             'Endothelial'))
Idents(ATAC_Seurat_females) = factor(x = Idents(ATAC_Seurat_females), 
                                                 levels = c('Germ cells',
                                                            'coelEpi',  'sPAX8', 'preGC_I_OSR1',
                                                            'ovarianSurf','preGC_II_KITLG', 'preGC_III_GJA1',
                                                            'Oi',  'Gi', 
                                                            'PV', 'Endothelial', 'Immune'))

# convert to Seurat object to CellDataSet format 
females_cds <- as.cell_data_set(x = ATAC_Seurat_females)

input_cds_females <- monocle3::detect_genes(females_cds)

# Ensure there are no peaks included with zero reads
input_cds_females <- input_cds_females[Matrix::rowSums(exprs(input_cds_females)) != 0,] 

# Make cicero object
females_cicero <- make_cicero_cds(input_cds_females, reduced_coordinates = reducedDims(input_cds_females)$UMAP)
females_cicero

# Get the chromosome sizes from the Seurat object
genome_females <- seqlengths(ATAC_Seurat_females)

# Use chromosome 1 to save some time
# Omit this step to run on the whole genome
#genome <- genome[1]

# Convert chromosome sizes to a dataframe
genome_females.df <- data.frame("chr" = names(genome_females), "length" = genome_females)

# Run cicero
conns_females <- run_cicero(females_cicero, genomic_coords = genome_females.df, sample_num = 100)

temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

ccans_females <- generate_ccans(conns_females) # --> COACCESSIBILITY CUTOFF: 0.25

write.csv(x = conns_females, file = paste0(outdir, experiment_prefix, "_conns_intragonadal.csv"))
write.csv(x = ccans_females, file = paste0(outdir, experiment_prefix, "_ccans_intragonadal.csv"))

# How many CCANs do I have with a co-accessibility threshold of 0.25 ? 
length(unique(ccans_females$CCAN)) # 4177 

# Add gene info to CCANS peaks 

adata_var_females = read.csv(paste0(outdir, experiment_prefix, "adata_var_for_cicero.csv"))
library(dplyr)
adata_var_females_reduced = adata_var_females %>% dplyr::select(c("peaks_formatted", "peak_width", "exon", 
                                                              "gene_id", "gene", "gene_name", "annotation", 
                                                              "promoter", "tss_distance", "ENCODE_blacklist"))
adata_var_females_reduced$Peak <- adata_var_females_reduced$peaks_formatted
ccans_females_anno <- merge(ccans_females, adata_var_females_reduced, by = "Peak")

# Compute average number of peaks per CCAN
avg_peaks_CCAN <- mean(table(ccans_females_anno$CCAN)) # 13 

############# Find CCANs enriched in individual cell states ###########
# 1. calculated the fraction of cells of each cell type that have signal at a peak and assumed
# that the distribution of reads per cell across cell types is close to uniform

write.csv(x = ccans_females_anno, file = paste0(outdir, experiment_prefix, "_ccans_annotated_intragonadal.csv"))

#### PLOT CCANs with SNPs #### 
ccans_females_anno = read.csv(paste0(outdir, experiment_prefix, "_ccans_annotated.csv"))

# Plot networks with coverage plots
# Add fragment files 
fragments <- CreateFragmentObject(
  path = paste0(outdir, "merged/fragments.tsv.gz"),
  cells = colnames(ATAC_Seurat_females), 
  validate.fragments = TRUE
)

Fragments(ATAC_Seurat_females) <- fragments
ATAC_Seurat_females

# Add CICERO connections to Seurat object
links_females <- ConnectionsToLinks(conns = conns_females, ccans = ccans_females_anno, threshold = 0)
Idents(ATAC_Seurat_females) <- ATAC_Seurat_females@meta.data$cell_type
table(Idents(ATAC_Seurat_females))
ATAC_Seurat_females_mesenchymal = subset(ATAC_Seurat_females, idents = c('Oi', 'Gi', 'M_MGP', 'M_MullDuct_LGR5',
                                                                         'M_prog_ISL1'))
Idents(ATAC_Seurat_females_mesenchymal) = factor(x = Idents(ATAC_Seurat_females_mesenchymal), 
                                               levels = c('Oi', 'Gi', 'M_prog_ISL1', 'M_MGP', 'M_MullDuct_LGR5'
                                                          ))
table(Idents(ATAC_Seurat_females_mesenchymal))

Links(ATAC_Seurat_females) <- links_females
Links(ATAC_Seurat_females_mesenchymal) <- links_females


# Add Annoation 
gene.ranges <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# convert to UCSC style
seqlevelsStyle(gene.ranges) <- "UCSC"
genome(gene.ranges) <- "hg38"

# set gene annotations
Annotation(ATAC_Seurat_females) <- gene.ranges
Annotation(ATAC_Seurat_females_mesenchymal) <- gene.ranges


library(stringr)
retrieve_plot_CCAN <- function(ccons, ccans, gene_of_interest, cutoff){
  # Retrieve CCAN of interest
  ccan_gene <- ccans %>% dplyr::filter(gene_name == gene_of_interest)
  ccan_number <- ccan_gene$CCAN[1]
  ccan_full <- ccans %>% dplyr::filter(CCAN == ccan_number)
  print(ccan_number)
  # Get parameters for plotting
  chromosome <- str_split(ccan_full$Peak[1], pattern = "-")[[1]][1]
  start <- as.numeric(str_split(ccan_full$Peak[1], pattern = "-")[[1]][2])
  end <- as.numeric(str_split(ccan_full$Peak[nrow(ccan_full)], pattern = "-")[[1]][3])
  print(paste0(chromosome, "-", start, "-", end))
  return(ccan_full)
  
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

ccan_LHX9_female <- retrieve_plot_CCAN(conns_females, ccans_females_anno, "LHX9", 0.25)
ccan_GATA2_female <- retrieve_plot_CCAN(conns_females, ccans_females_anno, "GATA2", 0.25)
ccan_GATA4_female <- retrieve_plot_CCAN(conns_females, ccans_females_anno, "GATA4", 0.25)
ccan_ARX_female <- retrieve_plot_CCAN(conns_females, ccans_females_anno, "ARX", 0.25)
ccan_NR2F1_female <- retrieve_plot_CCAN(conns_females, ccans_females_anno, "NR2F1", 0.25)
ccan_LHX2_female <- retrieve_plot_CCAN(conns_females, ccans_females_anno, "LHX2", 0.25)

females_colors <- c('#fcbdc4', '#94714e', # Gi and Oi
                    '#e6dec8', '#f5eecb',  '#edcba8')
pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/females_GATA2_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_females_mesenchymal, region = "chr3-128476000-128496000", window = 2000)
p & scale_fill_manual(values = females_colors)
dev.off()

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/females_GATA4_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_females_mesenchymal, region = "chr8-11660781-11790330", window = 2000)
p & scale_fill_manual(values = females_colors)
dev.off()

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/females_ARX2_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_females_mesenchymal, region = "chrX-25004418-25021216", window = 2000)
p & scale_fill_manual(values = females_colors)
dev.off()

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/females_NR2F1_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_females_mesenchymal, region = "chr5-93584943-93598705", window = 2000)
p & scale_fill_manual(values = females_colors)
dev.off()

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/females_LHX9_coverage.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 7) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_females_mesenchymal, region = "chr1-197910000-197930000", window = 2000)
p & scale_fill_manual(values = females_colors)
dev.off()

retrieve_plot_CCAN_fromPeak <- function(ccons, ccans, peak_of_interest, cutoff){
  # Retrieve CCAN of interest
  ccan_peak <- ccans %>% dplyr::filter(peaks_formatted == peak_of_interest)
  ccan_number <- ccan_peak$CCAN[1]
  ccan_full <- ccans %>% dplyr::filter(CCAN == ccan_number)
  print(ccan_number)
  #return(ccan_full)
  # Get parameters for plotting
  chromosome <- str_split(ccan_full$Peak[1], pattern = "-")[[1]][1]
  start <- as.numeric(str_split(ccan_full$Peak[1], pattern = "-")[[1]][2])
  print(start)
  end <- as.numeric(str_split(ccan_full$Peak[nrow(ccan_full)], pattern = "-")[[1]][3])
  print(end)
  # Plot CCAN 
  plot_connections(ccons, chromosome, start, end,
                   gene_model = gene_anno, 
                   viewpoint = peak_of_interest,
                   coaccess_cutoff = cutoff, 
                   connection_width = 1.2, 
                   connection_color = "black",
                   peak_color = "black",
                   alpha_by_coaccess = TRUE, 
                   gene_model_color = "blue",
                   collapseTranscripts = "longest")
  
}

############## PCOS GWAS #################

pdf(file = "/home/jovyan/MULTIOME/figures/SNP_PCOS_2_connections.pdf",   # The directory you want to save the file in
    width = 3.5, # The width of the plot in inches
    height = 3.2) # The height of the plot in inches
plot_connections(conns_females, 'chr9', 123872613, 124369810,
                 gene_model = gene_anno, 
                 viewpoint = "chr9-124182485-124183268",
                 coaccess_cutoff = 0, 
                 connection_width = 1.2, 
                 connection_color = "gray",
                 peak_color = "darkgray",
                 alpha_by_coaccess = TRUE, 
                 gene_model_color = "blue",
                 collapseTranscripts = "longest")
dev.off()

pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/SNP_PCOS_connections.pdf",   # The directory you want to save the file in
    width = 3.5, # The width of the plot in inches
    height = 3.2) # The height of the plot in inches
retrieve_plot_CCAN_fromPeak(conns_females, ccans_females_anno, "chr16-52313541-52314547", 0)
dev.off()

female_colors <- c('#ff0000', # germs
                   '#366b36', '#4b944a', # sup
                   '#d9439a', '#ffb5ca', '#ff6390',  '#f582c5', # preGC
                   '#fcbdc4', '#94714e', # Gi and Oi
                   '#ffc266', '#e36a1e',  '#ffb485')
pdf(file = "/home/jovyan/MULTIOME_april2021/figures_females/females_TOX3_coverage.pdf",   # The directory you want to save the file in
    width = 5.5, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches
p <- CoveragePlot(ATAC_Seurat_females, region = "chr16-52068254-52826640", 
             window = 10000, links = FALSE, annotation = FALSE)
p & scale_fill_manual(values = female_colors)
dev.off()


