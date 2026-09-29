# Somatic cell gene regulatory networks

Code for "Testis somatic cell gene regulatory networks underlie dramatic feedback responses to the depletion of male germ cells"

# Abstract

Mammalian germ line maintenance is mediated through crosstalk between germ cells and their somatic cell environment. We developed a mouse model with one fertile and one infertile testis within the same animal to compare gene expression and chromatin accessibility across fertility states. Sertoli cells and peritubular myoid cells displayed the greatest transcriptional and epigenetic differences in genes related to growth factors, extracellular matrix, and fibrosis. We identified transcription factors, target regions, and genes consistently altered between infertile and fertile testes. Key factors such as Sox9 in Sertoli cells were inferred to control downstream genes including Dhh. As validation, knockout of a predicted enhancer in the first intron of Dhh reduced its expression, whereas knockout of the full gene impaired spermatogenesis. Comparison of fertile and infertile human patients indicated conserved gene expression changes consistent with mouse. This work presents a comprehensive gene regulatory network analysis of somatic cells, providing a systematic map of growth factor regulation in the mammalian testis. 

# Methods

Gene counts were analyzed with Seurat v4 (67) for clustering, integration and differential gene expression and Monocle version 3 (68) for clustering. Filtering cutoffs were set as follows: minimum genes per cell = 100, maximum genes = 8,000, minimum UMIs per cell was set on a per-sample basis after inspecting a rank plot of UMIs per cell ranging between 103 and 103.6 UMI/cell. For cells with ATAC data, filtering criteria of minimum fragments per cell of 100 and max of 10,000 were also applied. Cells with over 30% mitochondrial reads were also excluded. The RNA assay was processed using Seurat’s NormalizeData and ScaleData using 2,000 variable features before being integrated using Seurat’s FindIntegrationAnchors and IntegrateData functions with CCA integration, 1:30 dimensions and otherwise default parameters. The first 30 dimensions were used for integration of the RNA assay of all samples (regardless of associated modality) and joint UMAP creation. Clusters were found using monocle (cluster_cells, resolution = 3 x 10-5, k = 12, partition_qval = 0.05). Differentially expressed genes were identified via Seurat’s FindMarkers function with default parameters. Gene modules were assayed with Seurat’s AddModuleScore function and for the androgen response used MSigDB Hallmark Androgen Response gene set (HALLMARK_ANDROGEN_RESPONSE) obtained via msigdbr. Scores were averaged per biological replicate/testis and compared between paired fertile and infertile samples using a paired Wilcoxon signed-rank test.

For the snMultiome samples, the ATAC modality was processed independently of the RNA assay. After cell type annotation using the RNA assay, peaks were called using Signac (69) with MACS2 separately for each cell type. First a shared peak list was created by merging all samples and then that shared peak list was used to create Seurat objects for each sample individually. FindTopFeatures with min.cutoff = “q0”, RunTFIDF, RunSVD were run. ATAC samples were integrated using FindIntegrationAnchors, reduction = “rlsi” and 2:30 dimensions. Then LSI embeddings were integrated using 1:30 dimensions. Finally, the integrated ATAC assay was joined to the integrated RNA assay by cell name. scRNAseq and snMultiome were integrated as described above via their RNA assays. For WT, W/Wv, Dhh-KO, Dhh-ΔE and Stra8-KO samples, cells were isolated and sequenced as described above, except that WT data was supplemented with additional cells from (33) and (36). All samples, including experimental fertile and infertile, were processed as above except that Harmony integration was used (70) to generate a single object before splitting. Clustering and cell type annotations were re-assigned to the integrated object. 

Differentially accessible peaks were identified via Seurat’s FindMarkers function with default options (i.e., the Wilcoxon rank-sum test). TxDb.Mmusculus.UCSC.mm10.knownGene were used with ChIPseeker (71) to annotate peaks by relative locations to genes (promoter, intergenic, etc). All peak visualization was performed using Signac’s CoveragePlot. Motif accessibility was calculated via ChromVAR(72) using the JASPAR2022(73) set using Core collection and tax_group = vertebrates, and then using Signac’s AddMotifs, FindMotifs and RunChromVAR with the appropriate genomes. Promoter activity was calculated with Signac’s GeneActivity function using 2000bp upstream of the TSS and 0bp downstream. Pearson’s correlations were performed against gene expression on a per-cell basis within each cell type. Ingenuity Pathway Analysis (QIAGEN Inc, https://www.qiagenbioinformatics.com/products/ ingenuity-pathway-analysis) was used for all pathway analyses. In all cases gene lists were used with a minimum fold change cutoff of ± 1.5 and p-adjusted value of ≤ 0.05. 

The GRN analysis was conducted using SCENIC+ v1.0a1 (38) following the standard vignettes with minor modifications. First, custom cisTarget databases were created by extracting the genomic sequences for each peak feature in the snATAC-seq component of the multiome data, using the mm10 genome, and then compiling the databases using the create_cistarget_motif_databases.py tool (https://github.com/aertslab/create_cisTarget_databases). Briefly, 45 topics were selected for the snATAC-seq data following topic modeling and were binarized with both the otsu method and with ntop=3000. The motifs-v10-nr.mgi-m0.00001-o0.0.tbl database was used for the mouse motif annotation. The search space was defined as 0–500 kb. Regulons were filtered with the following parameters: rho_threshold = 0.03, min_regions_per_gene = 0 and min_target_genes = 10. All other parameters were maintained as the defaults. 

To assess ligand-receptor interactions between cell populations, we implemented the Ligand-Receptor Analysis Framework (LIANA) v0.1.12 (74), which infers cell-cell communication using a consensus of 16 cell signaling database resources and 5 CCC methods (Natmi, Connectome, LogFC Mean, SingleCellSignalR, CellphoneDB) with default parameters. Both the germ cell and somatic cell annotations were included across all samples with the RNA assay. We considered the consensus rank generated via Roust Rank Aggregation as the significance p-value to predict the intercellular crosstalk between each pair based on the expression level of known receptors and ligands in the respective clusters and filtered interactions to those with p-value < 0.05. Then, given the SCENIC+ output, for each cell type, we computed the differential regulon activity or gene expression for a) the target gene expression, b) the transcription factor expression, c) the gene-based regulon activity, d) the region-based regulon activity. We then filtered the target genes for those that interact in a significant receptor-ligand interaction as above and where the transcription factor expression and regulon activity are significantly differential between conditions (p-adjusted < 0.05) and with concordant directionality. This yielded a list of regulons and their downstream signaling genes that are inferred to be differentially regulated by condition.  

# Packages

dplyr #version 1.1.4
Seurat #version 5.1.0
SeuratWrappers #version 0.3.1
SeuratData #version 0.2.2.9001
ggplot2 #version 3.5.1
patchwork #version 1.3.0
cowplot #version 1.1.3
viridis #version 0.6.5
monocle3 #version 1.3.7
scCustomize #version 2.1.2
cetcolor #version 0.2.0
stringr #version 1.5.1
pheatmap #version 1.0.12
tidyverse #version 2.0.0
RColorBrewer #version 1.1-3
future #version 1.34.0
Polychrome #version 1.5.1
reshape2 #version 1.4.4
MASS  #version 7.3-61
plotly #version 4.10.4
scales #version 1.3.0
Signac #version 1.15.0
EnsDb.Hsapiens.v86 #version 2.99.0
EnsDb.Mmusculus.v79 #version 2.99.0
EnsDb.Rnorvegicus.v79 #version 2.99.0
BSgenome.Mmusculus.UCSC.mm10 #version 1.4.3
JASPAR2020 #version 0.99.10
TFBSTools #version 1.44.0 
ggforce #version 0.5.0 
eulerr #version 7.0.2
