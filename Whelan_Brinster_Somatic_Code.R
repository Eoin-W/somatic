# R version 4.4.1

library(dplyr) #version 1.1.4
library(Seurat) #version 5.1.0
library(SeuratWrappers) #version 0.3.1
library(SeuratData) #version 0.2.2.9001
library(ggplot2) #version 3.5.1
library(patchwork) #version 1.3.0
library(cowplot) #version 1.1.3
library(viridis) #version 0.6.5
library(monocle3) #version 1.3.7
library(scCustomize) #version 2.1.2
library(cetcolor) #version 0.2.0
library(stringr) #version 1.5.1
library(pheatmap) #version 1.0.12
library(tidyverse) #version 2.0.0
library(RColorBrewer) #version 1.1-3
library(future) #version 1.34.0
library(Polychrome) #version 1.5.1
library(reshape2) #version 1.4.4
library(MASS)  #version 7.3-61
library(plotly) #version 4.10.4
library(scales) #version 1.3.0
library(Signac) #version 1.15.0
library(EnsDb.Hsapiens.v86) #version 2.99.0
library(EnsDb.Mmusculus.v79) #version 2.99.0
library(EnsDb.Rnorvegicus.v79) #version 2.99.0
library(BSgenome.Mmusculus.UCSC.mm10) #version 1.4.3
library(JASPAR2020) #version 0.99.10
library(TFBSTools) #version 1.44.0 
library(ggforce) #version 0.5.0 
library(eulerr) #version 7.0.2



# the 10x hdf5 file contains both data types. 

inputdata.10x.1 <- Read10X_h5("/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1GFPpos/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.2 <- Read10X_h5("/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1GFPnegA/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.3 <- Read10X_h5("/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1GFPnegB/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.4 <- Read10X_h5("/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1CTL/outs/filtered_feature_bc_matrix.h5")

inputdata.10x.5 <- Read10X_h5("/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2GFPpos/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.6 <- Read10X_h5("/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2GFPnegA/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.7 <- Read10X_h5("/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2GFPnegB/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.8 <- Read10X_h5("/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2CTL/outs/filtered_feature_bc_matrix.h5")

inputdata.10x.9 <- Read10X_h5("/data/ewhelan/somatic_multiome/AMHR2selection/somatic_CTL_AMHR2pos/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.10 <- Read10X_h5("/data/ewhelan/somatic_multiome/AMHR2selection/somatic_GFPtrans_AMHR2pos/outs/filtered_feature_bc_matrix.h5")

inputdata.10x.11 <- Read10X_h5("/data/ewhelan/somatic_multiome/AMHR2selection/somatic_CTL_AMHR2posMACS/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.12 <- Read10X_h5("/data/ewhelan/somatic_multiome/AMHR2selection/somatic_GFPtrans_AMHR2posMACS/outs/filtered_feature_bc_matrix.h5")

inputdata.10x.13 <- Read10X_h5("/data/ewhelan/somatic_multiome/AMHR2selection/somatic_CTL_AMHR2neg/outs/filtered_feature_bc_matrix.h5")
inputdata.10x.14 <- Read10X_h5("/data/ewhelan/somatic_multiome/AMHR2selection/somatic_GFPtrans_AMHR2neg/outs/filtered_feature_bc_matrix.h5")





# extract RNA and ATAC data
rna_counts1 <- inputdata.10x.1$`Gene Expression`
atac_counts1 <- inputdata.10x.1$Peaks

rna_counts2 <- inputdata.10x.2$`Gene Expression`
atac_counts2 <- inputdata.10x.2$Peaks

rna_counts3 <- inputdata.10x.3$`Gene Expression`
atac_counts3 <- inputdata.10x.3$Peaks

rna_counts4 <- inputdata.10x.4$`Gene Expression`
atac_counts4 <- inputdata.10x.4$Peaks

rna_counts5 <- inputdata.10x.5$`Gene Expression`
atac_counts5 <- inputdata.10x.5$Peaks

rna_counts6 <- inputdata.10x.6$`Gene Expression`
atac_counts6 <- inputdata.10x.6$Peaks

rna_counts7 <- inputdata.10x.7$`Gene Expression`
atac_counts7 <- inputdata.10x.7$Peaks

rna_counts8 <- inputdata.10x.8$`Gene Expression`
atac_counts8 <- inputdata.10x.8$Peaks

rna_counts9 <- inputdata.10x.9$`Gene Expression`
atac_counts9 <- inputdata.10x.9$Peaks

rna_counts10 <- inputdata.10x.10$`Gene Expression`
atac_counts10 <- inputdata.10x.10$Peaks

rna_counts11 <- inputdata.10x.11$`Gene Expression`
atac_counts11 <- inputdata.10x.11$Peaks

rna_counts12 <- inputdata.10x.12$`Gene Expression`
atac_counts12 <- inputdata.10x.12$Peaks

rna_counts13 <- inputdata.10x.13$`Gene Expression`
atac_counts13 <- inputdata.10x.13$Peaks

rna_counts14 <- inputdata.10x.14$`Gene Expression`
atac_counts14 <- inputdata.10x.14$Peaks

names.data <- c("multiome1GFPpos", "multiome1GFPnegA", "multiome1GFPnegB", "multiome1CTL", 
                "multiome2GFPpos", "multiome2GFPnegA", "multiome2GFPnegB", "multiome2CTL",
                "multiome3CTLAMHR2posFACS", "multiome3transplantedAMHR2posFACS", 
                "multiome4CTLAMHR2posMACS", "multiome4transplantedAMHR2posMACS",
                "multiome3CTLAMHR2negFACS", "multiome3transplantedAMHR2negFACS")


names.data.all <- c("W54GFP1", "W54GFP2", "W54control", 
                    "WWv2GFP1", "WWv2GFP2", "WWv2control",
                    "WWv3GFPpositive", "WWv3GFPnegative", "WWv3control", 
                    "WWv4GFP1", "WWv4GFP2", "WWv4GFP3", "WWv4control", 
                    "WWv5GFPpositive", "WWv5GFPnegative", "WWv5control", 
                    "WWv6GFPpositive", "WWv6GFPnegative", "WWv6control", 
                    "WWv7GFPpositive", "WWv7GFPnegative", "WWv7control",
                    "multiome1GFPpos", "multiome1GFPnegA", "multiome1GFPnegB", "multiome1CTL", 
                    "multiome2GFPpos", "multiome2GFPnegA", "multiome2GFPnegB", "multiome2CTL",
                    "multiome3CTLAMHR2posFACS", "multiome3transplantedAMHR2posFACS", 
                    "multiome4CTLAMHR2posMACS", "multiome4transplantedAMHR2posMACS",
                    "multiome3CTLAMHR2negFACS", "multiome3transplantedAMHR2negFACS")

table(c(names.data, names.data.all))



list.data <- c(rna_counts1, rna_counts2, rna_counts3,
               rna_counts4, rna_counts5, rna_counts6, 
               rna_counts7, rna_counts8,
               rna_counts9, rna_counts10, rna_counts11, rna_counts12, rna_counts13, rna_counts14)




inflection.point <- c(1:length(list.data))
log_lib_size_at_inflection <- c(1:length(list.data))

par(mfrow = c(3,3)) #This sets up how many graphs per page, in this case 2x1. Change as appropriate (e.g. 9 samples would be c(3,3))


for(y in 1:length(list.data)){
  print(names.data[y])
  expression.df <- as.data.frame(list.data[y])
  
  umi_per_barcode <- colSums(expression.df)
  barcode_rank <- rank(-umi_per_barcode)
  
  log_lib_size <- log10(umi_per_barcode)
  
  o <- order(barcode_rank)
  log_lib_size <- log_lib_size[o]
  barcode_rank <- barcode_rank[o]
  
  rawdiff <- diff(log_lib_size)/diff(barcode_rank)
  inflection <- which(rawdiff == min(rawdiff[500:10000], na.rm=TRUE))
  

  
  plot(barcode_rank, log_lib_size, xlim=c(1,8000),
       pch = 20,
       main = paste(names.data[y], "cells", inflection, sep="_", "UMIs", round(10^log_lib_size[inflection]))
  )
  
  abline(v=inflection, col="blue", lwd=2)
  
  abline(h=log_lib_size[inflection], col="green", lwd=2)
  
  inflection.point[y] <- inflection
  log_lib_size_at_inflection[y] <- log_lib_size[inflection]
}


# Create Seurat object
multiome.mouse.1 <- CreateSeuratObject(counts = rna_counts1, project = "multiome1GFPpos")
multiome.mouse.1[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.1, pattern = "^mt")

multiome.mouse.2 <- CreateSeuratObject(counts = rna_counts2, project = "multiome1GFPnegA")
multiome.mouse.2[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.2, pattern = "^mt")

multiome.mouse.3 <- CreateSeuratObject(counts = rna_counts3, project = "multiome1GFPnegB")
multiome.mouse.3[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.3, pattern = "^mt")

multiome.mouse.4 <- CreateSeuratObject(counts = rna_counts4, project = "multiome1CTL")
multiome.mouse.4[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.4, pattern = "^mt")

multiome.mouse.5 <- CreateSeuratObject(counts = rna_counts5, project = "multiome2GFPpos")
multiome.mouse.5[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.5, pattern = "^mt")

multiome.mouse.6 <- CreateSeuratObject(counts = rna_counts6, project = "multiome2GFPnegA")
multiome.mouse.6[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.6, pattern = "^mt")

multiome.mouse.7 <- CreateSeuratObject(counts = rna_counts7, project = "multiome2GFPnegB")
multiome.mouse.7[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.7, pattern = "^mt")

multiome.mouse.8 <- CreateSeuratObject(counts = rna_counts8, project = "multiome2CTL")
multiome.mouse.8[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.8, pattern = "^mt")

multiome.mouse.9 <- CreateSeuratObject(counts = rna_counts9, project = "multiome3CTL_AMHR2posFACS")
multiome.mouse.9[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.9, pattern = "^mt")

multiome.mouse.10 <- CreateSeuratObject(counts = rna_counts10, project = "multiome3transplanted_AMHR2posFACS")
multiome.mouse.10[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.10, pattern = "^mt")

multiome.mouse.11 <- CreateSeuratObject(counts = rna_counts11, project = "multiome4CTL_AMHR2posMACS")
multiome.mouse.11[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.11, pattern = "^mt")

multiome.mouse.12 <- CreateSeuratObject(counts = rna_counts12, project = "multiome4transplanted_AMHR2posMACS")
multiome.mouse.12[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.12, pattern = "^mt")

multiome.mouse.13 <- CreateSeuratObject(counts = rna_counts13, project = "multiome3CTL_AMHR2negFACS")
multiome.mouse.13[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.13, pattern = "^mt")

multiome.mouse.14 <- CreateSeuratObject(counts = rna_counts14, project = "multiome3transplanted_AMHR2negFACS")
multiome.mouse.14[["percent.mt"]] <- PercentageFeatureSet(multiome.mouse.14, pattern = "^mt")



# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts1 <- StringToGRanges(rownames(atac_counts1), sep = c(":", "-"))
grange.use1 <- seqnames(grange.counts1) %in% standardChromosomes(grange.counts1)
atac_counts1 <- atac_counts1[as.vector(grange.use1), ]

grange.counts2 <- StringToGRanges(rownames(atac_counts2), sep = c(":", "-"))
grange.use2 <- seqnames(grange.counts2) %in% standardChromosomes(grange.counts2)
atac_counts2 <- atac_counts2[as.vector(grange.use2), ]

grange.counts3 <- StringToGRanges(rownames(atac_counts3), sep = c(":", "-"))
grange.use3 <- seqnames(grange.counts3) %in% standardChromosomes(grange.counts3)
atac_counts3 <- atac_counts3[as.vector(grange.use3), ]

grange.counts4 <- StringToGRanges(rownames(atac_counts4), sep = c(":", "-"))
grange.use4 <- seqnames(grange.counts4) %in% standardChromosomes(grange.counts4)
atac_counts4 <- atac_counts4[as.vector(grange.use4), ]

grange.counts5 <- StringToGRanges(rownames(atac_counts5), sep = c(":", "-"))
grange.use5 <- seqnames(grange.counts5) %in% standardChromosomes(grange.counts5)
atac_counts5 <- atac_counts5[as.vector(grange.use5), ]

grange.counts6 <- StringToGRanges(rownames(atac_counts6), sep = c(":", "-"))
grange.use6 <- seqnames(grange.counts6) %in% standardChromosomes(grange.counts6)
atac_counts6 <- atac_counts6[as.vector(grange.use6), ]

grange.counts7 <- StringToGRanges(rownames(atac_counts7), sep = c(":", "-"))
grange.use7 <- seqnames(grange.counts7) %in% standardChromosomes(grange.counts7)
atac_counts7 <- atac_counts7[as.vector(grange.use7), ]

grange.counts8 <- StringToGRanges(rownames(atac_counts8), sep = c(":", "-"))
grange.use8 <- seqnames(grange.counts8) %in% standardChromosomes(grange.counts8)
atac_counts8 <- atac_counts8[as.vector(grange.use8), ]

grange.counts9 <- StringToGRanges(rownames(atac_counts9), sep = c(":", "-"))
grange.use9 <- seqnames(grange.counts9) %in% standardChromosomes(grange.counts9)
atac_counts9 <- atac_counts9[as.vector(grange.use9), ]

grange.counts10 <- StringToGRanges(rownames(atac_counts10), sep = c(":", "-"))
grange.use10 <- seqnames(grange.counts10) %in% standardChromosomes(grange.counts10)
atac_counts10 <- atac_counts10[as.vector(grange.use10), ]

grange.counts11 <- StringToGRanges(rownames(atac_counts11), sep = c(":", "-"))
grange.use11 <- seqnames(grange.counts11) %in% standardChromosomes(grange.counts11)
atac_counts11 <- atac_counts11[as.vector(grange.use11), ]

grange.counts12 <- StringToGRanges(rownames(atac_counts12), sep = c(":", "-"))
grange.use12 <- seqnames(grange.counts12) %in% standardChromosomes(grange.counts12)
atac_counts12 <- atac_counts12[as.vector(grange.use12), ]

grange.counts13 <- StringToGRanges(rownames(atac_counts13), sep = c(":", "-"))
grange.use13 <- seqnames(grange.counts13) %in% standardChromosomes(grange.counts13)
atac_counts13 <- atac_counts13[as.vector(grange.use13), ]

grange.counts14 <- StringToGRanges(rownames(atac_counts14), sep = c(":", "-"))
grange.use14 <- seqnames(grange.counts14) %in% standardChromosomes(grange.counts14)
atac_counts14 <- atac_counts14[as.vector(grange.use14), ]




annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "GRCm38"


atac.frags <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1GFPpos/outs/atac_fragments.tsv.gz")
atac.frags2 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1GFPnegA/outs/atac_fragments.tsv.gz")
atac.frags3 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1GFPnegB/outs/atac_fragments.tsv.gz")
atac.frags4 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Somatic1/multiomeSomatic1CTL/outs/atac_fragments.tsv.gz")
atac.frags5 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2GFPpos/outs/atac_fragments.tsv.gz")
atac.frags6 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2GFPnegA/outs/atac_fragments.tsv.gz")
atac.frags7 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2GFPnegB/outs/atac_fragments.tsv.gz")
atac.frags8 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/Multiome10/multiomeSomatic2CTL/outs/atac_fragments.tsv.gz")
atac.frags9 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/AMHR2selection/somatic_CTL_AMHR2pos/outs/atac_fragments.tsv.gz")
atac.frags10 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/AMHR2selection/somatic_GFPtrans_AMHR2pos/outs/atac_fragments.tsv.gz")
atac.frags11 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/AMHR2selection/somatic_CTL_AMHR2posMACS/outs/atac_fragments.tsv.gz")
atac.frags12 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/AMHR2selection/somatic_GFPtrans_AMHR2posMACS/outs/atac_fragments.tsv.gz")
atac.frags13 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/AMHR2selection/somatic_CTL_AMHR2neg/outs/atac_fragments.tsv.gz")
atac.frags14 <- CreateFragmentObject(path = "/data/ewhelan/somatic_multiome/AMHR2selection/somatic_GFPtrans_AMHR2neg/outs/atac_fragments.tsv.gz")



#to make the counts matrices we need a merged file so let's make that first. 

chrom_assay1 <- CreateChromatinAssay(
  counts = atac_counts1,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay2 <- CreateChromatinAssay(
  counts = atac_counts2,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags2,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay3 <- CreateChromatinAssay(
  counts = atac_counts3,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags3,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay4 <- CreateChromatinAssay(
  counts = atac_counts4,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags4,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay5 <- CreateChromatinAssay(
  counts = atac_counts5,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags5,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay6 <- CreateChromatinAssay(
  counts = atac_counts6,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags6,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay7 <- CreateChromatinAssay(
  counts = atac_counts7,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags7,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay8 <- CreateChromatinAssay(
  counts = atac_counts8,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags8,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay9 <- CreateChromatinAssay(
  counts = atac_counts9,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags9,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay10 <- CreateChromatinAssay(
  counts = atac_counts10,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags10,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay11 <- CreateChromatinAssay(
  counts = atac_counts11,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags11,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay12 <- CreateChromatinAssay(
  counts = atac_counts12,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags12,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay13 <- CreateChromatinAssay(
  counts = atac_counts13,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags13,
  # min.cells = 10,
  annotation = annotations
)

chrom_assay14 <- CreateChromatinAssay(
  counts = atac_counts14,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags14,
  # min.cells = 10,
  annotation = annotations
)


multiome.mouse.1[["ATAC"]] <- chrom_assay1
multiome.mouse.2[["ATAC"]] <- chrom_assay2
multiome.mouse.3[["ATAC"]] <- chrom_assay3
multiome.mouse.4[["ATAC"]] <- chrom_assay4
multiome.mouse.5[["ATAC"]] <- chrom_assay5
multiome.mouse.6[["ATAC"]] <- chrom_assay6
multiome.mouse.7[["ATAC"]] <- chrom_assay7
multiome.mouse.8[["ATAC"]] <- chrom_assay8

multiome.mouse.9[["ATAC"]] <- chrom_assay9
multiome.mouse.10[["ATAC"]] <- chrom_assay10
multiome.mouse.11[["ATAC"]] <- chrom_assay11
multiome.mouse.12[["ATAC"]] <- chrom_assay12
multiome.mouse.13[["ATAC"]] <- chrom_assay13
multiome.mouse.14[["ATAC"]] <- chrom_assay14

list.of.seurats <- c(multiome.mouse.1,multiome.mouse.2, multiome.mouse.3, multiome.mouse.4, 
                     multiome.mouse.5, multiome.mouse.6, multiome.mouse.7, multiome.mouse.8,
                     multiome.mouse.9, multiome.mouse.10, multiome.mouse.11, multiome.mouse.12,
                     multiome.mouse.13, multiome.mouse.14)

for(current.loop in 1:length(list.of.seurats)){
  # current.loop <- 1
  current.seurat <- list.of.seurats[[current.loop]]
  cell.names <- Cells(current.seurat)
  current.name <- names.data[current.loop]
  cell.names.fixed <- sub("-1$", "", cell.names)
  cell.names.fixed <- paste0(current.name, "_", cell.names.fixed)
  current.seurat.fixed <- RenameCells(object = current.seurat, new.names = cell.names.fixed)
  print(head(Cells(current.seurat.fixed)))
  print(table((current.seurat.fixed$orig.ident)))
  current.seurat.fixed -> list.of.seurats[[current.loop]]
}


combined <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

setwd("/data/ewhelan/Robjects")
save(combined, file="combined.Robj")


VlnPlot(combined, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = T, split.by = "orig.ident", pt.size = 0) + NoLegend()


#now make new counts objects based on the combined file

DefaultAssay(combined) <- "ATAC"

counts <- FeatureMatrix(
  fragments = atac.frags,
  features = granges(combined)
)

counts2 <- FeatureMatrix(
  fragments = atac.frags2,
  features = granges(combined)
)

counts3 <- FeatureMatrix(
  fragments = atac.frags3,
  features = granges(combined)
)

counts4 <- FeatureMatrix(
  fragments = atac.frags4,
  features = granges(combined)
)

counts5 <- FeatureMatrix(
  fragments = atac.frags5,
  features = granges(combined)
)

counts6 <- FeatureMatrix(
  fragments = atac.frags6,
  features = granges(combined)
)

counts7 <- FeatureMatrix(
  fragments = atac.frags7,
  features = granges(combined)
)

counts8 <- FeatureMatrix(
  fragments = atac.frags8,
  features = granges(combined)
)

counts9 <- FeatureMatrix(
  fragments = atac.frags9,
  features = granges(combined)
)

counts10 <- FeatureMatrix(
  fragments = atac.frags10,
  features = granges(combined)
)

counts11 <- FeatureMatrix(
  fragments = atac.frags11,
  features = granges(combined)
)

counts12 <- FeatureMatrix(
  fragments = atac.frags12,
  features = granges(combined)
)

counts13 <- FeatureMatrix(
  fragments = atac.frags13,
  features = granges(combined)
)

counts14 <- FeatureMatrix(
  fragments = atac.frags14,
  features = granges(combined)
)

save(counts, file="counts1.Robj")
save(counts2, file="counts2.Robj")
save(counts3, file="counts3.Robj")
save(counts4, file="counts4.Robj")
save(counts5, file="counts5.Robj")
save(counts6, file="counts6.Robj")
save(counts7, file="counts7.Robj")
save(counts8, file="counts8.Robj")
save(counts9, file="counts9.Robj")
save(counts10, file="counts10.Robj")
save(counts11, file="counts11.Robj")
save(counts12, file="counts12.Robj")
save(counts13, file="counts13.Robj")
save(counts14, file="counts14.Robj")



atac.assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags,
  # min.cells = 10,
  annotation = annotations
)


atac.assay2 <- CreateChromatinAssay(
  counts = counts2,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags2,
  # min.cells = 10,
  annotation = annotations
)


atac.assay3 <- CreateChromatinAssay(
  counts = counts3,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags3,
  # min.cells = 10,
  annotation = annotations
)


atac.assay4 <- CreateChromatinAssay(
  counts = counts4,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags4,
  # min.cells = 10,
  annotation = annotations
)


atac.assay5 <- CreateChromatinAssay(
  counts = counts5,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags5,
  # min.cells = 10,
  annotation = annotations
)


atac.assay6 <- CreateChromatinAssay(
  counts = counts6,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags6,
  # min.cells = 10,
  annotation = annotations
)


atac.assay7 <- CreateChromatinAssay(
  counts = counts7,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags7,
  # min.cells = 10,
  annotation = annotations
)


atac.assay8 <- CreateChromatinAssay(
  counts = counts8,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags8,
  # min.cells = 10,
  annotation = annotations
)

atac.assay9 <- CreateChromatinAssay(
  counts = counts9,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags9,
  # min.cells = 10,
  annotation = annotations
)

atac.assay10 <- CreateChromatinAssay(
  counts = counts10,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags10,
  # min.cells = 10,
  annotation = annotations
)

atac.assay11 <- CreateChromatinAssay(
  counts = counts11,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags11,
  # min.cells = 10,
  annotation = annotations
)

atac.assay12 <- CreateChromatinAssay(
  counts = counts12,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags12,
  # min.cells = 10,
  annotation = annotations
)

atac.assay13 <- CreateChromatinAssay(
  counts = counts13,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags13,
  # min.cells = 10,
  annotation = annotations
)

atac.assay14 <- CreateChromatinAssay(
  counts = counts14,
  sep = c(":", "-"),
  genome = 'GRCm38',
  fragments = atac.frags14,
  # min.cells = 10,
  annotation = annotations
)



mouse1.atac <- CreateSeuratObject(counts = atac.assay, assay = "peaks", project = "multiome1GFPpos")
mouse2.atac <- CreateSeuratObject(counts = atac.assay2, assay = "peaks", project = "multiome1GFPnegA")
mouse3.atac <- CreateSeuratObject(counts = atac.assay3, assay = "peaks", project = "multiome1GFPnegB")
mouse4.atac <- CreateSeuratObject(counts = atac.assay4, assay = "peaks", project = "multiome1CTL")
mouse5.atac <- CreateSeuratObject(counts = atac.assay5, assay = "peaks", project = "multiome2GFPpos")
mouse6.atac <- CreateSeuratObject(counts = atac.assay6, assay = "peaks", project = "multiome2GFPnegA")
mouse7.atac <- CreateSeuratObject(counts = atac.assay7, assay = "peaks", project = "multiome2GFPnegB")
mouse8.atac <- CreateSeuratObject(counts = atac.assay8, assay = "peaks", project = "multiome2CTL")
mouse9.atac <- CreateSeuratObject(counts = atac.assay9, assay = "peaks", project = "multiome3CTL_AMHR2posFACS")
mouse10.atac <- CreateSeuratObject(counts = atac.assay10, assay = "peaks", project = "multiome3transplanted_AMHR2posFACS")
mouse11.atac <- CreateSeuratObject(counts = atac.assay11, assay = "peaks", project = "multiome4CTL_AMHR2posMACS")
mouse12.atac <- CreateSeuratObject(counts = atac.assay12, assay = "peaks", project = "multiome4transplanted_AMHR2posMACS")
mouse13.atac <- CreateSeuratObject(counts = atac.assay13, assay = "peaks", project = "multiome3CTL_AMHR2negFACS")
mouse14.atac <- CreateSeuratObject(counts = atac.assay14, assay = "peaks", project = "multiome3transplanted_AMHR2negFACS")


list.of.atac <- c(mouse1.atac, mouse2.atac, mouse3.atac, mouse4.atac, mouse5.atac, mouse6.atac, mouse7.atac, mouse8.atac,
                  mouse9.atac, mouse10.atac, mouse11.atac, mouse12.atac, mouse13.atac, mouse14.atac)

for(current.sample in 1:(length(list.of.atac))){
  current.atac <- list.of.atac[[current.sample]]
  print(cat(names.data[current.sample], " before filter cells = ", length(Cells(current.atac))))
  current.atac <- subset(
    x = current.atac,
    subset = nCount_peaks < 1* 10^5 &
      nCount_peaks > 100
  )
  print(cat(names.data[current.sample], " after filter cells = ", length(Cells(current.atac))))
  current.atac <- FindTopFeatures(current.atac, min.cutoff = 10)
  current.atac <- RunTFIDF(current.atac)
  current.atac <- RunSVD(current.atac)
  current.atac -> list.of.atac[[current.sample]]
}

### integrate ATAC samples ###

merge.mouse.atac <- merge(
  x = list.of.atac[[1]],
  y = list.of.atac[2:length(list.of.atac)]
)

save(merge.mouse.atac, file="merge.mouse.atac.Robj")
# load("merge.mouse.atac.Robj")
#can run this if you want to see the merged file. 
merge.mouse.atac <- FindTopFeatures(merge.mouse.atac, min.cutoff = 10)
merge.mouse.atac <- RunTFIDF(merge.mouse.atac)
merge.mouse.atac <- RunSVD(merge.mouse.atac)
merge.mouse.atac <- RunUMAP(merge.mouse.atac, reduction = "lsi", dims = 2:30)
DimPlot(merge.mouse.atac, group.by = "orig.ident", raster = FALSE)

# 
# 
# list.of.mouse.multiome <- list(mouse1.atac, mouse2.atac, mouse3.atac)

DefaultAssay(merge.mouse.atac) <- "peaks"

integration.anchors <- FindIntegrationAnchors(
  object.list = list.of.atac,
  anchor.features = rownames(merge.mouse.atac),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merge.mouse.atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
DimPlot(integrated, raster = FALSE)

save(integrated, file = "integrated.Robj")


###
###
### RNA INTEGRATION


### STANDARD INTEGRATION ###

DefaultAssay(multiome.mouse.1) <- "RNA"
multiome.mouse.1[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.2) <- "RNA"
multiome.mouse.2[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.3) <- "RNA"
multiome.mouse.3[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.4) <- "RNA"
multiome.mouse.4[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.5) <- "RNA"
multiome.mouse.5[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.6) <- "RNA"
multiome.mouse.6[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.7) <- "RNA"
multiome.mouse.7[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.8) <- "RNA"
multiome.mouse.8[["ATAC"]] <- NULL

DefaultAssay(multiome.mouse.9) <- "RNA"
multiome.mouse.9[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.10) <- "RNA"
multiome.mouse.10[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.11) <- "RNA"
multiome.mouse.11[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.12) <- "RNA"
multiome.mouse.12[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.13) <- "RNA"
multiome.mouse.13[["ATAC"]] <- NULL
DefaultAssay(multiome.mouse.14) <- "RNA"
multiome.mouse.14[["ATAC"]] <- NULL

current.list <- c(multiome.mouse.1, multiome.mouse.2, multiome.mouse.3,
                  multiome.mouse.4, multiome.mouse.5, multiome.mouse.6,
                  multiome.mouse.7, multiome.mouse.8, 
                  multiome.mouse.9, multiome.mouse.10, multiome.mouse.11, 
                  multiome.mouse.12, multiome.mouse.13, multiome.mouse.14)


names.data <- c(                "multiome1GFPpos", "multiome1GFPnegA", "multiome1GFPnegB", "multiome1CTL", 
                                "multiome2GFPpos", "multiome2GFPnegA", "multiome2GFPnegB", "multiome2CTL", 
                                "multiome3CTL_AMHR2posFACS", "multiome3transplanted_AMHR2posFACS", 
                                "multiome4CTL_AMHR2posMACS", "multiome4transplanted_AMHR2posMACS",
                                "multiome3CTL_AMHR2negFACS", "multiome3transplanted_AMHR2negFACS")
sample.names.simple <- c(                "multiome1GFPpos", "multiome1GFPnegA", "multiome1GFPnegB", "multiome1CTL", 
                                         "multiome2GFPpos", "multiome2GFPnegA", "multiome2GFPnegB", "multiome2CTL", 
                                         "multiome3CTL_AMHR2pos", "multiome3transplanted_AMHR2pos", 
                                         "multiome4CTL_AMHR2pos", "multiome4transplanted_AMHR2pos",
                                         "multiome3CTL_AMHR2neg", "multiome3transplanted_AMHR2neg")
replicate <- c("mouse5", "mouse5", "mouse5", "mouse5", "mouse6", "mouse6", "mouse6", "mouse6",
               "mouse7", "mouse7", "mouse8", "mouse8", "mouse7", "mouse7")
selection <- c("GFPpos", "GFPneg", "GFPneg","GFPneg", "GFPpos", "GFPneg", "GFPneg","GFPneg", 
               "AMHR2pos", "AMHR2pos", "AMHR2neg", "AMHR2neg", "AMHR2pos", "AMHR2pos")
treatment <- c("GFPtransplanted", "GFPtransplanted", "GFPtransplanted", "Control", "GFPtransplanted", "GFPtransplanted", "GFPtransplanted", "Control",
               "Control", "GFPtransplanted", "Control", "GFPtransplanted", "Control", "GFPtransplanted")



for(current.sample in 1:length(current.list)){
  print(current.sample)
  current.seurat <- current.list[[current.sample]]
  # current.seurat[['species']] <- "mouse"
  current.seurat[['sample']] <- sample.names.simple[current.sample]
  current.seurat[['technology']] <- "GEX+ATAC"
  current.seurat[['selection']] <- selection[current.sample]
  current.seurat[['treatment']] <- treatment[current.sample]
  DefaultAssay(current.seurat) <- "RNA"
  current.list[[current.sample]] <- current.seurat
}


for (i in 1:length(x = current.list)) {
  #i=1
  current.list[[i]] <- NormalizeData(object = current.list[[i]],
                                     #normalization.method = "LogNormalize",
                                     #scale.factor = 10000,
                                     #margin = 1,
                                     verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]],
                                            #selection.method =  "vst", #"mean.var.plot",
                                            nfeatures = 2000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list,
                                          #reference = 1,
                                          dims = 1:30
)
save(current.anchors, file="current.anchors.Robj")

current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:30)
save(current.integrated, file="current.integrated.Robj")

DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = 30,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:30)

DimPlot(object = current.integrated, reduction = "umap", raster = FALSE, group.by = "orig.ident")

mouse.combined <- current.integrated
save(mouse.combined, file = "mouse.combined.Robj")

DefaultAssay(mouse.combined) <- "RNA"
VlnPlot(mouse.combined, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 0,
        group.by = "orig.ident" )

#have a look at undifferentiated spermatogonia
FeaturePlot(mouse.combined, features = c("Etv5", "Id4", "Sdc4", "Ret"), min.cutoff = 0, order = T)

#merge the datasets

cells.RNA <- (Cells(mouse.combined))
cells.ATAC <- (Cells(integrated))
cells.overlap <- (intersect(cells.RNA, cells.ATAC))
cat("cells.RNA", length(cells.RNA), 
    "cells.ATAC", length(cells.ATAC), 
    "cells.overlap", length(cells.overlap), 
    sep = " ")

mouse.multiome.integrated <- subset(mouse.combined, cells = cells.overlap)
mouse.multiome.ATAC <- subset(integrated, cells = cells.overlap)
mouse.multiome.integrated[["peaks"]] <- mouse.multiome.ATAC[["peaks"]]
mouse.multiome.integrated@reductions$lsi <- mouse.multiome.ATAC@reductions$integrated_lsi

mouse.multiome.integrated <- subset(x = mouse.multiome.integrated, 
                                    # subset = nFeature_RNA > 10^3.2 
                                    
                                    subset = percent.mt < 30 
                                    # & nCount_RNA > 10^3.2 
                                    & nCount_peaks < 10^5 &
                                      nCount_peaks > 1000)
# & nFeature_RNA < max.features) 

VlnPlot(mouse.multiome.integrated, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        pt.size = 0,
        group.by = "orig.ident" )

n.pcs <- 25

DefaultAssay(mouse.multiome.integrated) <- "peaks"
mouse.multiome.integrated <- FindMultiModalNeighbors(mouse.multiome.integrated, reduction.list = list("pca", "lsi"), dims.list = list(1:n.pcs, 2:n.pcs))
mouse.multiome.integrated <- RunUMAP(mouse.multiome.integrated, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")
DefaultAssay(mouse.multiome.integrated) <- "RNA"
mouse.multiome.integrated <- RunUMAP(object = mouse.multiome.integrated, reduction.name = "umap.rna", reduction = "pca", dims = 1:n.pcs)
DefaultAssay(mouse.multiome.integrated) <- "peaks"
mouse.multiome.integrated <- RunUMAP(object = mouse.multiome.integrated, reduction.name = "umap.atac", reduction = "lsi", dims = 2:n.pcs)


save(mouse.multiome.integrated, file = "mouse.multiome.integrated.2023.08.08.Robj")
load("mouse.multiome.integrated.2023.08.08.Robj")

p1 <- DimPlot(mouse.multiome.integrated, reduction = "umap.rna", group.by = "seurat_clusters", label = T)+ggtitle("RNA") #+NoLegend()
p2 <- DimPlot(mouse.multiome.integrated, reduction = "umap.atac", group.by = "seurat_clusters", label = T)+ggtitle("ATAC")#+NoLegend()
p3 <- DimPlot(mouse.multiome.integrated, reduction = "umap.wnn", group.by = "seurat_clusters", label = T)+ggtitle("merged (WNN)")#+NoLegend()

# DimPlot(mouse.multiome.integrated$, reduction = "umap.wnn", group.by = "seurat_clusters", label = T, split.by = "orig.ident")

# FeaturePlot(combined3, reduction = "umap.wnn", features = c("Etv5", "Id4"),
#             # max.cutoff = 2.5,
#             order = T)& scale_color_viridis(option="magma", direction = -1)


print(p1 | p2 | p3)

DefaultAssay(mouse.multiome.integrated) <- "RNA"

gene.of.interest <- c("Sox9")
p4 <- FeaturePlot(mouse.multiome.integrated, features = gene.of.interest, min.cutoff = 0, order = T, reduction = "umap.rna")
p5 <- FeaturePlot(mouse.multiome.integrated, features = gene.of.interest, min.cutoff = 0, order = T, reduction = "umap.atac")
p6 <- FeaturePlot(mouse.multiome.integrated, features = gene.of.interest, min.cutoff = 0, order = T, reduction = "umap.wnn")


print(p4 | p5 | p6)

FeaturePlot(mouse.multiome.integrated, features = c("Prm1", "Tnp1"), split.by = "treatment", order = T, reduction = "umap.atac")
FeaturePlot(somatic.integrated, features = c("Prm1", "Tnp1"), split.by = "treatment", min.cutoff = 0, order = F)

genes.spermatogonia <- c("Sdc4", "Etv5", "Uchl1", "Sycp1", "Crabp1", "Zbtb16", "Sall4") #spermatogonia
genes.diff.spermatogonia <- c("Sohlh1", "Kit", "Prdm9", "Stra8") #differentiating cells
genes.spermatocytes1 <- c("Piwil1", "Pttg1", "Insl6", "Spag6", "Tbpl1", "Tex101", "Spo11", "Mei1", "Piwil2")
genes.spermatocytes2 <- c("Sycp3", "Tdrd5", "Tbpl1", "Hormad1", "H2afx", "Mns1") #spermatocytes
genes.spermatids <- c("Acrv1", "Spaca1", "Tsga8", "Tssk1") #spermatids
genes.elongating <- c("Prm1", "Prm2", "Tnp1", "Tnp2", "Hspa1l", "Pgk2") #elongating spermatids
genes.somatic <- c("Cyp17a1", "Acta2", "Vwf", "Cd74", "Sox9", "Thy1")
#Leydig, myoid, endothelial, macrophage, sertoli, innate lymph
genes.leydig <- c("Cyp17a1", "Cyp11a1", "Star", "Hsd3b1")
genes.myoid <- c("Acta2", "Myh11", "Myl6", "Pdgfrb")
genes.endothelial <- c("Vwf", "Tie1", "Tek")
genes.macrophage <- c("Apoe", "Dab2", "Cd74", "Adgre1", "Amhr2")
genes.sertoli <- c("Clu", "Amhr2", "Sox9", "Ctsl", "Rhox8")
genes.innatelymph <- c("Id2", "Il7r", "Rora", "Thy1", "Ccl5", "Cd52")

gene.list.SSC <- c("Etv5", "Id4", "Nanos3", "Gfra1", "Sdc4", "Zbtb16")

gene.SSC.progenitor.early.diff <- c("Gfra1", "Etv5", "Sdc4", "Zbtb16", "Kit", "Sohlh1")
late.diff.prelep.lep <- c("Stra8", "Prdm9","Mei1", "Sycp3", "Spo11", "Hormad1")
pachytene.diplotene.early.round <- c("Piwil1", "Topaz1", "Mlh3", "Ccna1", "Acvr1", "Spaca1")
late.round.spermatids <- c("Pgk2", "Tssk1", "Prm1", "Tnp1")
# 
# Hermann1 <- c("Esrp1", "Nanos2", "Zic1", "Nanos3", "Sox3", "Upp1", "Lmo1", "Loxl2", "Galnt12")
# Hermann2 <- c("Pth1r", "Snx16", "Tbx1", "Dmc1", "Meiob", "Rad51ap2", "Adad2", "B3galnt1", "Ccnb1ip1")
# Hermann3 <- c("Dyx1c1", "4933414I15Rik", "Gm35584", "4930511A02Rik", "Prss42", "Rhog", "Gm10354", "Ssxb2", "Speer4e")
# Hermann4 <- c("4930513O06Rik", "Acot10", "Adam29", "1700027A15Rik", "1700080E11Rik", "Cyp2a12")
# 
# leptotene <- c("Mlh1", "Spo11", "Terf1", "Ube2b", "Xrcc5") #from gene ontology
# totalseqRat <- c("Itgb3", "Thy1", "Cd27", "Itgb1")
# 
# FGFs <- c("Fgf1", "Fgf2", "Fgf3", "Fgf4", "Fgf5",
#           "Fgf6", "Fgf7", "Fgf8", "Fgf9", "Fgf10",
#           "Fgf11", "Fgf12", "Fgf13", "Fgf14", #"Fgf15",
#           "Fgf16", "Fgf17", "Fgf18", "Fgf19", "Fgf20",
#           "Fgf21", "Fgf22", "Fgf23")

# germ.cell.various <- (c("Dazl", "Kit", "Ddx4", "Sycp3", "Prm1", "Tex101"))

T.cells <- c("Cd3e", "Cd8a", "Cd4", "Cd28", "Cd45", "Cd25", "Cd56") 
B.cells <- c("Cd19", "Cd79a")
Tregs <- c("Foxp3")



genes.spermatogonia <- c("Sdc4", "Etv5", "Uchl1", "Sycp1", "Crabp1", "Zbtb16", "Sall4") #spermatogonia
genes.diff.spermatogonia <- c("Sohlh1", "Kit", "Prdm9", "Stra8") #differentiating cells
genes.spermatocytes1 <- c("Piwil1", "Pttg1", "Insl6", "Spag6", "Tbpl1", "Tex101", "Spo11", "Mei1", "Piwil2")
genes.spermatocytes2 <- c("Sycp3", "Tdrd5", "Tbpl1", "Hormad1", "H2afx", "Mns1") #spermatocytes
genes.spermatids <- c("Acrv1", "Spaca1", "Tsga8", "Tssk1") #spermatids
genes.elongating <- c("Prm1", "Prm2", "Tnp1", "Tnp2", "Hspa1l", "Pgk2") #elongating spermatids
genes.somatic <- c("Cyp17a1", "Acta2", "Vwf", "Cd74", "Sox9", "Thy1")
#Leydig, myoid, endothelial, macrophage, sertoli, innate lymph
genes.leydig <- c("Cyp17a1", "Cyp11a1", "Star", "Hsd3b1")
genes.myoid <- c("Acta2", "Myh11", "Myl6", "Pdgfrb")
genes.endothelial <- c("Vwf", "Tie1", "Tek")
genes.macrophage <- c("Apoe", "Dab2", "Cd74", "Adgre1", "Amhr2")
genes.sertoli <- c("Clu", "Amhr2", "Sox9", "Ctsl", "Rhox8")
genes.innatelymph <- c("Id2", "Il7r", "Rora", "Thy1", "Ccl5", "Cd52")

DefaultAssay(mouse.multiome.integrated) <- "RNA"
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.spermatogonia),name="Spermatogonia")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.diff.spermatogonia),name="Diff.Spg")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.spermatocytes2),name="Spermatocytes")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.spermatids),name="Round")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.elongating),name="Elongating")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.leydig),name="Leydig")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.myoid),name="Myoid")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.endothelial),name="Endothelial")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.sertoli),name="Sertoli")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(genes.macrophage),name="Macrophages")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(T.cells),name="Tcells")
mouse.multiome.integrated <- AddModuleScore(mouse.multiome.integrated,features = list(B.cells),name="Bcells")


plot_grid(ncol = 5,
          DimPlot(mouse.multiome.integrated, label = T, reduction = "umap.atac") + NoLegend(),
          FeaturePlot(mouse.multiome.integrated,features = "Spermatogonia1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Diff.Spg1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Spermatocytes1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Round1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Elongating1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Leydig1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Sertoli1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Myoid1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Endothelial1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Macrophages1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Tcells1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(mouse.multiome.integrated,features = "Bcells1",reduction = "umap.atac", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme()
)

FeaturePlot(mouse.multiome.integrated, features = c("Dcn", "Loxl1", "Col1a1", "Acta2"), reduction = "umap.atac", min.cutoff = 0, order = F)

DefaultAssay(mouse.multiome.integrated) <- "RNA"
Idents(mouse.multiome.integrated) <- mouse.multiome.integrated[["seurat_clusters"]]
DimPlot(mouse.multiome.integrated, label = T, reduction = "umap.atac")
FeaturePlot(mouse.multiome.integrated, features = "nCount_RNA", max.cutoff = 10000, reduction = "umap.atac")

FeaturePlot(mouse.multiome.integrated, features = c("Dcn", "Loxl1", "Col1a1", "Acta2"), reduction = "umap.atac", min.cutoff = 0, order = F)

DefaultAssay(mouse.multiome.integrated) <- "RNA"
Idents(mouse.multiome.integrated) <- mouse.multiome.integrated[["seurat_clusters"]]
DimPlot(mouse.multiome.integrated, label = T, reduction = "umap.atac")
# FeaturePlot(mouse.multiome.integrated, features = "nCount_RNA", max.cutoff = 10000, reduction = "umap.atac")

new.cluster.ids <- character(length(levels(Idents(mouse.multiome.integrated))))

spermatogonia.ids <- c(31, 13, 38)+1
diff.spermatogonia.ids <- c(14, 19)+1
meiotic.ids <- c(27, 23,21, 25)+1
spermatid.ids <- c(17, 28, 8, 16)+1
telocytes.ids <-c(0, 11)+1
myoidintermediate.ids <- c(6,3,22,2)+1
Leydig.ids <- c(9, 10, 32, 1, 4)+1
Sertoli.ids <-c(24, 26)+1
myoid.ids <- c(35, 7, 30, 18)+1
# lepzyg.ids <- c(27, 20, 24)
Tcells.ids <-c(12)+1
macrophages.ids <-c(5, 20, 33, 29)+1
endotheliala.ids <- c(37, 15, 34, 36)+1

new.cluster.ids[spermatogonia.ids] <- "spermatogonia"
new.cluster.ids[diff.spermatogonia.ids] <- "diff. spermatogonia"
new.cluster.ids[meiotic.ids] <- "meiotic"
new.cluster.ids[spermatid.ids] <- "spermatids"
new.cluster.ids[telocytes.ids] <- "telocytes"
new.cluster.ids[myoidintermediate.ids] <- "myoid intermediate"
new.cluster.ids[Leydig.ids] <- "Leydig"
new.cluster.ids[Sertoli.ids] <- "Sertoli"
new.cluster.ids[myoid.ids] <- "myoid"
# new.cluster.ids[lepzyg.ids] <- "Leptotene/Zygotene"
new.cluster.ids[Tcells.ids] <- "T cells"
new.cluster.ids[macrophages.ids] <- "macrophages"
new.cluster.ids[endotheliala.ids] <- "endothelial"



new.cluster.ids

names(x = new.cluster.ids) <- levels(x = mouse.multiome.integrated)
mouse.multiome.integrated <- RenameIdents(object = mouse.multiome.integrated, new.cluster.ids)
mouse.multiome.integrated$cell.type <- Idents(mouse.multiome.integrated)

mouse.multiome.integrated$cell.type <- factor(x = mouse.multiome.integrated$cell.type, levels = c("spermatogonia",
                                                                                                  "diff. spermatogonia",
                                                                                                  "meiotic",
                                                                                                  "spermatids",
                                                                                                  "telocytes",
                                                                                                  "myoid intermediate",
                                                                                                  "myoid",
                                                                                                  "Leydig",
                                                                                                  "Sertoli",
                                                                                                  "T cells",
                                                                                                  "macrophages",
                                                                                                  "endothelial"
))

DimPlot(mouse.multiome.integrated, reduction = "umap.atac", label = F)





FeaturePlot(mouse.multiome.integrated, features = c("Amhr2", "Sox9", "Cyp17a1"), 
            min.cutoff = 0, order = T, 
            # split.by = "treatment",
            reduction = "umap.wnn")


FeaturePlot(mouse.multiome.integrated, features = "nCount_RNA", reduction = "umap.rna",max.cutoff = 10000)

DimPlot(mouse.multiome.integrated, split.by = "treatment", group.by = "cell.type")



gene.of.interest <- c("Sox9")
CoveragePlot(mouse.multiome.integrated,
             extend.downstream = 5000,
             extend.upstream = 5000,
             region = "chr6-25437229-25471723", 
             features = "Cd9", 
             assay = 'peaks', 
             group.by = "cell.type",
             expression.assay = 'RNA', peaks = F)


fragments.list <- (Fragments(mouse.multiome.germ))
head(fragments.list[[1]])
head(Annotation(mouse.multiome.germ))





seqlevelsStyle(BSgenome.Mmusculus.UCSC.mm10) <- "UCSC"  #"NCBI"


DefaultAssay(mouse.multiome.integrated) <- "peaks"
mouse.multiome.integrated <- RegionStats(mouse.multiome.integrated, genome = BSgenome.Mmusculus.UCSC.mm10, verbose = T)


mouse.multiome.integrated <- LinkPeaks(
  object = mouse.multiome.integrated,
  peak.assay = "peaks",
  # peak.slot = ""
  expression.assay = "RNA",
  genes.use = c("Tex12")
)

gene.of.interest <- "Mael"

CoveragePlot(mouse.multiome.integrated,
             extend.downstream = 5000,
             extend.upstream = 5000,
             region = gene.of.interest, 
             features = gene.of.interest, 
             assay = 'peaks', 
             group.by = "cell.type",
             expression.assay = 'RNA', peaks = F)


sertoli.subset <- subset(mouse.multiome.integrated, idents = "Sertoli")
leydig.subset <- subset(mouse.multiome.integrated, idents = "Leydig")

gene.of.interest <- "Serpina5"

CoveragePlot(leydig.subset,
             extend.downstream = 5000,
             extend.upstream = 5000,
             region = gene.of.interest, 
             features = gene.of.interest, 
             assay = 'peaks', 
             group.by = "treatment",
             expression.assay = 'RNA', peaks = F)


DefaultAssay(sertoli.subset) <- "peaks"
Idents(sertoli.subset) <- sertoli.subset[["treatment"]]
da_peaks <- FindMarkers(
  object = sertoli.subset,
  ident.1 = "Control",
  ident.2 = "GFPtransplanted",
  only.pos = FALSE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

sertoli.subset <- LinkPeaks(
  object = sertoli.subset,
  peak.assay = "peaks",
  # peak.slot = ""
  expression.assay = "RNA",
  genes.use = c("Ext1")
)

CoveragePlot(sertoli.subset,
             extend.downstream = 130000,
             extend.upstream = 10000,
             region = "chr15-53239238-53240312", #top.da.peak[5], 
             features = "Ext1",
             assay = 'peaks', 
             group.by = "treatment",
             expression.assay = 'RNA', peaks = F)





Idents(mouse.multiome.integrated) <- mouse.multiome.integrated[["cell.type"]]
DimPlot(mouse.multiome.integrated, reduction = "umap.wnn", label = T)
DefaultAssay(mouse.multiome.integrated) <- "RNA"
FeaturePlot(mouse.multiome.integrated, reduction = "umap.wnn", features = "Stra8", order = T)

DefaultAssay(mouse.multiome.integrated) <- "peaks"
# pfm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(collection = c("CORE"), tax_group = 'vertebrates', all_versions = FALSE)
# )

pfm2 <- readJASPARMatrix("/venice/brinsterlab/output_10X/Jaspar_nonredundant_curated.txt", 
                         matrixClass="PFM")

mouse.multiome.integrated <- AddMotifs(
  object = mouse.multiome.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm2 
)

mouse.multiome.integrated <- RunChromVAR(
  object = mouse.multiome.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10
)




top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])


enriched.motifs <- FindMotifs(
  object = multiome.integrated,
  features = top.da.peak
)

enriched.motifs$gene.name <- str_to_title(enriched.motifs$motif.name)

enriched.motifs



# MA1603.1 = Dmrt1

DefaultAssay(mouse.multiome.integrated) <- "peaks"
mouse.multiome.integrated <- Footprint(
  object = mouse.multiome.integrated,
  motif.name = c("Dmrt1"), #this will work with the TF gene name but it is much faster with the motif name 
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# plot the footprint data for each group of cells
PlotFootprint(mouse.multiome.integrated, features = c("Dmrt1"))


gene_of_interest <- "Dmrt1"

DefaultAssay(mouse.multiome.integrated) <- 'chromvar'
# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = mouse.multiome.integrated,
  features = gene.of.interest,
  # order = T,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  reduction = "umap",
  pt.size = 0.1
)









#save files
setwd("/venice/brinsterlab/seurat.object.saves")
DimPlot(mouse.multiome.germ)
# save(mouse.multiome.germ, file = "mouse.multiome.germ.2023.03.27.Robj")


DefaultAssay(mouse.multiome.germ) <- "ATAC"
GeneActivityMatrix <- GeneActivity(mouse.multiome.germ)

mouse.multiome.germ[['promoter']] <- CreateAssayObject(counts = GeneActivityMatrix)
mouse.multiome.germ <- NormalizeData(
  object = mouse.multiome.germ,
  assay = 'promoter',
  normalization.method = 'LogNormalize',
  scale.factor = median(mouse.multiome.germ$nCount_promoter)
)
DefaultAssay(mouse.multiome.germ) <- "promoter"
FeaturePlot(mouse.multiome.germ, reduction = "umap.wnn", features = "Sdc4", order = T) + scale_color_viridis()

# save(mouse.multiome.germ, file = "mouse.multiome.germ.2023.03.28.Robj")

somatic.markers.mouse <- FindMarkers(somatic.integrated, ident.1 = c(17, 15))
write.csv(somatic.markers.mouse, file="somatic.markers.mouse.csv")



### ANALYSIS FOR PAPER ###






setwd("/Users/ewhelan/Desktop/Somatic Project/somatic_samples")
load("somatic.integrated.germ4.Robj")
# load("somatic.integrated.somatic4.Robj")
setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")
load("somatic.integrated.somatic6.Robj")


#calculate the average values
nCount_values.germ <- somatic.integrated.germ4$nCount_RNA
nCount_values.somatic <- somatic.integrated.somatic6$nCount_RNA
nCount_values <- c(nCount_values.germ, nCount_values.somatic)
average_nCount_RNA <- mean(nCount_values, na.rm = TRUE)
average_nCount_RNA

nFeature_values.germ <- somatic.integrated.germ4$nFeature_RNA
nFeature_values.somatic <- somatic.integrated.somatic6$nFeature_RNA
nFeature_values <- c(nFeature_values.germ, nFeature_values.somatic)
average_nFeature_RNA <- mean(nFeature_values, na.rm = TRUE)
average_nFeature_RNA

percent.mito_values.germ <- somatic.integrated.germ4$percent.mito
percent.mito_values.somatic <- somatic.integrated.somatic6$percent.mito
percent.mito_values <- c(percent.mito_values.germ, percent.mito_values.somatic)
average_percent.mito <- mean(percent.mito_values, na.rm = TRUE)
average_percent.mito





FeaturePlot(somatic.integrated.somatic6, features = c("Tcf21", "Loxl1", "Acta2", "Cyp17a1",
                                                      "Dcn", "Col1a1", "Star", "Hsd3b1"), ncol = 4)  & DarkTheme() & scale_color_viridis(option = "C")

Idents(somatic.integrated.somatic6) <- "cell.type3"

FeaturePlot(somatic.integrated.germ4, order = T,features = c("Agt", "Tnf", "Il6", "Ifng", "Tgfb1", "Egf",
                                                             "Hgf", "Vegfa", "Igf1", "Fgf2", "Tgfb3", "Fgf10",
                                                             "Nrg1", "Angpt2"), ncol = 4)  & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(somatic.integrated.somatic6, order = T,features = c("Foxb2", "Cck", "Tspan15"), split.by = "treatment")  & DarkTheme() & scale_color_viridis(option = "C")


col.vector <- c("#b12625",
                "#566b30",
                "#d9a528",
                "#cd8163",
                "#88cdea",
                "#7f8080",
                "#d4a1ca",
                "#7c57a4",
                "#6ec6a8",
                "#9bcb3c",
                "#ee6463",
                "#5570b6",
                "#f06ca8",
                "#f78c1e"
)


col.vector2 <- c( "chocolate4",
                  "cornsilk4",
                  "darkseagreen3",
                  "plum4",
                  "paleturquoise",
                  "darkgreen",
                  "blue4",
                  "gold2"
)

DimPlot(somatic.integrated.germ4, cols = col.vector2, split.by = "treatment")
DimPlot(somatic.integrated.germ4, cols = col.vector2)
DimPlot(somatic.integrated.somatic6, split.by = "treatment", cols = col.vector, label = F)
DimPlot(somatic.integrated.somatic6,  cols = col.vector, label = T)

DimPlot(somatic.integrated.somatic6, group.by = "seurat_clusters")



FeaturePlot(somatic.integrated.germ4, features = c("Adam3", "Ropn1l", "Tnp1", "Odf1"))
VlnPlot(somatic.integrated.somatic6, features = c("Anp32b", "Ccl5", "C1qb", "Abcb1a", "Tpm1", "Col1a2", "Fabp3"), pt.size = 0)
VlnPlot(somatic.integrated.somatic6, features = c("Bex1", "Mxra5", "Dcx", "Vcxc"), pt.size = 0)




somatic.integrated.somatic6$cell.type3 <- factor(somatic.integrated.somatic6$cell.type3, levels = c(
  "Sertoli", 
  "Leydig_cells",
  "Peritubular_myoid",
  "Perivascular_smooth_muscle",
  "Mesenchymal_progenitors",
  "Epithelial_cells",
  "Endothelial", 
  "Macrophages_peritubular",
  "Macrophages_interstitial",
  "Dendritic_cells",
  "T_cells_immunomodulatory",
  "T_cells_effector",
  "B_cells"
))
Idents(somatic.integrated.somatic6) <- somatic.integrated.somatic6$cell.type3
DimPlot(somatic.integrated.somatic6, cols = col.vector, label = F, split.by = "treatment")

setwd("~/Desktop/Somatic Project/Manuscript/Figures/InitialFigures")

# Open a PDF device
pdf("cell_type_rasterized_plots.pdf", width = 12, height = 8)

# Use tryCatch to ensure dev.off() is executed properly
tryCatch({
  # Create an empty base plot
  base_plot <- DimPlot(somatic.integrated.somatic6, cols = col.vector, label = FALSE) +
    ggtitle("All Cell Types Rasterized Separately") +
    theme_minimal()  # Optional: Adjust theme as needed
  
  # Rasterize each cell type separately
  for (cell_type in unique(somatic.integrated.somatic6$cell.type3)) {
    # Subset data for the current cell type
    subset_data <- subset(somatic.integrated.somatic6, cell.type3 == cell_type)
    
    # Add the rasterized layer for the current cell type
    base_plot <- base_plot +
      geom_point_rast(data = subset_data@meta.data, aes(x = subset_data@reductions$umap@cell.embeddings[, 1],
                                                        y = subset_data@reductions$umap@cell.embeddings[, 2]),
                      alpha = 0.8, size = 0.5)  # Adjust alpha and size as needed
  }
  
  # Print the plot into the PDF
  print(base_plot)
}, error = function(e) {
  # Catch any errors and print a message
  message("An error occurred: ", e$message)
}, finally = {
  # Ensure the PDF device is properly closed
  dev.off()
})







Idents(somatic.integrated.somatic6) <- "cell.type3"
DimPlot(somatic.integrated.somatic6)

myoid.cells <- subset(somatic.integrated.somatic6, idents = "Peritubular_myoid")
Idents(myoid.cells) <- "treatment"

myoid.markers <- FindMarkers(myoid.cells, ident.1 = "control", ident.2 = "GFP transplanted")
write.csv(myoid.markers, file="myoid.markers.mouse.csv")


sertoli.cells <- subset(somatic.integrated.somatic6, idents = "Sertoli")
Idents(sertoli.cells) <- "treatment"

sertoli.markers <- FindMarkers(sertoli.cells, ident.1 = "control", ident.2 = "GFP transplanted")
write.csv(sertoli.markers, file="sertoli.markers.mouse.csv")
DimPlot(somatic.integrated.somatic6, label = T)
FeaturePlot(somatic.integrated.somatic6, features = c("Kcnma1", "Socs3"), order = T,split.by = "treatment", min.cutoff = 0)

VlnPlot(somatic.integrated.somatic6, features = c("Kcnma1", "Socs3"), split.by = "treatment", pt.size = 0)

#load human data
load("/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/human.multiome.adult.prepubertal.01.11.2025.Robj")
DimPlot(human.multiome.integrated2, reduction = "umap.wnn")


human.multiome.integrated2 <- UpdatePath(human.multiome.integrated2, "/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/fragpath/H257B_multi/atac_fragments.tsv.gz")



DefaultAssay(human.multiome.integrated2) <- "peaks"
frags <- Fragments(human.multiome.integrated2)  # get list of fragment objects
Fragments(human.multiome.integrated2) <- NULL  # remove fragment information from assay

# create a vector with all the new paths, in the correct order for your list of fragment objects
# In this case we only have 1




new.paths <- list(
  "/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/fragpath/HS99/atac_fragments.tsv.gz",
  "/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/fragpath/HS102/atac_fragments.tsv.gz",
  "/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/fragpath/HS103/atac_fragments.tsv.gz",
  "/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/fragpath/HS104/atac_fragments.tsv.gz",
  "/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/fragpath/HS107/atac_fragments.tsv.gz",
  "/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/fragpath/HS112/atac_fragments.tsv.gz"
)

for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}

Fragments(human.multiome.integrated2) <- frags # assign updated list back to the object
Fragments(human.multiome.integrated2)


CoveragePlot(human.multiome.integrated2,
             extend.downstream = 5000,
             extend.upstream = 5000,
             ymax = 100,
             region = "SOX9", 
             features = "SOX9", 
             assay = 'peaks', 
             group.by = "cell.type",split.by = "stage",
             expression.assay = 'RNA', peaks = F)


DefaultAssay(human.multiome.integrated2) <- "RNA"
Idents(human.multiome.integrated2) <- "cell.type"
all.markers <- FindAllMarkers(human.multiome.integrated2)


all.markers.by.treatment <- FindAllMarkers(human.multiome.integrated2, group.by = "stage")


spermatid_markers <- FindMarkers(human.multiome.integrated2, ident.1 = "Spermatids", min.pct = 0.2 )

filtered_genes <- spermatid_markers %>%
  filter(avg_log2FC > 1, p_val < 0.05)

# View the filtered results
blacklist <- rownames(filtered_genes)


Idents(human.multiome.integrated2) <- "cell.type"
DimPlot(human.multiome.integrated2)

get_counts <- function(filtered_data, cell_type) {
  counts <- table(ifelse(filtered_data$avg_log2FC > 0, "positive", "negative"))
  return(data.frame(CellType = cell_type, Positive = counts["positive"], Negative = counts["negative"]))
}
result_list <- list()

setwd("~/Desktop/multiomeHuman/DEGs_filtered")

cell_types <- c("Mesenchymal_progenitors", "Sertoli_cells", "Endothelial", "Peritubular_myoid", "Macrophages", "T_cells", "Leydig_cells")
for (cell_type in cell_types) {
  # cell_type <- "Leydig_cells"
  print(cell_type)
  cells <- subset(human.multiome.integrated2, idents = cell_type)
  Idents(cells) <- "stage"
  adult_vs_prepubertal <- FindMarkers(cells, ident.1 = "adult")
  filtered <- adult_vs_prepubertal[!rownames(adult_vs_prepubertal) %in% blacklist, ]
  filtered <- filtered %>% filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05)
  write.csv(filtered, file=paste0(cell_type, ".csv"))
  # Get the counts for positive and negative log2FC
  result_list[[cell_type]] <- get_counts(filtered, cell_type)
}

# Combine all the results into a single data frame
final_results <- do.call(rbind, result_list)

# View the final results
print(final_results)


FeaturePlot(human.multiome.integrated2, split.by = "stage", features = "ZFPM2", reduction = "umap.wnn") & DarkTheme() & scale_color_viridis(option = "C")





#venn diagrams



load("/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/humanSomaticMultiome/SCENICPlus_Output/eGRN_AUC.h5seurat")
hfile <- LoadH5Seurat("/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/humanSomaticMultiome/SCENICPlus_Output/eGRN_AUC.h5seurat")


VlnPlot(human.multiome.integrated2, features = "PSAP", group.by = "cell.type", split.by = "stage")


load("/Users/ewhelan/Documents/scRNAseq/humanSomaticMultiome/human.multiome.adult.prepubertal.01.11.2025.Robj")

Idents(somatic.integrated.somatic6) <- "cell.type3"
DimPlot(somatic.integrated.somatic6)
sertoli.mouse <- subset(somatic.integrated.somatic6, idents = "Sertoli")


Idents(human.multiome.integrated2) <- "cell.type"
DimPlot(human.multiome.integrated2, reduction = "umap.wnn")
sertoli.human <- subset(human.multiome.integrated2, idents = "Sertoli_cells")


Idents(sertoli.mouse) <- "treatment"
DefaultAssay(sertoli.mouse) <- "RNA"

Idents(sertoli.human) <- "stage"
DefaultAssay(sertoli.human) <- "RNA"

gene.of.interest <- c("Psap", "Hbegf", "Inha", "Lgals1")
gene.of.interest <- c("Notch3", "Antxr1", "Col27a1", "Vcl")
gene.of.interest <- c("Notch3", "Antxr1", "Col27a1", "Vcl")
gene.of.interest1 <- c("Psap", "Dhh", "Lama1", "Cd63", #"Lpar3", "F3", "Amfr", "Ptprd", 
                       "Arf6", "Afdn", "Ripk1", "Bambi", "Adam10", #"Plxnb2", "Bsg",
                       "Myl9","Canx")
gene.of.interest2 <- c("Timp3", "Calm3", "Nucb2", "Glg1", "Lamc1", #"Fzd3", 
                       "Hspa4", #"Lman1", 
                       "Lamp1",  "Smap1", "Mylk", "Mcfd2", 
                       # "Lrpap1",
                       "Vasn")
# gene.of.interest2 <- c("Psap", "Dhh")

p1 <- Stacked_VlnPlot(sertoli.mouse, features = gene.of.interest1, pt.size = 0, ncol = 1)
p2 <- Stacked_VlnPlot(sertoli.human, features = toupper(gene.of.interest1), pt.size = 0, ncol = 1)
p3 <- Stacked_VlnPlot(sertoli.mouse, features = gene.of.interest2, pt.size = 0, ncol = 1)
p4 <- Stacked_VlnPlot(sertoli.human, features = toupper(gene.of.interest2), pt.size = 0, ncol = 1)

plot_grid(p1, p2, p3, p4, ncol = 4)



in.common.genes <- read.csv("/Users/ewhelan/Desktop/humanSomaticGRN/in_common_factors.csv")


gene.of.interest <- c("Psap", "Dhh")
# gene.of.interest2 <- c("Psap", "Dhh")

p1 <- Stacked_VlnPlot(sertoli.mouse, features = gene.of.interest, pt.size = 0, ncol = 1)
p2 <- Stacked_VlnPlot(sertoli.human, features = toupper(gene.of.interest), pt.size = 0, ncol = 1)

plot_grid(p1, p2)

for(current.loop in 1:46){
  # current.loop <- 1
  gene.of.interest.mouse <- in.common.genes[(current.loop*4-3):(current.loop*4),2]
  gene.of.interest.human <- in.common.genes[(current.loop*4-3):(current.loop*4),1]
  
  p1 <- Stacked_VlnPlot(sertoli.mouse, features = gene.of.interest.mouse, pt.size = 0, ncol = 1)
  p2 <- Stacked_VlnPlot(sertoli.human, features = gene.of.interest.human, pt.size = 0, ncol = 1)
  
  print(plot_grid(p1, p2))
}



regulon_genes <- read_csv("/Users/ewhelan/Desktop/humanSomaticAnalysis/regulon_genes_mouse_human.csv")



Idents(human.multiome.integrated2) <- "cell.type"
Idents(somatic.integrated.somatic6) <- "cell.type"
DefaultAssay(somatic.integrated.somatic6) <- "RNA"
avgexp_human = AverageExpression(human.multiome.integrated2, return.seurat = T, add.ident = 'stage')
avgexp_mouse = AverageExpression(somatic.integrated.somatic6, return.seurat = T, add.ident = 'treatment')

DoHeatmap(avgexp_mouse, features = (unique(c(regulon_genes$Mouse_Sertoli_downstream, str_to_title(regulon_genes$Human_Sertoli_downstream)))))
DoHeatmap(avgexp_human, features = (unique(c(regulon_genes$Mouse_Sertoli_downstream_UC, regulon_genes$Human_Sertoli_downstream))))

table(Idents(human.multiome.integrated2))
human.sertoli <- subset(human.multiome.integrated2, idents = "Sertoli_cells")
table(Idents(somatic.integrated.somatic6))
mouse.sertoli <- subset(somatic.integrated.somatic6, idents = "Sertoli")

Idents(human.sertoli) <- "stage"
table(Idents(human.sertoli))
Idents(mouse.sertoli) <- "treatment"
table(Idents(mouse.sertoli))

data.sertoli.human.prepubertal <- as.matrix(GetAssayData(human.sertoli, slot = "data")[, WhichCells(human.sertoli, ident = "prepubertal")])
data.sertoli.human.adult <- as.matrix(GetAssayData(human.sertoli, slot = "data")[, WhichCells(human.sertoli, ident = "adult")])

data.sertoli.mouse.control <- as.matrix(GetAssayData(mouse.sertoli, slot = "data")[, WhichCells(mouse.sertoli, ident = "control")])
data.sertoli.mouse.GFP <- as.matrix(GetAssayData(mouse.sertoli, slot = "data")[, WhichCells(mouse.sertoli, ident = "GFP transplanted")])

rownames(data.sertoli.mouse.control) <- toupper(rownames(data.sertoli.mouse.control))
rownames(data.sertoli.mouse.GFP) <- toupper(rownames(data.sertoli.mouse.GFP))
length(rownames(data.sertoli.mouse.control))
common.gene.names <- (intersect(rownames(data.sertoli.mouse.control), rownames(data.sertoli.human.prepubertal)))

data.sertoli.mouse.control.intersect <- data.sertoli.mouse.control[common.gene.names,]
data.sertoli.mouse.GFP.intersect <- data.sertoli.mouse.GFP[common.gene.names,]
data.sertoli.human.prepubertal.intersect <- data.sertoli.human.prepubertal[common.gene.names,]
data.sertoli.human.adult.intersect <- data.sertoli.human.adult[common.gene.names,]

av.counts.sertoli.mouse.control.intersect <- apply(data.sertoli.mouse.control.intersect, 1, mean)
av.counts.sertoli.mouse.GFP.intersect <- apply(data.sertoli.mouse.GFP.intersect, 1, mean)
av.counts.sertoli.human.prepubertal.intersect <- apply(data.sertoli.human.prepubertal.intersect, 1, mean)
av.counts.sertoli.human.adult.intersect <- apply(data.sertoli.human.adult.intersect, 1, mean)


av.counts.all <- cbind(av.counts.sertoli.mouse.control.intersect,
                       av.counts.sertoli.mouse.GFP.intersect,
                       av.counts.sertoli.human.prepubertal.intersect,
                       av.counts.sertoli.human.adult.intersect
)

dim(av.counts.all)

av.counts.all.df <- as.data.frame(av.counts.all)


gene.name.list.all <- rownames(av.counts.all.df)



###
hs_column <- av.counts.all.df  
ncg_column <- av.counts.all.df 

correlation_values <- cor(hs_column, ncg_column, use = "pairwise.complete.obs")


print(correlation_values)
rsquared <- correlation_values^2
print(rsquared)

AdultData <- av.counts.all.df
PrenatalData <- av.counts.all.df
correlation_matrix <- cor(AdultData, PrenatalData)
print(correlation_matrix)
correltation_matrix_r <- correlation_matrix
correlation_matrix <- correlation_matrix^2



# Convert the correlation matrix to a tibble
correlation_df <- as_tibble(as.data.frame(rsquared), rownames = "Cell_Type_1")

# Factor the levels to maintain the original order
correlation_df$Cell_Type_1 <- factor(correlation_df$Cell_Type_1, levels = rownames(correlation_matrix))
correlation_df <- gather(correlation_df, key = "Cell_Type_2", value = "Correlation", -Cell_Type_1)

# Factor the levels for Cell_Type_2
correlation_df$Cell_Type_2 <- factor(correlation_df$Cell_Type_2, levels = colnames(correlation_matrix))

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(correlation_df, aes(x = Cell_Type_1, y = Cell_Type_2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +  # Add correlation values
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  # scale_fill_gradient2(low = "#451d5b", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title = "Correlation Heatmap",
       x = "Cell Type (In Vivo)",
       y = "Cell Type (In Vitro)") +
  theme(axis.text.x = element_text(angle = -22, hjust = 0))  # Rotate x-axis labels at a 45-degree angle



print(heatmap_plot)

# regulon_genes_all <- unique(c(regulon_genes$Human_Sertoli_TFs,
#                        regulon_genes$Mouse_Sertoli_TFs_UC))

regulon_genes_all <- unique(c(regulon_genes$Human_Sertoli_downstream,
                              regulon_genes$Mouse_Sertoli_downstream_UC))

genes.of.interest <- intersect(regulon_genes_all, common.gene.names)


av.counts.sertoli.mouse.control.intersect.df <- as.data.frame(av.counts.sertoli.mouse.control.intersect)
av.counts.sertoli.mouse.GFP.intersect.df <- as.data.frame(av.counts.sertoli.mouse.GFP.intersect)
av.counts.sertoli.human.prepubertal.intersect.df <- as.data.frame(av.counts.sertoli.human.prepubertal.intersect)
av.counts.sertoli.human.adult.intersect.df <- as.data.frame(av.counts.sertoli.human.adult.intersect)

gene.name.list.ncg <- rownames(av.counts.sertoli.mouse.control.intersect.df)


av.counts.sertoli.mouse.control.intersect.df$gene.name <- gene.name.list.ncg
av.counts.sertoli.mouse.GFP.intersect.df$gene.name <- gene.name.list.ncg
av.counts.sertoli.human.prepubertal.intersect.df$gene.name <- gene.name.list.ncg
av.counts.sertoli.human.adult.intersect.df$gene.name <- gene.name.list.ncg

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff


av.counts.sertoli.mouse.control.intersect.shortlist <- av.counts.sertoli.mouse.control.intersect.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)

av.counts.sertoli.mouse.GFP.intersect.shortlist <- av.counts.sertoli.mouse.GFP.intersect.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)

av.counts.sertoli.human.prepubertal.intersect.shortlist <- av.counts.sertoli.human.prepubertal.intersect.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)

av.counts.sertoli.human.adult.intersect.shortlist <- av.counts.sertoli.human.adult.intersect.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)


heatmap.matrix <- cbind(av.counts.sertoli.mouse.control.intersect.shortlist$av.counts.sertoli.mouse.control.intersect,
                        av.counts.sertoli.mouse.GFP.intersect.shortlist$av.counts.sertoli.mouse.GFP.intersect,
                        av.counts.sertoli.human.prepubertal.intersect.shortlist$av.counts.sertoli.human.prepubertal.intersect,
                        av.counts.sertoli.human.adult.intersect.shortlist$av.counts.sertoli.human.adult.intersect)

rownames(heatmap.matrix) <- rownames(av.counts.sertoli.mouse.control.intersect.shortlist)
colnames(heatmap.matrix) <- c("mouse.control", "mouse.GFP", "human.prepubertal", "human.adult")



# Load necessary library


# Scale the data by row (genes)
scaled_matrix <- t(scale(t(heatmap.matrix)))

# Generate the heatmap
pheatmap(scaled_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "none",  # Already scaled manually
         color = colorRampPalette(c("blue", "white", "red"))(50))  # Color gradient


# Compute log2 fold change (log2FC) while avoiding division by zero
log2FC_human <- log2((heatmap.matrix[, "human.adult"] + 1e-6) / (heatmap.matrix[, "human.prepubertal"] + 1e-6))
log2FC_mouse <- log2((heatmap.matrix[, "mouse.control"] + 1e-6) / (heatmap.matrix[, "mouse.GFP"] + 1e-6))

# Create a matrix for heatmap
log2FC_matrix <- as.matrix(data.frame(Human_log2FC = log2FC_human, Mouse_log2FC = log2FC_mouse))
rownames(log2FC_matrix) <- rownames(heatmap.matrix)  # Keep gene names

# Replace NA, NaN, or Inf with zero
log2FC_matrix[is.na(log2FC_matrix) | is.infinite(log2FC_matrix)] <- 0

# Generate the heatmap
log2FC_matrix.capped <- pmax(pmin(log2FC_matrix, 5), -5)

# Generate the heatmap
pheatmap(log2FC_matrix.capped,
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # No need to cluster conditions
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scale
         main = "Log2 Fold Change Heatmap")


# Define conditions
same_sign <- (log2FC_matrix[, "Human_log2FC"] * log2FC_matrix[, "Mouse_log2FC"]) > 0
different_sign <- (log2FC_matrix[, "Human_log2FC"] * log2FC_matrix[, "Mouse_log2FC"]) < 0
similar_values <- (abs(log2FC_matrix[, "Human_log2FC"]) < 1) & (abs(log2FC_matrix[, "Mouse_log2FC"]) < 1)

# Exclude similar values from same_sign and different_sign calculations
same_sign <- same_sign & !similar_values
different_sign <- different_sign & !similar_values

# Compute percentages
percent_same_sign <- mean(same_sign, na.rm = TRUE) * 100
percent_different_sign <- mean(different_sign, na.rm = TRUE) * 100
percent_similar <- mean(similar_values, na.rm = TRUE) * 100

# Print results
cat("Percentage with the same sign (excluding similar values):", round(percent_same_sign, 2), "%\n")
cat("Percentage with different signs (excluding similar values):", round(percent_different_sign, 2), "%\n")
cat("Percentage with similar values (< ±1):", round(percent_similar, 2), "%\n")

### try again with just mouse and make a rounded dotplot


data.sertoli.mouse.control <- as.matrix(GetAssayData(mouse.sertoli, slot = "data")[, WhichCells(mouse.sertoli, ident = "control")])
data.sertoli.mouse.GFP <- as.matrix(GetAssayData(mouse.sertoli, slot = "data")[, WhichCells(mouse.sertoli, ident = "GFP transplanted")])


av.counts.sertoli.mouse.control <- apply(data.sertoli.mouse.control, 1, mean)
av.counts.sertoli.mouse.GFP <- apply(data.sertoli.mouse.GFP, 1, mean)



av.counts.all <- cbind(av.counts.sertoli.mouse.control,
                       av.counts.sertoli.mouse.GFP
)

dim(av.counts.all)

av.counts.all.df <- as.data.frame(av.counts.all)



# Convert data to long format
av.counts.all.df$gene <- rownames(av.counts.all.df)  # Ensure gene column exists
plot_data <- av.counts.all.df %>%
  filter(gene %in% gene.list) %>%
  pivot_longer(cols = -gene, names_to = "Treatment", values_to = "Expression")

# Explicitly map each treatment to **fixed radii**
plot_data <- plot_data %>%
  mutate(Ring = case_when(
    Treatment == "av.counts.sertoli.mouse.control" ~ 5,  # Control at radius 5
    Treatment == "av.counts.sertoli.mouse.GFP" ~ 10      # GFP at radius 10
  ))

# Ensure no unwanted NA or zero values in Ring
plot_data <- plot_data %>%
  filter(!is.na(Ring) & Ring > 0)

# Create circular dot plot with fixed radial positions
ggplot(plot_data, aes(x = gene, y = Ring, size = Expression, fill = Expression)) +
  geom_point(shape = 21, color = "black") +  # Black outline for clarity
  scale_size(range = c(2, 8)) +  # Adjust dot size range
  scale_fill_gradient(low = "blue", high = "red") +  # Expression intensity
  coord_radial(theta = "x", clip = "off") +  # Proper radial layout
  theme_minimal() +
  theme(
    axis.title = element_blank(), 
    axis.text.y = element_blank(),  # Hide radial labels
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", vjust = 1)
  ) +
  labs(title = "Circular Dot Plot of Gene Expression", size = "Expression Level")


### ### ### ### ### ### #
### ### ### ### ### ### #
### redo the dotplots ###
### ### ### ### ### ### #
### ### ### ### ### ### #
### ### ### ### ### ### #
### ### ### ### ### ### #
### ### ### ### ### ### #
### ### ### ### ### ### #
### ### ### ### ### ### #


setwd("/Users/ewhelan/Desktop/Somatic Project/somatic_samples")
load("somatic.integrated.germ4.Robj")
DimPlot(somatic.integrated.germ4)

data.SSCs <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "SSCs")])
data.Progenitors <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "Progenitors")])
data.DiffSpermatogonia <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "DiffSpermatogonia")])
data.PrelepSpermatocytes <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "PrelepSpermatocytes")])
data.EarlySpermatocytes <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "EarlySpermatocytes")])
data.LateSpermatocytes <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "LateSpermatocytes")])
data.RoundSpermatids <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "RoundSpermatids")])
data.ElongatingSpermatids <- as.matrix(GetAssayData(somatic.integrated.germ4, slot = "data")[, WhichCells(somatic.integrated.germ4, ident = "ElongatingSpermatids")])





av.counts.SSCs <- apply(data.SSCs, 1, mean)
av.counts.Progenitors <- apply(data.Progenitors, 1, mean)
av.counts.DiffSpermatogonia <- apply(data.DiffSpermatogonia, 1, mean)
av.counts.PrelepSpermatocytes <- apply(data.PrelepSpermatocytes, 1, mean)
av.counts.EarlySpermatocytes <- apply(data.EarlySpermatocytes, 1, mean)
av.counts.LateSpermatocytes <- apply(data.LateSpermatocytes, 1, mean)
av.counts.RoundSpermatids <- apply(data.RoundSpermatids, 1, mean)
av.counts.ElongatingSpermatids <- apply(data.ElongatingSpermatids, 1, mean)


av.counts.all <- cbind(av.counts.SSCs,
                       av.counts.Progenitors,
                       av.counts.DiffSpermatogonia,
                       av.counts.PrelepSpermatocytes,
                       av.counts.EarlySpermatocytes,
                       av.counts.LateSpermatocytes,
                       av.counts.RoundSpermatids,
                       av.counts.ElongatingSpermatids
)

dim(av.counts.all)

av.counts.all.df <- as.data.frame(av.counts.all)


gene.name.list.all <- rownames(av.counts.all.df)



###
hs_column <- av.counts.all.df  
ncg_column <- av.counts.all.df 

correlation_values <- cor(hs_column, ncg_column, use = "pairwise.complete.obs")


print(correlation_values)
rsquared <- correlation_values^2
print(rsquared)

AdultData <- av.counts.all.df
PrenatalData <- av.counts.all.df
correlation_matrix <- cor(AdultData, PrenatalData)
print(correlation_matrix)
correltation_matrix_r <- correlation_matrix
correlation_matrix <- correlation_matrix^2



# Convert the correlation matrix to a tibble
correlation_df <- as_tibble(as.data.frame(rsquared), rownames = "Cell_Type_1")

# Factor the levels to maintain the original order
correlation_df$Cell_Type_1 <- factor(correlation_df$Cell_Type_1, levels = rownames(correlation_matrix))
correlation_df <- gather(correlation_df, key = "Cell_Type_2", value = "Correlation", -Cell_Type_1)

# Factor the levels for Cell_Type_2
correlation_df$Cell_Type_2 <- factor(correlation_df$Cell_Type_2, levels = colnames(correlation_matrix))

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(correlation_df, aes(x = Cell_Type_1, y = Cell_Type_2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +  # Add correlation values
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  # scale_fill_gradient2(low = "#451d5b", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title = "Correlation Heatmap",
       x = "Cell Type (In Vivo)",
       y = "Cell Type (In Vitro)") +
  theme(axis.text.x = element_text(angle = -22, hjust = 0))  # Rotate x-axis labels at a 45-degree angle

#i am here

print(heatmap_plot)

av.counts.all.df <- as.data.frame(av.counts.all)
# Define genes of interest
genes_of_interest <- c("Sdc4", "Fgfr3", "Sohlh2", "Sycp3", "Acrv1", "Prm1")

# Convert row names to a column
av.counts.all.df <- av.counts.all.df %>%
  rownames_to_column(var = "gene")

# Filter for selected genes
dotplot_data <- av.counts.all.df %>%
  filter(gene %in% genes_of_interest) %>%
  pivot_longer(cols = -gene, names_to = "cell_type", values_to = "expression")

# Factorize cell types to retain original order
dotplot_data$cell_type <- factor(dotplot_data$cell_type, 
                                 levels = colnames(av.counts.all.df)[-1]) 

DimPlot(somatic.integrated.germ4)
all.markers.somatic <- FindAllMarkers(somatic.integrated.germ4, logfc.threshold = 0, min.pct = 0)



# Define selected genes
# selected_genes <- c("Sdc4", "Fgfr3", "Sohlh2", "Sycp3", "Acrv1", "Prm1")
selected_genes <- c("Sdc4", "Itga9", "Itgb1", "Bmpr1a", "Bmpr1b", "Bmpr2", "Robo1", "Sdc1", "Cd93", "Cd47", "Itga1", "Itga6", "Dag1") #leydig genes
selected_genes <- c("Fgfr1", "Aplp1", "Aplp2", "Lrp10", "Ncstn", "Notch2", "Lrp6", "Ar", "Sort1", "Itga1", "Itga9",
                    "Sdc1", "Ptch1", "Cdon", "Ryr2", "Itgav", "Itgb5", "Tgfbr3", "Acvr1", "Itgb1") #sertoli
selected_genes <- c("Mttp", "App", "Itga6", "Cx3cr1", "Sdc4", "Ldlr", "Dag1", "GPc4", "Sdc1", "Fgfr2", "Fgfr1", "Ptprs",
                    "Mcam", "Gfra1", "Robo1", "Ptprf", "Cd93") #myoid
av.counts.all.df <- as.data.frame(av.counts.all)

av_counts_long <- av.counts.all.df %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell_type", values_to = "expression") %>%
  filter(gene %in% selected_genes)


av_counts_long$cell_type <- gsub("av.counts.", "", av_counts_long$cell_type)


wide_pct_expr <- all.markers.somatic %>%
  filter(gene %in% selected_genes) %>%
  select(gene, cluster, pct.1) %>%
  pivot_wider(names_from = cluster, values_from = pct.1) #%>%

pct_expr_long <- wide_pct_expr %>%
  pivot_longer(-gene, names_to = "cell_type", values_to = "percent_expressed")


dim(av_counts_long)
dim(pct_expr_long)

pct_expr_long <- pct_expr_long %>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))



dotplot_data <- av_counts_long %>%
  left_join(pct_expr_long, by = c("gene", "cell_type"))

cell.types.list <- c("ElongatingSpermatids",
                     "RoundSpermatids",
                     "LateSpermatocytes",
                     "EarlySpermatocytes",
                     "PrelepSpermatocytes",
                     "DiffSpermatogonia",
                     "Progenitors",
                     
                     "SSCs" )


dotplot_data$cell_type <- factor(dotplot_data$cell_type, levels = cell.types.list)
dotplot_data$gene <- factor(dotplot_data$gene, levels = selected_genes)

ggplot(dotplot_data, aes(x = gene, y = cell_type, size = percent_expressed, color = expression)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(0.5, 6)) +  # Dot size based on percent expressed
  scale_color_gradient2(low = "#474747", mid = "red", high = "yellow", midpoint = 0.3) +  
  theme_minimal() +
  labs(title = "Radial Dot Plot of Selected Genes", x = "", y = "") +
  coord_radial(inner.radius = 0.5, r.axis.inside = TRUE, start = 0, end = pi) +  # Circular layout
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Rotate gene labels
  #       axis.title = element_blank(),
  #       panel.grid = element_blank())
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank())



DefaultAssay(somatic.integrated.somatic6) <- "RNA"
VlnPlot(somatic.integrated.somatic6, split.by = "treatment", features = c("Kitl"))
# "Gdnf", "Fgf2", "Fgf11", "Csf1", "Egf", "Csf1"
table(somatic.integrated.somatic6$experiment)








Idents(somatic.integrated.somatic6) <- somatic.integrated.somatic6$cell.type3
table(Idents(somatic.integrated.somatic6))

Idents(somatic.integrated.somatic6) <- factor(somatic.integrated.somatic6$cell.type3, 
                                              levels = c("Sertoli",
                                                         "Epithelial_cells", 
                                                         "Leydig_cells",
                                                         "Peritubular_myoid",
                                                         "Perivascular_smooth_muscle", 
                                                         "Mesenchymal_progenitors", 
                                                         "Endothelial", 
                                                         "Macrophages_peritubular", 
                                                         "Macrophages_interstitial", 
                                                         "Dendritic_cells",
                                                         "T_cells_immunomodulatory",
                                                         "T_cells_effector", 
                                                         "B_cells"
                                                         
                                              ))

list.of.genes <- c("Sox9",
                   #"Ctsl", 
                   "Amhr2", "Clu","Wt1", #sertoli
                   "Epcam", "Krt8", "Krt18","Cftr",
                   "Cyp17a1", 
                   "Star",  "Hsd3b1",#"Actg2",
                   "Col4a3", "Etv1",
                   
                   "Acta2", #myoid
                   "Myh11", 
                   #"Des", 
                   "Pdgfrb",
                   "Mcam", "Col1a1", "Lama2",
                   
                   
                   
                   "Igf1",
                   "Tcf21", #mesenchymal prog
                   
                   
                   "Vwf","Tie1", "Tek", #endothelial
                   
                   "Ptprc",
                   
                   "Ccl4",  "Creb5",
                   
                   "Cd74","H2-Ab1","Csf1r", "Slamf9", #macrophage
                   "Adgre1", 
                   
                   "Ccl8", "F13a1", 
                   "Itgax", #dendritic
                   "Flt3",
                   "Clec9a",
                   
                   
                   #Tcell
                   
                   #"Cd3e", 
                   "Il2ra","Cd3d", 
                   
                   "Cd8a", "Cd4", 
                   
                   #"Tbx21",
                   "Ifng",
                   
                   
                   
                   "Ighm" #B
                   
                   
                   
                   
)


length(list.of.genes)

p1 <- DotPlot(somatic.integrated.somatic6, features = rev(list.of.genes))+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=-0)) + 
  coord_flip()+ scale_colour_gradient2(low = "yellow", mid = "orange2", high = "red3")


# + scale_colour_gradient2(
#   name = waiver(),
#   
#   low = muted("white"),
#   mid = "orange",
#   high = muted("darkred"),
#   # midpoint = 0,
#   space = "Lab",
#   na.value = "grey50",
#   transform = "identity",
#   guide = "colourbar",
#   aesthetics = "colour"
# )

myoid.markers <- FindMarkers(somatic.integrated.somatic6, ident.1 = "Peritubular_myoid")
myoid.markers$pct.diff <- (myoid.markers$pct.1)-(myoid.markers$pct.2)
peritubular.macrophage.markers <- FindMarkers(somatic.integrated.somatic6, ident.1 = "Macrophages_peritubular", ident.2 = "Macrophages_interstitial")
peritubular.macrophage.markers$pct.diff <- (peritubular.macrophage.markers$pct.1)-(peritubular.macrophage.markers$pct.2)
list.of.genes2 <- c(
  # #"undiff"
  "Etv5",
  # 
  "Id4",
  # 
  "Ret",
  "Gfra1",
  "Sdc4",
  # "Mxra7", 
  # "Mpped2", "Tubb5", "Zbtb16", 
  # "Fosl2", 
  # "Rims1",
  "Sohlh1",
  "Crabp1",
  # #diff
  # 
  # "Tbl1x", "Hs6st2", "Ppil1",  "Marcksl1", "Cenpa", "Glis1",
  # "Igf1r","Itih3", "Pygo1","Tcea3", "Ccdc88c", "Usp26",  "Kmt2a",
  # "Gdi1", 
  # 
  # "Sh3rf3", "Jag2", "Colgalt1", "Snhg11", "Asph", "Rhob",
  
  "Elavl2",
  "Kit", 
  # 
  # #early meiotic
  # 
  # "Rpa2", 
  "Dazl",
  
  #"Trank1",
  "Dmc1", "Prss50",
  # "Pbx3", 
  "Prdm9", "Mei1", 
  # "Figla", 
  "Stra8",
  # "Fbxo47","Ccnb3", "Peg3",
  # "Sycp1", 
  "Tex101",
  # "Hfm1",
  
  "Sycp3",
  
  
  "Rad51ap2", 
  # "Ccdc152",
  "Meiob", 
  # "Gml", 
  "Hormad1",
  
  #late meiotic
  
  #"Atr",  #"Tesmin",
  "Piwil2",
  
  # "Marcks",
  # "Msh3", 
  "Spo11",
  # "Ptchd3",
  "Slc25a21",
  # "Dnm3", 
  "Tdrd5", 
  "Piwil1", 
  
  
  # "Abca15", "Abca17", "Gdpd4", "Fbxo43",
  
  "Abca14", "Ccdc178", "Tbpl1",
  "Mns1", 
  
  #round
  # "Pde4d", 
  # "Aox3", #"Syt6",
  # "Tll1",
  "Catsper3", "Fam24a", 
  # "Catsper1", "Ccdc27",  "Acr", 
  "Spaca1"  , "Spag6" , 
  # "Iqca1l", 
  "Acrv1" , "Saxo1" ,
  # "Pipox", "Slc39a2",  "Glt6d1", 
  # "Cd46", "Dyrk4", "Hemgn", "Crb1", "Tmco5b", 
  # "Ttll2",   
  #elongating spermatids
  "Sppl2c", 
  "Tnp1", 
  # "Tnp2", 
  "Prm1", 
  # "Smcp", "Spata3", "M1ap",  "Akap12", 
  "Tppp2", #"Pgk2", #"Odf1", "St6galnac2",
  "Spem1")


length(list.of.genes2)


p2 <- DotPlot(somatic.integrated.germ4, 
              features = rev(list.of.genes2))+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=-0)) + 
  coord_flip()+ scale_colour_gradient2(low = "white", mid = "orange2", high = "red3")


plot_grid(p1, p2)



p <- DotPlot(somatic.integrated.germ4, features = rev(c(
  "Aplp1", "Lrp6", "Sdc4", "Grin2d",
  "Cd46", "Met", "Igf1r", "Igf2r", "Insr", 
  "Itgb5", "Mag", "Itgb1", "Sdc1", "Cd93", "Itga9", 
  "Dag1", "Itga6", "Itga1", "Rpsa"
)))


p + 
  coord_flip() +
  scale_color_viridis_c(option = "C") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#redo the dot plots 

p <- DotPlot(somatic.integrated.germ4, features = rev(c(
  "Aplp1", "Lrp6", "Sdc4", "Grin2d",
  "Cd46", "Met", "Igf1r", "Igf2r", "Insr", 
  "Itgb5", "Mag", "Itgb1", "Sdc1", "Cd93", "Itga9", 
  "Dag1", "Itga6", "Itga1", "Rpsa"
)))



dotplot_leydig <- DotPlot(somatic.integrated.germ4, features = rev(c(
  "Sdc4", "Itga9", "Itgb1", 
  "Bmpr1a", "Bmpr1b", "Bmpr2", "Robo1", "Sdc1", "Cd93", "Cd47",
  # 'Itga9', "Itgb1", 
  "Itga6", "Dag1"
)))

leydig.genes <- c("Bmpr1a", "Bmpr1b", "Bmpr2",
                  "Sdc4", "Itga9", "Itgb1", 
                  "Robo1", "Sdc1", "Cd93", "Cd47",
                  'Itga1', #"Itgb1",
                  "Itga6", "Dag1")



# DotPlot(somatic.integrated.germ4, 
#               features = rev(leydig.genes))+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=-0)) + 
#   coord_flip()+ scale_colour_gradient2(low = "blue", mid = "white", midpoint = 0, high = "red3")


DotPlot(somatic.integrated.germ4,
        features = rev(leydig.genes)) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = -0)) +
  coord_flip() +
  scale_color_gradientn(
    colours = c("lightblue", "white", "yellow", "green", "darkgreen"),
    values = scales::rescale(c(-1, 0.2, 0.5, 1, 2))  
  )


myoid.genes <- c("Sdc1",  "Cd93", "Mcam","Ptprf", "Ptprs",
                 
                 "Sdc4",  "Dag1","Itga6","Robo1",  "Gfra1", 
                 "Fgfr1", "Gpc4",  
                 "Ldlr","App","Fgfr2","Cx3cr1","Mttp")


DotPlot(somatic.integrated.germ4,
        features = rev(myoid.genes)) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = -0)) +
  coord_flip() +
  scale_color_gradientn(
    colours = c("lightblue", "white", "yellow", "green", "darkgreen"),
    values = scales::rescale(c(-1, 0, 0.5, 1, 2)) 
  )

macrophage.genes <- c("Chl1", "Igf1r", "Igf2r",  "Fgfr2",
                      "Sdc2","Ldlr", "Sort1",
                      "Itgb1",  "Art1",
                      "Cd47",   
                      "Grin2d", "Cftr" )


DotPlot(somatic.integrated.germ4,
        features = rev(unique(macrophage.genes))) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = -0)) +
  coord_flip() +
  # scale_color_viridis(option = "D")
  scale_color_gradientn(
    colours = c("lightblue", "white", "yellow", "green", "darkgreen"),
    values = scales::rescale(c(-1, 0.2, 0.5, 1, 2)) 
  )



sertoli.genes <- c("Fgfr1", "Aplp1", "Aplp2", "Lrp10", "Ncstn", "Notch2",
                   "Lrp6", "Rpsa", "Lrp1", "Ar", "Sort1", "Itga1", "Itga9",
                   "Sdc1", "Ptch1", "Smo", "Cdon", "Ryr2", "Itgav", "Itgb5", 
                   "Tgfbr3", "Acvr1", "Itgb1"
)


DotPlot(somatic.integrated.germ4,
        features = rev(unique(sertoli.genes))) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = -0)) +
  coord_flip() +
  scale_color_gradientn(
    colours = c("lightblue", "white", "yellow", "green", "darkgreen"),
    values = scales::rescale(c(-1, 0, 1, 2.5)) 
  )


mesenchymal.genes <- c("Aplp1", "Lrp6", "Sdc4"  ,"Grin2d",
                       "Cd46","Met", "Igf1r", "Igf2r", "Insr", 
                       "Itgb5", "Mag", "Itgb1", "Sdc1", "Cd93", "Itga9", 
                       "Dag1", "Itga6",  "Itga1", "Rpsa"
)


DotPlot(somatic.integrated.germ4,
        features = rev(unique(mesenchymal.genes))) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = -0)) +
  scale_size_continuous(range = c(-1, 8)) +
  coord_flip() +
  scale_color_gradientn(
    # colours = c("lightblue", "white", "orange2", "red3"),
    colours = c("lightblue", "white", "yellow", "green", "darkgreen"),
    values = scales::rescale(c(-1, 0.4, 1, 1.5, 3))  
  )



setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")
load("somatic.integrated.somatic6.Robj")

Idents(somatic.integrated.somatic6) <- "treatment"
table(Idents(somatic.integrated.somatic6))
somatic.integrated.somatic6 <- RenameIdents(somatic.integrated.somatic6, 
                                            'GFP transplanted' = "fertile",
                                            'control' = "infertile")
somatic.integrated.somatic6$fertility <- Idents(somatic.integrated.somatic6)
DimPlot(somatic.integrated.somatic6)
DefaultAssay(somatic.integrated.somatic6) <- "RNA"
Idents(somatic.integrated.somatic6) <- "treatment"
somatic.integrated.somatic6 <- RenameIdents(somatic.integrated.somatic6, 
                                            "GFP transplanted" = "fertile", 
                                            "control" = "infertile"
                                            
)
somatic.integrated.somatic6$fertility <- Idents(somatic.integrated.somatic6)

somatic.integrated.somatic6$fertility <- factor(somatic.integrated.somatic6$fertility,
                                                levels = c(
                                                  "fertile", 
                                                  "infertile"
                                                ))



Sertoli.cells <- subset(somatic.integrated.somatic6, idents = "Sertoli")
DefaultAssay(Sertoli.cells) <- "RNA"
p1 <- VlnPlot(Sertoli.cells, features = c("Amhr2", "Sox9", "Dhh"), split.by = "treatment", pt.size = 0, cols = c( "#c05127", "#48b1e0"))



Idents(somatic.integrated.somatic6) <- "cell.type3"
DimPlot(somatic.integrated.somatic6, label = T)
table(Idents(somatic.integrated.somatic6))
cell.types.mouse <- c("Sertoli",
                      "Epithelial_cells",
                      "Leydig_cells",
                      "Peritubular_myoid",
                      "Perivascular_smooth_muscle",
                      "Mesenchymal_progenitors",
                      "Endothelial",
                      "Macrophages_peritubular",
                      "Macrophages_interstitial",
                      "Dendritic_cells",
                      "T_cells_immunomodulatory",
                      "T_cells_effector",
                      "B_cells"
)
mouse.sertoli.subset <- subset(somatic.integrated.somatic6, idents = "Sertoli")
mouse.mesenchymal.subset <- subset(somatic.integrated.somatic6, idents = "Mesenchymal_progenitors")
mouse.myoid.subset <- subset(somatic.integrated.somatic6, idents = "Peritubular_myoid")
mouse.epithelial.subset <- subset(somatic.integrated.somatic6, idents = "Epithelial_cells")
mouse.smooth.muscle.subset <- subset(somatic.integrated.somatic6, idents = "Perivascular_smooth_muscle")
mouse.Leydig.subset <- subset(somatic.integrated.somatic6, idents = "Leydig_cells")
mouse.macro_inter.subset <- subset(somatic.integrated.somatic6, idents = "Macrophages_interstitial")
mouse.macro_peri.subset <- subset(somatic.integrated.somatic6, idents = "Macrophages_peritubular")
mouse.B.subset <- subset(somatic.integrated.somatic6, idents = "B_cells")
mouse.Tmod.subset <- subset(somatic.integrated.somatic6, idents = "T_cells_immunomodulatory")
mouse.Teff.subset <- subset(somatic.integrated.somatic6, idents = "T_cells_effector")
mouse.dendritic.subset <- subset(somatic.integrated.somatic6, idents = "Dendritic_cells")
mouse.endothelial.subset <- subset(somatic.integrated.somatic6, idents = "Endothelial")

#calculate fold changes for figures



DefaultAssay(mouse.sertoli.subset) <- "RNA"
sertoli.tfs <- c("Nfia", "Nfe2l2", "Sox9",
                 "Sparc", "App", "Psap", "Col6a6",
                 "Dhh", "Gstm7", "Fndc5", "Inha", "Jam2")
Idents(mouse.sertoli.subset) <- "experiment"
table(Idents(mouse.sertoli.subset))
mouse.sertoli.subset.multiome <- subset(mouse.sertoli.subset, idents = "multiome")
Idents(mouse.sertoli.subset.multiome) <- "fertility"
Idents(mouse.sertoli.subset) <- "fertility"
table(Idents(mouse.sertoli.subset.multiome))
sertoli.markers.tfs <- FindMarkers(mouse.sertoli.subset.multiome, ident.1 = "infertile", features = unique(sertoli.tfs))
sertoli.markers.tfs$pct.diff <- sertoli.markers.tfs$pct.1 - sertoli.markers.tfs$pct.2
sertoli.markers.tfs$gene <- rownames(sertoli.markers.tfs)

sertoli.tfs <- c("Nfia", "Nfe2l2", "Sox9",
                 "Sparc", "App", "Psap", "Col6a6",
                 "Dhh", "Gstm7", "Fndc5", "Inha", "Jam2")


sertoli.markers.tfs <- sertoli.markers.tfs %>%
  filter(gene %in% sertoli.tfs) %>%
  mutate(gene = factor(gene, levels = rev(sertoli.tfs)))  # rev() to keep top gene at top of plot


ggplot(sertoli.markers.tfs, aes(x = "Sertoli", y = gene)) +
  geom_point(aes(size = pct.diff, color = avg_log2FC)) +
  scale_color_gradientn(
    colors = c("lightblue", "white", "orange2", "red3"),
    values = scales::rescale(c(-1, 0, 0.7, 1.5), from = c(-1, 1.5)),
    limits = c(-1, 1.5)
  ) +
  scale_size_continuous(range = c(10, 20)) +
  theme_minimal() +
  labs(x = NULL, y = "Gene", color = "avg_log2FC", size = "pct.diff") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())





leydig.tfs <- c("Rora","Esr1",  "Nfia", "Esrrg","Bmp1", "Adam12",  "Lgals1", "Slit3", "Col4a5", "Lama4")
Idents(mouse.Leydig.subset) <- "fertility"
table(Idents(mouse.Leydig.subset))
Idents(mouse.Leydig.subset) <- "experiment"
table(Idents(mouse.Leydig.subset))
mouse.Leydig.subset.multiome <- subset(mouse.Leydig.subset, idents = "multiome")
Idents(mouse.Leydig.subset.multiome) <- "fertility"
leydig.markers.tfs <- FindMarkers(mouse.Leydig.subset.multiome, ident.1 = "infertile", features = leydig.tfs)
leydig.markers.tfs$pct.diff <- leydig.markers.tfs$pct.1 - leydig.markers.tfs$pct.2
leydig.markers.tfs$gene <- rownames(leydig.markers.tfs)


leydig.markers.tfs <- leydig.markers.tfs %>%
  filter(gene %in% leydig.tfs) %>%
  mutate(gene = factor(gene, levels = rev(leydig.tfs)))  # rev() to keep top gene at top of plot


ggplot(leydig.markers.tfs, aes(x = "Sertoli", y = gene)) +
  geom_point(aes(size = pct.diff, color = avg_log2FC)) +
  scale_color_gradientn(
    colors = c("lightblue", "white", "orange2", "red3"),
    values = scales::rescale(c(-1, 0, 0.3, 1.2), from = c(-1, 1.2)),
    limits = c(-1, 1.2)
  ) +
  scale_size_continuous(range = c(10, 20)) +
  theme_minimal() +
  labs(x = NULL, y = "Gene", color = "avg_log2FC", size = "pct.diff") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



mouse.macrophage.subset <- subset(somatic.integrated.somatic6, idents = c("Macrophages_interstitial", "Macrophages_peritubular"))


macrophages.tfs <- c("Hes1",
                     "Nr2f2",
                     "Ebf1",
                     "Junb",
                     "Maf",
                     "Mitf",
                     "Mef2a", "Mef2c", 
                     "Hlf","Creb5", 
                     "Ikzf1",   
                     "Alcam", "Igf1", "Pf4", "Psap","F13a1", "Pdgfb",
                     "Sirpa", "Hspa1a", "Hsp90aa1")
Idents(mouse.macrophage.subset) <- "fertility"
table(Idents(mouse.macrophage.subset))
macrophages.markers.tfs <- FindMarkers(mouse.macrophage.subset, ident.1 = "infertile", features = macrophages.tfs)
macrophages.markers.tfs$pct.diff <- macrophages.markers.tfs$pct.1 - macrophages.markers.tfs$pct.2
macrophages.markers.tfs$gene <- rownames(macrophages.markers.tfs)



macrophages.markers.tfs <- macrophages.markers.tfs %>%
  filter(gene %in% macrophages.tfs) %>%
  mutate(gene = factor(gene, levels = rev(macrophages.tfs)))  # rev() to keep top gene at top of plot


ggplot(macrophages.markers.tfs, aes(x = "macrophages", y = gene)) +
  geom_point(aes(size = pct.diff, color = avg_log2FC)) +
  scale_color_gradientn(colors = c("darkblue", "lightblue", "white", "orange2", "red3"),
                        values = scales::rescale(c(-2,
                                                   -1,
                                                   0,
                                                   1,
                                                   2))) +
  scale_size_continuous(range = c(10, 20)) +  # Make small values small, large values BIG
  theme_minimal() +
  labs(x = NULL, y = "Gene", color = "avg_log2FC", size = "pct.diff") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())





mesenchymal.tfs <- c("Wt1", 
                     "Glis3", 
                     "Irf1", 
                     "Jund", 
                     "Creb5", 
                     "Atf3", 
                     "Klf2",
                     "Egr1", 
                     "Cebpd",
                     "Klf4",
                     "Fosb", 
                     "Fos", 
                     "Junb",
                     "Fhl2",
                     "Pknox2", 
                     "App", 
                     "Igfbp4",
                     "Cxcl12",
                     "Hspa1a",
                     "C3",
                     "Dcn",
                     "Igf1",
                     "Col4a2",
                     "Col3a1",
                     "Col4a4",
                     "Col1a2",
                     "Lamc3",
                     "Lama2"
)

FeaturePlot(mouse.mesenchymal.subset.multiome, features = mesenchymal.tfs, raster = T)

Idents(mouse.mesenchymal.subset) <- "fertility"
table(Idents(mouse.mesenchymal.subset))
Idents(mouse.mesenchymal.subset) <- "fertility"
table(Idents(mouse.mesenchymal.subset))
Idents(mouse.mesenchymal.subset) <- "experiment"
table(Idents(mouse.mesenchymal.subset))
mouse.mesenchymal.subset.multiome <- subset(mouse.mesenchymal.subset, idents = "multiome")
Idents(mouse.mesenchymal.subset.multiome) <- "fertility"

mesenchymal.markers.tfs <- FindMarkers(mouse.mesenchymal.subset.multiome, ident.1 = "infertile", features = unique(mesenchymal.tfs))

mesenchymal.markers.tfs$pct.diff <- mesenchymal.markers.tfs$pct.1 - mesenchymal.markers.tfs$pct.2
mesenchymal.markers.tfs$gene <- rownames(mesenchymal.markers.tfs)


# Ensure genes are a factor in the specified order
mesenchymal.markers.tfs <- mesenchymal.markers.tfs %>%
  filter(gene %in% mesenchymal.tfs) %>%
  mutate(gene = factor(gene, levels = rev(mesenchymal.tfs)))  # rev() to keep top gene at top of plot

# Create the plot
ggplot(mesenchymal.markers.tfs, aes(x = "mesenchymal", y = gene)) +
  geom_point(aes(size = pct.diff, color = avg_log2FC)) +
  scale_color_gradientn(colors = c("darkblue", "lightblue", "white", "orange2", "red3"),
                        values = scales::rescale(c(-1,
                                                   -0.5,
                                                   0,
                                                   0.5,
                                                   1.2))) +
  scale_size_continuous(range = c(10, 20)) +  # Make small values small, large values BIG
  theme_minimal() +
  labs(x = NULL, y = "Gene", color = "avg_log2FC", size = "pct.diff") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())




# Your custom gene order
myoid.tfs <- c("Maf", "Egr1", "Nfib", 
               "Fhl2", "Frem1", "Fos", 
               "Mafb", "Runx1", "Cebpd",
               "Jun", "Irf1", "Rel", "Ar",
               "Klf4", "Irf8")
myoid.ligands <- c("Hgf", "Col1a1","Ntn1", "Nid1",
                   "Tgm2", "Lamc3", "Hspg2", "Lama2",
                   "Lamb1", "Lamc1", "Ncam1", "Fgf2", 
                   "Trf", "Hspa8", "Cx3cl1", "Spon1", "P4hb")
myoid.tfs<-c(myoid.tfs,myoid.ligands)


DefaultAssay(mouse.myoid.subset) <- "RNA"

Idents(mouse.myoid.subset) <- "experiment"
table(Idents(mouse.myoid.subset))
mouse.myoid.subset.multiome <- subset(mouse.myoid.subset, idents = "multiome")
Idents(mouse.myoid.subset.multiome) <- "fertility"
Idents(mouse.myoid.subset) <- "fertility"
table(Idents(mouse.myoid.subset.multiome))
myoid.markers.tfs <- FindMarkers(mouse.myoid.subset.multiome, ident.1 = "infertile", features = unique(myoid.tfs))
myoid.markers.tfs$pct.diff <- myoid.markers.tfs$pct.1 - myoid.markers.tfs$pct.2
myoid.markers.tfs$gene <- rownames(myoid.markers.tfs)

# Ensure genes are a factor in the specified order
myoid.markers.tfs <- myoid.markers.tfs %>%
  filter(gene %in% myoid.tfs) %>%
  mutate(gene = factor(gene, levels = rev(myoid.tfs)))  # rev() to keep top gene at top of plot

# Create the plot
ggplot(myoid.markers.tfs, aes(x = "myoid", y = gene)) +
  geom_point(aes(size = pct.diff, color = avg_log2FC)) +
  scale_color_gradientn(
    colors = c("darkblue", "lightblue", "white", "orange2", "red3", "darkred"),
    values = scales::rescale(c(-2.8 -2, -1.5, 1, 2, 3.6), from = c(-2.8, 3.6)),
    limits = c(-2.8, 3.6)
  ) +
  scale_size_continuous(range = c(10, 20)) +
  theme_minimal() +
  labs(x = NULL, y = "Gene", color = "avg_log2FC", size = "pct.diff") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())







#venn diagrams


sertoli.degs.mouse <- read.csv("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Sertoli_DEGs_mouse_1_infertile_2_fertile.csv")
sertoli.degs.human <- read.csv("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Sertoli_cells_DEGs_human_1_infertile_2_fertile.csv")
head(sertoli.degs.mouse)
head(sertoli.degs.human)

# Example data: use your actual data frames
# (sertoli.degs.human, sertoli.degs.mouse)

# Identify upregulated genes
human_up <- sertoli.degs.human$gene[sertoli.degs.human$avg_log2FC > 0]
mouse_up <- sertoli.degs.mouse$gene[sertoli.degs.mouse$avg_log2FC > 0]

# Identify downregulated genes
human_down <- sertoli.degs.human$gene[sertoli.degs.human$avg_log2FC < 0]
mouse_down <- sertoli.degs.mouse$gene[sertoli.degs.mouse$avg_log2FC < 0]

# Make gene names consistent (uppercase for comparison)
mouse_up_upper <- toupper(mouse_up)
mouse_down_upper <- toupper(mouse_down)


# Overlaps for upregulated
up_overlap <- intersect(human_up, mouse_up_upper)

# Overlaps for downregulated
down_overlap <- intersect(human_down, mouse_down_upper)

# Counts
length(human_up)        # total human upregulated
length(mouse_up_upper)  # total mouse upregulated
length(up_overlap)      # overlap upregulated

length(human_down)      # total human downregulated
length(mouse_down_upper)# total mouse downregulated
length(down_overlap)    # overlap downregulated

venn_list_up <- list(
  Human = human_up,
  Mouse = mouse_up_upper
)

venn_list_down <- list(
  Human = human_down,
  Mouse = mouse_down_upper
)


ggVennDiagram)
ggVennDiagram(venn_list_up, label_alpha = 0) + 
  scale_fill_gradient(low="white", high="red") + 
  labs(title="Upregulated Genes Overlap")




myoid.degs.mouse <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Peritubular_myoid_DEGs_mouse_1_infertile_2_fertile.csv')
myoid.degs.human <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Peritubular_myoid_DEGs_human_1_infertile_2_fertile.csv')
endothelial.degs.mouse <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Endothelial_DEGs_mouse_1_infertile_2_fertile.csv')
endothelial.degs.human <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Endothelial_cells_DEGs_human_1_infertile_2_fertile.csv')
leydig.degs.human <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Leydig_cells_DEGs_human_1_infertile_2_fertile.csv')
leydig.degs.mouse <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Leydig_cells_DEGs_mouse_1_infertile_2_fertile.csv')
macrophages.degs.mouse1 <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Macrophages_peritubular_DEGs_mouse_1_infertile_2_fertile.csv')
macrophages.degs.mouse2 <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Macrophages_interstitial_DEGs_mouse_1_infertile_2_fertile.csv')
macrophages.degs.mouse <- rbind(macrophages.degs.mouse1, macrophages.degs.mouse2)
macrophages.degs.mouse <- macrophages.degs.mouse[!duplicated(macrophages.degs.mouse$gene), ]

macrophages.degs.human <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Macrophages_DEGs_human_1_infertile_2_fertile.csv')
mesenchymal.degs.mouse <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Mesenchymal_progenitors_DEGs_mouse_1_infertile_2_fertile.csv')
mesenchymal.degs.human <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Mesenchymal_progenitors_DEGs_human_1_infertile_2_fertile.csv')
smooth.degs.mouse <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Perivascular_smooth_muscle_DEGs_mouse_1_infertile_2_fertile.csv')
smooth.degs.human <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/Perivascular_smooth_muscle_DEGs_human_1_infertile_2_fertile.csv')
t.cells.degs.mouse1 <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/T_cells_effector_DEGs_mouse_1_infertile_2_fertile.csv')
t.cells.degs.mouse2 <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/T_cells_immunomodulatory_DEGs_mouse_1_infertile_2_fertile.csv')
t.cells.degs.mouse <- rbind(t.cells.degs.mouse1, t.cells.degs.mouse2)
t.cells.degs.mouse <- t.cells.degs.mouse[!duplicated(t.cells.degs.mouse$gene), ]

t.cells.degs.human <- read.csv('/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/finalDEGs/T_cells_DEGs_human_1_infertile_2_fertile.csv')


### without filtering but using intersecting genes ###

# Create empty lists to store plots
up_plots <- list()
down_plots <- list()

# Loop over each cell type
for (cell in cell_types) {
  
  # Dynamically access data frames by name
  human_df <- get(paste0(cell, ".degs.human"))
  mouse_df <- get(paste0(cell, ".degs.mouse"))
  
  # Convert mouse gene names to uppercase for consistent filtering
  mouse_df$gene <- toupper(mouse_df$gene)
  human_df$gene <- toupper(human_df$gene)
  
  #Filter genes to only include those in intersect.genes
  human_df <- human_df[human_df$gene %in% intersect.genes, ]
  mouse_df <- mouse_df[mouse_df$gene %in% intersect.genes, ]
  
  # Identify upregulated genes
  human_up <- human_df$gene[human_df$avg_log2FC > 0]
  mouse_up <- mouse_df$gene[mouse_df$avg_log2FC > 0]
  
  # Identify downregulated genes
  human_down <- human_df$gene[human_df$avg_log2FC < 0]
  mouse_down <- mouse_df$gene[mouse_df$avg_log2FC < 0]
  
  # Note: mouse genes already uppercase
  mouse_up_upper <- mouse_up
  mouse_down_upper <- mouse_down
  
  # Calculate overlaps
  up_overlap <- length(intersect(human_up, mouse_up_upper))
  down_overlap <- length(intersect(human_down, mouse_down_upper))
  
  # Euler diagram data for upregulated
  fit_up <- euler(c(
    Human = length(human_up),
    Mouse = length(mouse_up_upper),
    "Human&Mouse" = up_overlap
  ))
  
  # Euler diagram data for downregulated
  fit_down <- euler(c(
    Human = length(human_down),
    Mouse = length(mouse_down_upper),
    "Human&Mouse" = down_overlap
  ))
  
  # Plot upregulated
  p_up <- plot(fit_up,
               fills = list(fill = c("#E41A1C", "#377EB8"), alpha = 0.6),
               labels = list(font = 2),
               main = paste("Upregulated:", cell))
  
  # Plot downregulated
  p_down <- plot(fit_down,
                 fills = list(fill = c("#4DAF4A", "#FF7F00"), alpha = 0.6),
                 labels = list(font = 2),
                 main = paste("Downregulated:", cell))
  
  # Add to lists
  up_plots[[cell]] <- p_up
  down_plots[[cell]] <- p_down
}

# Combine all upregulated plots
combined_up <- wrap_plots(up_plots, ncol = 4) +
  plot_annotation(title = "Upregulated Genes Overlaps (Filtered by intersect.genes)")

# Combine all downregulated plots
combined_down <- wrap_plots(down_plots, ncol = 4) +
  plot_annotation(title = "Downregulated Genes Overlaps (Filtered by intersect.genes)")

# Display
print(combined_up)
print(combined_down)

# Save
ggsave("combined_upregulated_eulerr_filtered.pdf", combined_up, width = 16, height = 8, dpi = 300)
ggsave("combined_downregulated_eulerr_filtered.pdf", combined_down, width = 16, height = 8, dpi = 300)

message("Filtered proportional Euler diagrams generated and saved!")


#just one set of differentially regulated genes not up or down



# Create an empty list to store plots
diff_plots <- list()

# Loop over each cell type
for (cell in cell_types) {
  
  # Dynamically access data frames
  human_df <- get(paste0(cell, ".degs.human"))
  mouse_df <- get(paste0(cell, ".degs.mouse"))
  
  # Convert mouse gene names to uppercase for consistent filtering
  mouse_df$gene <- toupper(mouse_df$gene)
  human_df$gene <- toupper(human_df$gene)
  
  # Filter genes to only include those in intersect.genes
  human_df <- human_df[human_df$gene %in% intersect.genes, ]
  mouse_df <- mouse_df[mouse_df$gene %in% intersect.genes, ]
  
  # Identify differentially regulated genes (either up or down)
  human_diff <- human_df$gene[human_df$avg_log2FC != 0]
  mouse_diff <- mouse_df$gene[mouse_df$avg_log2FC != 0]
  
  # Calculate overlaps
  diff_overlap <- length(intersect(human_diff, mouse_diff))
  
  print(cat(cell, " human_diff ", length(human_diff), " mouse_diff ", length(mouse_diff), 
            " diff_overlap ", diff_overlap,
            " width ", (length(human_diff)+length(mouse_diff)-(diff_overlap)),
            " percentage overlap ", (diff_overlap/(length(human_diff)+length(mouse_diff))), "%   "
  ))
  
  # Euler diagram data
  fit_diff <- euler(c(
    Human = length(human_diff),
    Mouse = length(mouse_diff),
    "Human&Mouse" = diff_overlap
  ))
  
  # Plot
  p_diff <- plot(fit_diff,
                 fills = list(fill = c("#66C2A5", "#FC8D62"), alpha = 0.6),
                 labels = list(font = 2),
                 main = paste("Differentially Regulated Genes:", cell))
  
  # Add to list
  diff_plots[[cell]] <- p_diff
}

# Combine all differential regulation plots
combined_diff <- wrap_plots(diff_plots, ncol = 4) +
  plot_annotation(title = "Differentially Regulated Genes Overlaps (Filtered by intersect.genes)")

# Display
print(combined_diff)

# Save
ggsave("combined_differentially_regulated_eulerr_filtered.pdf", combined_diff, width = 16, height = 8, dpi = 300)

message("Differentially regulated proportional Euler diagrams generated and saved!")






### load in multiome data ###


#load human data
load("/Users/ewhelan/Desktop/Somatic Project/somatic_samples/somatic.integrated.somatic4_withATAC.Robj")
load("/Users/ewhelan/Desktop/Somatic Project/somatic_samples/somatic.integrated.germ4_withATAC.Robj")
DimPlot(somatic.integrated.somatic4_ATAC_subset)



DefaultAssay(somatic.integrated.somatic4_ATAC_subset) <- "peaks"
frags <- Fragments(somatic.integrated.somatic4_ATAC_subset)  # get list of fragment objects
Fragments(somatic.integrated.somatic4_ATAC_subset) <- NULL  # remove fragment information from assay

# create a vector with all the new paths, in the correct order for your list of fragment objects

new.paths <- list(
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags2/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags3/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags4/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags5/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags6/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags7/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags8/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags9/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags10/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags11/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags12/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags13/atac_fragments.tsv.gz",
  "/Users/ewhelan/Library/CloudStorage/Box-Box/multiomics_with_david/somatic_project_files/frag_files/atac.frags14/atac_fragments.tsv.gz"
)

for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}

Fragments(somatic.integrated.somatic4_ATAC_subset) <- frags # assign updated list back to the object
Fragments(somatic.integrated.somatic4_ATAC_subset)


DefaultAssay(somatic.integrated.somatic4_ATAC_subset) <- "peaks"
somatic.integrated.somatic4_ATAC_subset <- NucleosomeSignal(somatic.integrated.somatic4_ATAC_subset)
somatic.integrated.somatic4_ATAC_subset <- TSSEnrichment(somatic.integrated.somatic4_ATAC_subset)


VlnPlot(somatic.integrated.somatic4_ATAC_subset, features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal"))


somatic.integrated.somatic4_ATAC_subset -> somatic.integrated.somatic7_ATAConly

cells.somatic7 <- Cells(somatic.integrated.somatic7_ATAConly)
length(cells.somatic7)
cells.somatic6 <- Cells(somatic.integrated.somatic6)
length(cells.somatic6)
intersecting.cells <- intersect(cells.somatic6, cells.somatic7)
length(intersecting.cells)

somatic.integrated.somatic7_ATAConly <- subset(somatic.integrated.somatic4_ATAC_subset, cells = intersecting.cells)
somatic.integrated.somatic6_ATAConly <- subset(somatic.integrated.somatic6, cells = intersecting.cells)



somatic.integrated.somatic7_ATAConly@reductions$umap <- somatic.integrated.somatic6_ATAConly@reductions$umap
Idents(somatic.integrated.somatic7_ATAConly) <- somatic.integrated.somatic6_ATAConly$cell.type3

somatic.integrated.somatic7_ATAConly$cell.type3 <- Idents(somatic.integrated.somatic7_ATAConly)


VlnPlot(somatic.integrated.somatic7_ATAConly, features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal"))



somatic.TSS <- somatic.integrated.somatic7_ATAConly$TSS.enrichment
average_somatic.TSS <- mean(somatic.TSS, na.rm = TRUE)
average_somatic.TSS

somatic.nucleosome <- somatic.integrated.somatic7_ATAConly$nucleosome_signal
average_somatic.nucleosome <- mean(somatic.nucleosome, na.rm = TRUE)
average_somatic.nucleosome

somatic.peaks.per.cells <- somatic.integrated.somatic7_ATAConly$nCount_peaks
average_somatic.peaks <- mean(somatic.peaks.per.cells, na.rm = TRUE)
average_somatic.peaks

DimPlot(somatic.integrated.somatic6_ATAConly, group.by = "cell.type3", label = T)
DimPlot(somatic.integrated.somatic7_ATAConly, group.by = "cell.type3", label = T)

CoveragePlot(somatic.integrated.somatic7_ATAConly,
             extend.downstream = 5000,
             extend.upstream = 5000,
             ymax = 100,
             region = "Sox9", 
             features = "Sox9", 
             assay = 'peaks', 
             group.by = "cell.type3",
             split.by = "treatment",
             expression.assay = 'RNA', peaks = F)

list.of.cell.types <- levels(somatic.integrated.somatic6_ATAConly$cell.type3)


for(cell in list.of.cell.types) {
  # cell <- "Sertoli"
  print(cell)
  current.subset <- subset(somatic.integrated.somatic7_ATAConly, idents = cell)
  
  Idents(current.subset) <- "treatment"
  print(table(Idents(current.subset)))
  current.markers.atac <- FindMarkers(current.subset, 
                                      # logfc.threshold = 1,
                                      ident.1 = "control", 
                                      ident.2 = "GFP transplanted")
  
  
  sig_cutoff <- 0.05
  logfc_cutoff <- 0.585
  num_DARs_up <- sum(current.markers.atac$p_val < sig_cutoff & (current.markers.atac$avg_log2FC) > logfc_cutoff, na.rm = TRUE)
  num_DARs_down <- sum(current.markers.atac$p_val < sig_cutoff & (current.markers.atac$avg_log2FC) < -logfc_cutoff, na.rm = TRUE)
  print(num_DARs_up)
  print(num_DARs_down)
}




Idents(somatic.integrated.somatic7_ATAConly) <- "cell.type3"
atac.sertoli <- subset(somatic.integrated.somatic7_ATAConly, idents = "Sertoli")

library(BSgenome.Mmusculus.UCSC.mm10)
seqlevelsStyle(BSgenome.Mmusculus.UCSC.mm10) <- "UCSC"  #"NCBI"
DefaultAssay(atac.sertoli) <- "peaks"
atac.sertoli <- RegionStats(atac.sertoli, genome = BSgenome.Mmusculus.UCSC.mm10, verbose = T)

Idents(atac.sertoli) <- "treatment"
atac.sertoli <- LinkPeaks(
  object = atac.sertoli,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Dhh", "Jam2", "Gstm7","Inha","Col6a6", "Fndc5")
)

CoveragePlot(atac.sertoli,
             extend.downstream = 1000,
             extend.upstream = 30000,
             # ymax = 100,
             region = "Jam2", 
             features = c("Sox9", "Jam2"), 
             assay = 'peaks', 
             group.by = "cell.type3",
             split.by = "treatment",
             expression.assay = 'RNA', peaks = T)

chr.vec <- "chr15"
pos.vec.start <- 98895202
pos.vec.end <- 98896542
gr <- GRanges(
  seqnames = chr.vec,
  ranges = IRanges(start = pos.vec.start, end = pos.vec.end)
)
CoveragePlot(atac.sertoli,
             extend.downstream = 3000,
             extend.upstream = 3000,
             # ymax = 100,
             region = "Dhh", 
             features = c("Sox9", "Dhh"), 
             region.highlight = gr,
             assay = 'peaks', 
             links = F,
             group.by = "cell.type3",
             split.by = "treatment",
             expression.assay = 'RNA', peaks = T) & scale_fill_manual(values = c("#48b1e0", "#c05127"))

CoveragePlot(atac.sertoli,
             extend.downstream = 30000,
             extend.upstream = 30000,
             # ymax = 100,
             region = "Gstm7", 
             features = c("Sox9", "Gstm7"), 
             assay = 'peaks', 
             group.by = "cell.type3",
             split.by = "treatment",
             expression.assay = 'RNA', peaks = T)


CoveragePlot(atac.sertoli,
             extend.downstream = 100000,
             extend.upstream = 100000,
             # ymax = 100,
             region = "Fndc5", 
             features = c("Sox9", "Fndc5"), 
             assay = 'peaks', 
             group.by = "cell.type3",
             split.by = "treatment",
             expression.assay = 'RNA', peaks = T)



CoveragePlot(atac.sertoli,
             extend.downstream = 100000,
             extend.upstream = 100000,
             # ymax = 100,
             region = "Inha", 
             features = c("Sox9", "Inha"), 
             assay = 'peaks', 
             group.by = "cell.type3",
             split.by = "treatment",
             expression.assay = 'RNA', peaks = T)


CoveragePlot(atac.sertoli,
             extend.downstream = 100000,
             extend.upstream = 100000,
             # ymax = 100,
             region = "Col6a6", 
             features = c("Sox9", "Col6a6"), 
             assay = 'peaks', 
             group.by = "cell.type3",
             split.by = "treatment",
             expression.assay = 'RNA', peaks = T)






# Read and preprocess
pathways <- read.csv("/Users/ewhelan/Desktop/somatic_figures/IPA/all.pathways.csv")

# Define desired plot order
cell_order <- c("Sertoli", "leydig", "Mesenchymal", "peritubular_myoid", "myoid_perivascular",
                "endothelial", "macrophage.peritubular", "macrophage.interstitial",
                "dendritic", "T_immuno", "T_effectors")

# Apply transformations
pathways <- pathways %>%
  mutate(
    cell.type = factor(cell.type, levels = cell_order),
    z_score = -z_score,
    p_value_signed = ifelse(z_score > 0, p_value, -p_value)
  )

# Consistent color scale
zlim <- range(pathways$z_score, na.rm = TRUE)

# Define which cell types should use pseudo-log x-axis
log_scaled <- c("peritubular_myoid", "dendritic", "T_effectors")

# Generate plots
plots <- pathways %>%
  arrange(cell.type) %>%
  group_by(cell.type) %>%
  group_split() %>%
  lapply(function(df) {
    p <- ggplot(df, aes(x = p_value_signed, y = reorder(pathway, p_value_signed), fill = z_score)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient2(
        low = "#48b1e0", high = "#c05127", mid = "white", midpoint = 0,
        limits = zlim
      ) +
      labs(title = as.character(df$cell.type[1]),
           x = "-log10(p-value)",
           y = NULL,
           fill = "z-score") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Add pseudo-log x-axis conditionally
    if (as.character(df$cell.type[1]) %in% log_scaled) {
      p <- p + scale_x_continuous(
        trans = pseudo_log_trans(base = 10),
        breaks = pretty_breaks(n = 5)
      )
    }
    
    return(p)
  })

# Combine with shared legend
wrap_plots(plots, ncol = 1, guides = "collect") &
  theme(legend.position = "right")



somatic.integrated.germ4$cell.type3 <- somatic.integrated.germ4$cell.type
merged.somatic <- merge(somatic.integrated.germ4, somatic.integrated.somatic6)
VlnPlot(merged.somatic, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size = 0)


DefaultAssay(somatic.integrated.somatic6) <- "RNA"
Idents(somatic.integrated.somatic6) <- "treatment"
somatic.integrated.fertile <- subset(somatic.integrated.somatic6, idents = "GFP transplanted")
somatic.integrated.infertile <- subset(somatic.integrated.somatic6, idents = "control")
p1 <- FeaturePlot(somatic.integrated.fertile, features = c("Piwil4", "Sycp3", "Tnp1", "Prm1"), ncol = 1) & scale_color_viridis(option = "C") & DarkTheme()
p2 <- FeaturePlot(somatic.integrated.infertile, features = c("Piwil4", "Sycp3", "Tnp1", "Prm1"), ncol = 1) & scale_color_viridis(option = "C") & DarkTheme()

p3 <- FeaturePlot(somatic.integrated.germ4, order = T, features = c("Piwil4", "Sycp3", "Tnp1", "Prm1"), ncol = 1) & scale_color_viridis(option = "C") & DarkTheme()
plot_grid(p3, p2, p1, ncol = 3)


p1 <- DimPlot(somatic.integrated.somatic6, split.by = "replicate", group.by = "treatment", cols = c( "#c05127", "#48b1e0"))
p3 <- DimPlot(somatic.integrated.germ4, split.by = "replicate", group.by = "treatment", cols = c( "#c05127","#48b1e0"))

plot_grid(p1, p3, ncol = 1)



### HUMAN ANALYSIS ###


library(dplyr) #version 1.1.4
library(patchwork) #version 1.2.0
library(gplots) #version 3.1.3.1
library(ggplot2) #version 3.5.1
library(monocle3) #version 1.3.1
library(cowplot) #version 1.1.3
library(viridis) #version 0.6.5
library(tidyverse) #version 2.0.0
library(forcats) #version 1.0.0
library(scCustomize) #version 2.1.2
library(stringr) #version 1.5.1
library(ComplexHeatmap) #version 2.14.0
library(circlize) #version 0.4.16
library(lattice) #version 0.22-6
library(scales) #version 1.3.0
library(ggridges) #version 0.5.6
library(reshape2) #version 1.4.4
library(ggrastr) #version 1.0.2
library(tibble) #version 3.2.1
library(scico) #version 1.5.0
library(RColorBrewer) #version 1.1-3
library(DoubletFinder) #version 2.0.3
library(Seurat) #version 5.1.0
library(SeuratWrappers) #version 0.3.1
library(SeuratObject) #version 5.0.2
library(Signac) #v1.15.0
library(ggrepel) #v0.96

# Single-cell analysis of developing and azoospermia human testicles reveals central role of Sertoli cells
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7655944/
# GEO here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149512


adult1counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/Adult_1.csv", row.names = 1)
adult2counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/Adult_2.csv", row.names = 1)
adult3counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/Adult_3.csv", row.names = 1)
adult4counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/Adult_4.csv", row.names = 1)
adult5counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/Adult_5.csv", row.names = 1)

age2counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/2YO.csv", row.names = 1)
age5counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/5YO.csv", row.names = 1)
age11counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/11YO.csv", row.names = 1)
age17counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/17YO.csv", row.names = 1)

AZFa_Delcounts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/AZFa_Del.csv", row.names = 1)

KS1counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/KS_1.csv", row.names = 1)
KS2counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/KS_2.csv", row.names = 1)
KS3counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/KS_3.csv", row.names = 1)

INOA_1counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/INOA_1.csv", row.names = 1)
INOA_2counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/INOA_2.csv", row.names = 1)
INOA_3counts <- read.csv("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO/GSE149512_RAW/INOA_3.csv", row.names = 1)

list.of.dataframes <- list(adult1counts, adult2counts, adult3counts, adult4counts, adult5counts,
                           age2counts, age5counts, age11counts, age17counts,
                           AZFa_Delcounts, KS1counts, KS2counts, KS3counts, 
                           INOA_1counts, INOA_2counts, INOA_3counts)

adult1 <- CreateSeuratObject(adult1counts, project = "adult1")
adult1$orig.ident <- "adult1"
adult2 <- CreateSeuratObject(adult1counts, project = "adult2")
adult2$orig.ident <- "adult2"
adult3 <- CreateSeuratObject(adult1counts, project = "adult3")
adult3$orig.ident <- "adult3"
adult4 <- CreateSeuratObject(adult1counts, project = "adult4")
adult4$orig.ident <- "adult4"
adult5 <- CreateSeuratObject(adult1counts, project = "adult5")
adult5$orig.ident <- "adult5"

age2 <- CreateSeuratObject(age2counts, project = "age2")
age2$orig.ident <- "age2"
age5 <- CreateSeuratObject(age5counts, project = "age5")
age5$orig.ident <- "age5"
age11 <- CreateSeuratObject(age11counts, project = "age11")
age11$orig.ident <- "age11"
age17 <- CreateSeuratObject(age17counts, project = "age17")
age17$orig.ident <- "age17"

AZFa_Del <- CreateSeuratObject(AZFa_Delcounts, project = "AZFa_Del")
AZFa_Del$orig.ident <- "AZFa_Del"

KS1 <- CreateSeuratObject(KS1counts, project = "KS1")
KS1$orig.ident <- "KS1"
KS2 <- CreateSeuratObject(KS2counts, project = "KS2")
KS2$orig.ident <- "KS2"
KS3 <- CreateSeuratObject(KS3counts, project = "KS3")
KS3$orig.ident <- "KS3"

INOA_1 <- CreateSeuratObject(INOA_1counts, project = "INOA_1")
INOA_1$orig.ident <- "INOA_1"
INOA_2 <- CreateSeuratObject(INOA_2counts, project = "INOA_2")
INOA_2$orig.ident <- "INOA_2"
INOA_3 <- CreateSeuratObject(INOA_3counts, project = "INOA_3")
INOA_3$orig.ident <- "INOA_3"



names.data2 <- c("adult1", "adult2", "adult3", "adult4", "adult5",
                 "age2", "age5", "age11", "age17",
                 "AZFa_Del", "KS1", "KS2", "KS3", 
                 "INOA_1", "INOA_2", "INOA_3")

human.list2 <- list.of.dataframes



for(a in 1:length(human.list2)){
  print(length(Cells(human.list2[[a]])))
}

inflection.point2 <- c(1:length(human.list2))
log_lib_size_at_inflection2 <- c(1:length(human.list2))


par(mfrow = c(3,3))

for(y in 1:length(human.list2)){
  # y=1
  print(names.data2[y])
  expression.df <- as.data.frame(human.list2[[y]])
  
  umi_per_barcode <- colSums(expression.df)
  barcode_rank <- rank(-umi_per_barcode)
  # plot(barcode_rank, umi_per_barcode, 
  #      #ylim=c(1,2000),
  #      xlim=c(1,10000))
  
  log_lib_size <- log10(umi_per_barcode)
  #plot(barcode_rank, log_lib_size, xlim=c(1,10000))
  
  o <- order(barcode_rank)
  log_lib_size <- log_lib_size[o]
  barcode_rank <- barcode_rank[o]
  
  rawdiff <- diff(log_lib_size)/diff(barcode_rank)
  inflection <- which(rawdiff == min(rawdiff[500:4000], na.rm=TRUE))
  
  
  plot(barcode_rank, log_lib_size, xlim=c(1,8000),
       pch = 20,
       main = paste(y, names.data2[y], "cells", inflection, sep="_", "UMIs", round(10^log_lib_size[inflection]))
  )
  
}

min.UMIs2 <- 10^log_lib_size_at_inflection2

#manually set minimum umi numbers

min.UMIs2[1] <- 10^3.7 #HS51
min.UMIs2[2] <- 10^3.9 #HS30H
min.UMIs2[3] <- 10^3.5 #HS27
min.UMIs2[4] <- 10^3.8 #HS26
min.UMIs2[5] <- 10^3.6 #HS25
min.UMIs2[6] <- 10^3.7 #HS22H
min.UMIs2[7] <- 10^4 #HS18H
min.UMIs2[8] <- 10^3.6 #HS14H
min.UMIs2[9] <- 10^3 #HS12H
min.UMIs2[10] <- 10^3.5 #HS3S
min.UMIs2[11] <- 10^3.6 #HS46
min.UMIs2[12] <- 10^3.5 #HS43
min.UMIs2[13] <- 10^3.6 #HS41FT
min.UMIs2[14] <- 10^3.3 #HS40
min.UMIs2[15] <- 10^3.6 #HS31

min.UMIs2[16] <- 10^3.8#"HS49_1", 
min.UMIs2[17] <- 10^3.8#"HS26H", 
# min.UMIs2[18] <- 10^3.8#"HS578_multi"

#other cutoffs

min.features <- 500
max.features <- 8000
max.mito <- 20




list.of.seurats <- c(adult1, adult2, adult3, adult4, adult5,
                     age2, age5, age11, age17,
                     AZFa_Del, KS1, KS2, KS3, 
                     INOA_1, INOA_2, INOA_3)

stage.data <-  c("adult", "adult", "adult", "adult", "adult",
                 "age2", "age5", "age11", "age17",
                 "AZFa_Del", "KS", "KS", "KS", 
                 "INOA", "INOA", "INOA")

# HS578_multi) #26w6d


for(current.sample in 1:length(list.of.seurats)){
  # current.sample <- 1
  current.seurat <- list.of.seurats[[current.sample]]
  # current.seurat$replicate <- names.data2[current.sample]
  # current.seurat$experiment <- "human in vivo"
  current.seurat$stage <- stage.data[current.sample]
  current.seurat[['percent.mito']] <- PercentageFeatureSet(current.seurat, pattern = "^MT.")
  list.of.seurats[[current.sample]] <- current.seurat
}

Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)



for(current.sample in 1:(length(list.of.seurats))){
  # current.sample <- 5
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat <- subset(x = current.seurat, 
                           subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mito < max.mito)
  # subset = percent.mito < max.mito & nCount_RNA > min.UMIs2[current.sample]) # 
  
  current.seurat <- NormalizeData(current.seurat)
  current.seurat <- CellCycleScoring(current.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  list.of.seurats[[current.sample]] <- current.seurat
}

number.of.samples <- length(list.of.seurats)
number.of.cells.per.sample <-  data.frame(
  SampleName = c(1:number.of.samples),
  CellNumber = c(1:number.of.samples),
  GeneNumber = c(1:number.of.samples),
  SampleNumber = c(1:number.of.samples),
  SampleOrigin = c(1:number.of.samples),
  medianUMI = c(1:number.of.samples),
  medianGenes = c(1:number.of.samples),
  medianMito =c(1:number.of.samples),
  cutoff = c(1:number.of.samples)
)

for(a in 1:length(list.of.seurats)){
  number.of.cells.per.sample[a,4]<-a
  number.of.cells.per.sample[a,1]<-names.data2[a]
  number.of.cells.per.sample[a,2]<-length(Cells(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,3]<-length(rownames(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,5]<-stage.data[a]
  number.of.cells.per.sample[a,6]<-median(list.of.seurats[[a]]$nCount_RNA)
  number.of.cells.per.sample[a,7]<-median(list.of.seurats[[a]]$nFeature_RNA)
  number.of.cells.per.sample[a,8]<-median(list.of.seurats[[a]]$percent.mito)
  number.of.cells.per.sample[a,9]<-min.UMIs2[a]
}

print(number.of.cells.per.sample)

Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)




all.merged4 <- Cells.merged

options(future.globals.maxSize = 94 * 1024^3) 
all.merged4 <- SCTransform(all.merged4, vars.to.regress = c("S.Score", "G2M.Score"))

all.merged4 <- RunPCA(all.merged4)
# all.merged4 <- RunUMAP(all.merged4, dims = 1:30)
# DimPlot(all.merged4, reduction = "umap", group.by = c("info", "Phase"))

all.merged4 <- IntegrateLayers(object = all.merged4, method = HarmonyIntegration, normalization.method = "SCT",new.reduction = "integrated.harmony", verbose = F)
all.merged4 <- IntegrateLayers(object = all.merged4, method = RPCAIntegration, normalization.method = "SCT",new.reduction = "integrated.rpca", verbose = F)

all.merged4 <- FindNeighbors(all.merged4, reduction = "integrated.harmony", dims = 1:30)
all.merged4 <- FindClusters(all.merged4, resolution = 0.6)
all.merged4 <- RunUMAP(all.merged4, dims = 1:30, reduction = "integrated.harmony", return.model = TRUE)
DimPlot(all.merged4, group.by = c("orig.ident", "stage"))


all.merged4 <- FindNeighbors(all.merged4, reduction = "integrated.rpca", dims = 1:30)
all.merged4 <- FindClusters(all.merged4, resolution = 0.6)
all.merged4 <- RunUMAP(all.merged4, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.rpca", return.model = TRUE)
DimPlot(all.merged4, reduction = "umap.rpca", group.by = c("orig.ident", "stage"))



FeaturePlot(all.merged4,reduction = "umap.rpca", features = c("DAZL", "FGFR3", "KIT", "SOHLH2", "STRA8", "MEIOB", "SYCP3", "ACVR1", "PRM1", "SOX9", "VWF", "ACTA2"))
setwd("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO")
# save(all.merged4, file = "Zhao_NOA_merged_all.merged4.Robj")
load("Zhao_NOA_merged.Robj")










obj <- Cells.merged
# obj <- all.merged4
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
# obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj))
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("orig.ident"))

options(future.globals.maxSize = 90000 * 1024^2)

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "rpca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
FeaturePlot(obj, features = c("TFAP2C", "KIT", "DDX4", "PIWIL4", "SOHLH2", "MEIOB", "PRM1", "SOX9", "NR5A1"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(obj, features = c("CYP17A1", "TCF21", "VWF", "SOX9"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(obj, features = c("SOX9"), split.by = "stage", reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
DimPlot(obj, group.by = c("stage"), reduction = "umap.rpca")


FeaturePlot(obj, features = "nCount_RNA", reduction = "umap.rpca", max.cutoff = 4000)
VlnPlot(obj, features = "nCount_RNA", group.by = "orig.ident")
human_normal_NOA_subset <- subset(obj, nCount_RNA>4300)
# human_normal_NOA_subset <- obj
for(n.dims in 20:30){
  human_normal_NOA_subset <- FindNeighbors(human_normal_NOA_subset, reduction = "integrated.rpca", dims = 1:n.dims)
  # human_normal_NOA_subset <- FindClusters(human_normal_NOA_subset, resolution = 0.6, cluster.name = "rpca_clusters")
  human_normal_NOA_subset <- RunUMAP(human_normal_NOA_subset, reduction = "integrated.rpca", dims = 1:n.dims, reduction.name = "umap.rpca")
  print(DimPlot(human_normal_NOA_subset, group.by = c("stage", "rpca_clusters"),label = T, reduction = "umap.rpca")+ggtitle(n.dims))
}

n.dims
human_normal_NOA_subset <- FindNeighbors(human_normal_NOA_subset, reduction = "integrated.rpca", dims = 1:30)
# human_normal_NOA_subset <- FindClusters(human_normal_NOA_subset, resolution = 0.6, cluster.name = "rpca_clusters")
human_normal_NOA_subset <- RunUMAP(human_normal_NOA_subset, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
print(DimPlot(human_normal_NOA_subset, group.by = c("stage", "rpca_clusters"),label = T, reduction = "umap.rpca"))#+ggtitle(n.dims))



FeaturePlot(human_normal_NOA_subset, features = c("ACRV1", "KIT", "DDX4", "PIWIL4", "SOHLH2", "MEIOB", "PRM1", "SOX9", "NR5A1"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(human_normal_NOA_subset, features = c("VWF", "CD74", "CD79A", "CYP17A1", "STAR", "ACTA2", "TCF21", "CD3E", "SOX9"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(human_normal_NOA_subset, features = c("TCF21", "LOXL1", "ACTA2", "CYP17A1",
                                                  "DCN", "COL1A1", "MYL6", "HSD3B1"), ncol = 4, reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(human_normal_NOA_subset, features = c("CD3D", "LYZ", "VWF", "MCAM", "STEAP4", "ACTA2", "IGF2"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")


FeaturePlot(human_normal_NOA_subset, features = c("INSL3", "STAR", "CYP17A1", "HSD3B1", "HSD3B2", "CYP11A1", "LHCGR"),reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")


Idents(human_normal_NOA_subset) <- "orig.ident"
table(Idents(human_normal_NOA_subset))



Idents(human_normal_NOA_subset) <- "orig.ident"
human_normal_NOA_subset <- RenameIdents(human_normal_NOA_subset, 
                                        'adult1' = "adult",
                                        'adult2' = "adult",
                                        'adult3' = "adult",
                                        'adult4' = "adult",
                                        'adult5' = "adult",
                                        'age2' = "prepubertal",
                                        'age5' = "prepubertal"  ,
                                        'age11' = "prepubertal",
                                        'age17' = "prepubertal",
                                        
                                        'AZFa_Del' = "infertile",
                                        'KS1' = "infertile",
                                        'KS2' = "infertile",
                                        'KS3' = "infertile",
                                        'INOA_1' = "infertile",
                                        'INOA_2' = "infertile",
                                        'INOA_3' = "infertile")
human_normal_NOA_subset$treatment <- Idents(human_normal_NOA_subset)

Idents(human_normal_NOA_subset) <- "rpca_clusters"
p1 <- DimPlot(human_normal_NOA_subset, reduction = "umap.rpca", label = T)
p2 <- FeaturePlot(human_normal_NOA_subset, features = c("CYP17A1"),reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")

plot_grid(p1, p2)

human_normal_NOA_subset <- RenameIdents(human_normal_NOA_subset, 
                                        'adult1' = "adult",
                                        'adult2' = "adult",
                                        'adult3' = "adult",
                                        'adult4' = "adult",
                                        'adult5' = "adult",
                                        'age2' = "prepubertal",
                                        'age5' = "prepubertal"  ,
                                        'age11' = "prepubertal",
                                        'age17' = "prepubertal",
                                        
                                        'AZFa_Del' = "infertile",
                                        'KS1' = "infertile",
                                        'KS2' = "infertile",
                                        'KS3' = "infertile",
                                        'INOA_1' = "infertile",
                                        'INOA_2' = "infertile",
                                        'INOA_3' = "infertile")
human_normal_NOA_subset$treatment <- Idents(human_normal_NOA_subset)




Idents(human_normal_NOA_subset) <-  "rpca_clusters"
DimPlot(human_normal_NOA_subset, reduction = "umap.rpca", label = T, group.by = "rpca_clusters")

FeaturePlot(human_normal_NOA_subset, order = T, features = c("CD3D", "CD4", "CD8A",
                                                             "LYZ", 
                                                             "VWF", 
                                                             "MCAM",
                                                             "STEAP4", "RERGL", "COL4A1",
                                                             "ACTA2", 
                                                             "IGF2"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(human_normal_NOA_subset, order = T, features = c("FGFR3", 
                                                             "SDC4", 
                                                             "SOHLH2", 
                                                             "SYCP3",
                                                             "MEIOB", 
                                                             "OVOL1",
                                                             "TJP3",
                                                             "ACRV1",
                                                             "TNP1"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(human_normal_NOA_subset, order = T, features = c("SOX9", 
                                                             "CYP17A1", 
                                                             "STAR", 
                                                             "ACTA2",
                                                             "MCAM", 
                                                             "COL1A1",
                                                             "VWF",
                                                             "CSF1R",
                                                             "CD3D"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")



Idents(human_normal_NOA_subset) <- "rpca_clusters"
human_normal_NOA_subset2 <- JoinLayers(human_normal_NOA_subset)



cluster27 <- subset(human_normal_NOA_subset2, idents = 31)
DimPlot(cluster27, reduction = "umap.rpca", label = T, group.by = "rpca_clusters")
cluster27 <- JoinLayers(cluster27)
FeaturePlot(cluster27, features = c("CYP17A1"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")

cluster27 <- FindVariableFeatures(cluster27)
cluster27 <- FindNeighbors(cluster27, reduction = "integrated.rpca", dims = 1:20)
cluster27 <- FindClusters(cluster27, resolution = 0.6, cluster.name = "rpca_clusters")
cluster27 <- RunUMAP(cluster27, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")
DimPlot(cluster27, reduction = "umap.rpca")

cluster27 <- JoinLayers(cluster27)
VlnPlot(cluster27, features = "CYP17A1")
leydig.cells <- subset(cluster27, idents = 1)
not.leydig.cells <- Cells(subset(cluster27, idents = 1, invert = T))

leydig.cells.ids <- Cells(leydig.cells)

# human_normal_NOA_subset <- subset(human_normal_NOA_subset, idents = 13, invert = T)
human_normal_NOA_subset <- subset(human_normal_NOA_subset2, cells = not.leydig.cells, invert = T)
DimPlot(human_normal_NOA_subset, reduction = "umap", group.by = "cell.type")

human_normal_NOA_subset_original <- human_normal_NOA_subset

# human_normal_NOA_subset <- subset(human_normal_NOA_subset, cells = not.leydig.cells, invert = T)
human_normal_NOA_subset <- FindNeighbors(human_normal_NOA_subset, reduction = "integrated.rpca", dims = 1:30)
human_normal_NOA_subset <- FindClusters(human_normal_NOA_subset, resolution = 0.6, cluster.name = "rpca_clusters")
human_normal_NOA_subset <- RunUMAP(human_normal_NOA_subset, reduction = "integrated.rpca", dims = 1:21, reduction.name = "umap.rpca")
print(DimPlot(human_normal_NOA_subset, group.by = c("rpca_clusters"),label = T, reduction = "umap.rpca"))#+ggtitle(n.dims))




Idents(human_normal_NOA_subset) <- "rpca_clusters"
DimPlot(human_normal_NOA_subset, reduction = "umap", label = T, group.by = "rpca_clusters")
human_normal_NOA_subset2 <- JoinLayers(human_normal_NOA_subset)
cluster34markers <- FindMarkers((human_normal_NOA_subset2), ident.1 = "34")
cluster35markers <- FindMarkers((human_normal_NOA_subset2), ident.1 = "35")
cluster31markers <- FindMarkers((human_normal_NOA_subset2), ident.1 = "31")


FeaturePlot(human_normal_NOA_subset, order = T, features = c("SOX9", 
                                                             "EPCAM", 
                                                             "CYP17A1", 
                                                             "CPA3",
                                                             "KRT19", 
                                                             "CLDN10","LHCGR",
                                                             "COL1A1",
                                                             "PAX8",
                                                             "KRT8",
                                                             "CD3D"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")


FeaturePlot(human_normal_NOA_subset, reduction = "umap.rpca",features = c(
  "PRM1", "TNP1", "SOHLH1", "SOHLH2", "STRA8", "MEIOB", "MEIOC", "SYCP1", "PRDM9", "TEX101", "SPO11", "TDRD5", "FAM24A", 
  "MNS1", "ACRV1", "PIWIL1"),  order = T,label = T)
human_normal_NOA_subset <- RenameIdents(human_normal_NOA_subset, 
                                        "7" = "Sertoli_cells", 
                                        "27" = "Sertoli_cells", 
                                        
                                        
                                        "23" = "Undifferentiated_spermatogonia", 
                                        "16" = "Undifferentiated_spermatogonia", 
                                        
                                        "19" = "Differentiating_spermatogonia", 
                                        "29" = "Preleptotene_spermatocytes", 
                                        
                                        "28" = "Early_spermatocytes", 
                                        "30" = "Early_spermatocytes", 
                                        "17" = "Early_spermatocytes", 
                                        "25" = "Early_spermatocytes", 
                                        
                                        
                                        "22" = "Late_spermatocytes", 
                                        "4" = "Late_spermatocytes", 
                                        "9" = "Late_spermatocytes", 
                                        "12" = "Late_spermatocytes", 
                                        "15" = "Late_spermatocytes", 
                                        
                                        "26" = "Round_spermatids", 
                                        "10" = "Round_spermatids", 
                                        "5" = "Round_spermatids", 
                                        "24" = "Round_spermatids", 
                                        "8" = "Round_spermatids", 
                                        "6" = "Round_spermatids", 
                                        
                                        "21" = "Elongating_spermatids", 
                                        "11" = "Elongating_spermatids", 
                                        "20" = "Elongating_spermatids", 
                                        
                                        "1" = "Peritubular_myoid", 
                                        "2" = "Peritubular_myoid", 
                                        
                                        "13" = "Perivascular_smooth_muscle", 
                                        
                                        "33" = "T_cells", 
                                        "18" = "Endothelial_cells", 
                                        "32" = "Endothelial_cells", 
                                        
                                        
                                        "34" = "Elongating_spermatids", 
                                        "35" = "Mast_cells", #marked by CPA3
                                        
                                        "31" = "Leydig_cells",
                                        "14" = "Macrophages", 
                                        
                                        "3" = "Mesenchymal_progenitors", 
                                        "0" = "Mesenchymal_progenitors"
                                        
                                        
)

# Idents(human_normal_NOA_subset, cells = leydig.cells.ids) <- 'Leydig'
human_normal_NOA_subset$cell.type <- Idents(human_normal_NOA_subset)

DimPlot(human_normal_NOA_subset, reduction = "umap", label = T, group.by = "cell.type")

FeaturePlot(human_normal_NOA_subset, reduction = "umap", features = c("KIT", "SOHLH2", "FGFR3", "GFRA1"), order = T)
human_normal_NOA_subset$cell.type <- Idents(human_normal_NOA_subset)

setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")
# save(human_normal_NOA_subset, file = "human_normal_NOA.Robj")
load("human_normal_NOA.Robj")
human_normal_NOA_subset$cell.type <- factor(human_normal_NOA_subset$cell.type,
                                            levels = c(
                                              # "SSCs", 
                                              # "Progenitors", 
                                              "Undifferentiated_spermatogonia", 
                                              "Differentiating_spermatogonia",
                                              "Preleptotene_spermatocytes", 
                                              "Early_spermatocytes", 
                                              "Late_spermatocytes", 
                                              "Round_spermatids", 
                                              "Elongating_spermatids", 
                                              
                                              "Sertoli_cells", 
                                              "Leydig_cells",
                                              "Peritubular_myoid", 
                                              
                                              "Perivascular_smooth_muscle", 
                                              
                                              "Mesenchymal_progenitors",
                                              "Endothelial_cells", 
                                              "Mast_cells",
                                              # "Epithelial_cells", 
                                              "Macrophages", 
                                              "T_cells"
                                              
                                              
                                              
                                              
                                              
                                            ))



col.vector <- c( "chocolate4",
                 # "cornsilk4",
                 "darkseagreen3",
                 "plum4",
                 "paleturquoise",
                 "darkgreen",
                 
                 "blue4",
                 "gold2",
                 
                 "#b12625",
                 "#566b30",
                 "#d9a528",
                 "#cd8163",
                 "#88cdea",
                 # "lightgray",
                 "#7f8080",
                 "#d4a1ca",
                 "#7c57a4",
                 
                 
                 
                 "#5570b6"
                 
                 
)

Idents(human_normal_NOA_subset) <- "treatment"
human_normal_NOA_subset <- RenameIdents(human_normal_NOA_subset, "adult" = "fertile", 
                                        "prepubertal" = "fertile")

human_normal_NOA_subset$fertility <- Idents(human_normal_NOA_subset)



DimPlot(human_normal_NOA_subset, group.by = "cell.type", reduction = "umap", 
        split.by = "fertility", label = T,
        cols = col.vector)


setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")
save(human_normal_NOA_subset, file = "human_normal_NOA_15jun25.Robj")
# load("human_normal_NOA.Robj")
















Idents(human_normal_NOA_subset) <- "cell.type"
DimPlot(human_normal_NOA_subset)
leydig.subset <- subset(human_normal_NOA_subset, idents = "Leydig_cells")
Idents(leydig.subset) <- "treatment"
table(Idents(leydig.subset))
VlnPlot(leydig.subset,pt.size = 0.1, features = c("ESR1", "RORA", "NFIA", "ESRRG",  "ADAM12", "BMP1", "SLIT3", "COL4A5", "LAMA4"), group.by = "treatment")


Idents(human_normal_NOA_subset) <- "cell.type"
DotPlot(human_normal_NOA_subset, features = c("FGF1",
                                              "FGF2",
                                              "FGF3",
                                              "FGF4",
                                              "FGF5",
                                              "FGF6",
                                              "FGF7",
                                              "FGF8",
                                              "FGF9",
                                              "FGF10",
                                              "FGF11",
                                              "FGF12",
                                              "FGF13",
                                              "FGF14",
                                              "FGF15",
                                              "FGF16",
                                              "FGF17",
                                              "FGF18",
                                              "FGF19",
                                              "FGF20",
                                              "FGF21",
                                              "FGF22",
                                              "FGF23",
                                              "FGFR1",
                                              "FGFR2", 
                                              "FGFR3", 
                                              "FGFR4"
                                              
))+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=-0)) + scale_colour_gradient2(
  name = waiver(),
  
  low = muted("white"),
  mid = "lightblue",
  high = muted("darkblue"),
  # midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  transform = "identity",
  guide = "colourbar",
  aesthetics = "colour"
)


human_somatic <- subset(human_normal_NOA_subset, idents = c("spermatogonia", "earlySpermatocytes",
                                                            "lateSpermatocytes", "elongatingSpermatids", 
                                                            "roundSpermatids"), invert = T)

human_somatic$cell.type3 <- human_somatic$cell.type
human_somatic$species <- "human"

DimPlot(human_somatic, reduction = "umap.rpca")
setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")
load("somatic.integrated.somatic6.Robj")

DefaultAssay(somatic.integrated.somatic6) <- "RNA"
gene.names.mouse <- rownames(somatic.integrated.somatic6)
gene.names.mouse[1:5]
gene.names.mouse.upper <- toupper(gene.names.mouse)
gene.names.mouse.upper[1:5]

mouse.assay.data <- somatic.integrated.somatic6@assays$RNA$data
rownames(mouse.assay.data) <- gene.names.mouse.upper
mouse.assay.counts <- somatic.integrated.somatic6@assays$RNA$counts
rownames(mouse.assay.counts) <- gene.names.mouse.upper

humanized.mouse.somatic <- CreateSeuratObject(counts = mouse.assay.data)

humanized.mouse.somatic$cell.type <- somatic.integrated.somatic6$cell.type3

humanized.mouse.somatic@assays$RNA$data <- humanized.mouse.somatic@assays$RNA$counts

humanized.mouse.somatic@assays$RNA$counts <- mouse.assay.counts


humanized.mouse.somatic$species <- "mouse"



Cells.merged <- merge(
  x = humanized.mouse.somatic,
  y = human_somatic)

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        raster = F,
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA")
)




obj <- Cells.merged
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
# obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj))
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")


options(future.globals.maxSize = 90000 * 1024^2)

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "rpca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
FeaturePlot(obj, features = c("TFAP2C", "KIT", "DDX4", "PIWIL4", "SOHLH2", "MEIOB", "PRM1", "SOX9", "NR5A1"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(obj, features = c("CYP17A1", "TCF21", "VWF", "SOX9"), reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(obj, features = c("SOX9"), split.by = "stage", reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
DimPlot(obj, group.by = c("cell.type"),split.by = "species",  label = T, reduction = "umap.rpca")



sertoli.subset <- subset(human_normal_NOA_subset, idents = c(13,28, 6, 9))
myoid.subset <- subset(human_normal_NOA_subset, idents = 13)
Idents(sertoli.subset) <- "treatment"

sertoli.subset <- FindNeighbors(sertoli.subset, reduction = "integrated.rpca", dims = 1:20)
sertoli.subset <- FindClusters(sertoli.subset, resolution = 0.8, cluster.name = "rpca_clusters")
sertoli.subset <- RunUMAP(sertoli.subset, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")



DimPlot(sertoli.subset, reduction = "umap.rpca", split.by = "treatment")
FeaturePlot(sertoli.subset, features = c("SOX9", "CLU", "AMH", "AMHR2"), split.by = "treatment", reduction = "umap.rpca")  & DarkTheme() & scale_color_viridis(option = "C")
Idents(sertoli.subset) <- "treatment"
sertoli.subset <- JoinLayers(sertoli.subset)
markers.infertile <- FindMarkers(sertoli.subset, ident.1 = "adult", ident.2 = "infertile")
setwd("~/Documents/scRNAseq/NOA_GEO_ZHAO")
write.csv(markers.infertile, file="markers.infertile.csv")



table(Idents(sertoli.subset))
sertoli.subset2 <- subset(sertoli.subset, idents = "prepubertal", invert = T)

table(Idents(somatic.integrated.somatic6))
mouse.Sertoli.subset <- subset(somatic.integrated.somatic6, idents = "Sertoli")
Idents(mouse.Sertoli.subset) <- "treatment"

DefaultAssay(mouse.Sertoli.subset) <- "RNA"

mouse.Sertoli.subset <- RenameIdents(mouse.Sertoli.subset, 
                                     'GFP transplanted' = "fertile",
                                     'control' = "infertile")

human.Sertoli.subset <- RenameIdents(sertoli.subset2, 
                                     'adult' = "fertile",
                                     'infertile' = "infertile")

gene.of.interest <- c("Sox9", "Nfia", "Nfe2l2") 

gene.of.interest <- c(#"Sparc", "App", "Psap")#, 
  "Dhh", "Inha", "Jam2")

gene.of.interest <- c("Sparc", "App", "Psap")

gene.of.interest <- c("Jam2", "App", "Psap", "Inha")


p1 <- VlnPlot(sertoli.subset2, pt.size = 0, features = toupper(gene.of.interest), ncol = 1)
p2 <- VlnPlot(mouse.Sertoli.subset, pt.size = 0, features = (gene.of.interest), ncol = 1)
plot_grid(p2, p1, ncol = 2)


DimPlot(mouse.Sertoli.subset)
DimPlot(human.Sertoli.subset)

mouse.Sertoli.fertile <- subset(mouse.Sertoli.subset, idents = "fertile")
mouse.Sertoli.infertile <- subset(mouse.Sertoli.subset, idents = "infertile")
human.Sertoli.fertile <- subset(human.Sertoli.subset, idents = "fertile")
human.Sertoli.infertile <- subset(human.Sertoli.subset, idents = "infertile")

proporiton.mouse.fertile <- (length(WhichCells(mouse.Sertoli.fertile, expression = Nfe2l2 > 0)))/length(Cells(mouse.Sertoli.fertile))
proportion.mouse.infertile <- (length(WhichCells(mouse.Sertoli.infertile, expression = Nfe2l2 > 0)))/length(Cells(mouse.Sertoli.infertile))

proportion.human.fertile <- (length(WhichCells(human.Sertoli.fertile, expression = NFE2L2 > 0)))/length(Cells(human.Sertoli.fertile))
proportion.human.infertile <- (length(WhichCells(human.Sertoli.infertile, expression = NFE2L2 > 0)))/length(Cells(human.Sertoli.infertile))

cat("mouse fertile", proporiton.mouse.fertile, "mouse infertile", proportion.mouse.infertile)
cat("human fertile", proportion.human.fertile, "human infertile", proportion.human.infertile)


gene.of.interest <- c("Fosb", "Egr2", "Zfp36l2") #, "Sparc", "App", "Psap", "Dhh", "Inha", "Jam2")
gene.of.interest <- c("Inhbb", "Inha", "Amh")
p1 <- VlnPlot(sertoli.subset2, features = toupper(gene.of.interest), ncol = 1)
p2 <- VlnPlot(mouse.Sertoli.subset, features = (gene.of.interest), ncol = 1)
plot_grid(p1, p2, ncol = 2)

gene.of.interest.m <- c("Sox9", "Psap", "Gstm7", "Dhh") 
gene.of.interest.h <- c("SOX9", "PSAP", "GSTM2", "DHH") 
p1 <- VlnPlot(sertoli.subset2, pt.size = 0, features = (gene.of.interest.h), ncol = 1)
p2 <- VlnPlot(mouse.Sertoli.subset, pt.size = 0, features = (gene.of.interest.m), ncol = 1)
plot_grid(p2, p1, ncol = 2)




DotPlot(human_normal_NOA_subset, features = c("FGF1",
                                              "FGF2",
                                              "FGF3",
                                              "FGF4",
                                              "FGF5",
                                              "FGF6",
                                              "FGF7",
                                              "FGF8",
                                              "FGF9",
                                              "FGF10",
                                              "FGF11",
                                              "FGF12",
                                              "FGF13",
                                              "FGF14",
                                              "FGF15",
                                              "FGF16",
                                              "FGF17",
                                              "FGF18",
                                              "FGF19",
                                              "FGF20",
                                              "FGF21",
                                              "FGF22",
                                              "FGF23",
                                              "FGFR1",
                                              "FGFR2", 
                                              "FGFR3", 
                                              "FGFR4"
))+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=-0)) + scale_colour_gradient2(
  name = waiver(),
  
  low = muted("white"),
  mid = "lightblue",
  high = muted("darkblue"),
  # midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  transform = "identity",
  guide = "colourbar",
  aesthetics = "colour"
)

Idents(somatic.integrated.somatic6) <- somatic.integrated.somatic6$cell.type3
table(Idents(somatic.integrated.somatic6))

# Idents(somatic.integrated.somatic6) <- factor(somatic.integrated.somatic6$cell.type3, levels = c("B_cells", "T_cells_effector", "T_cells_immunomodulatory", "Dendritic_cells", "Macrophages_interstitial", 
#                                               "Macrophages_peritubular", "Endothelial",  "Mesenchymal_progenitors", "Perivascular_smooth_muscle", 
#                                               "Peritubular_myoid", "Leydig_cells","Epithelial_cells", "Sertoli"))


DimPlot(human_normal_NOA_subset, group.by = "cell.type", reduction = "umap.rpca", 
        split.by = "fertility",
        cols = col.vector)

table(Idents(human_normal_NOA_subset))

human_normal_NOA_subset <- subset(human_normal_NOA_subset, idents = "Epithelial_cells", invert= T)

Idents(human_normal_NOA_subset) <- factor(human_normal_NOA_subset$cell.type, 
                                          levels = c("Sertoli_cells",
                                                     "Leydig_cells",
                                                     "Peritubular_myoid",
                                                     "Perivascular_smooth_muscle",
                                                     "Mesenchymal_progenitors",
                                                     "Endothelial_cells",
                                                     "Macrophages",
                                                     "T_cells",
                                                     
                                                     "SSCs",
                                                     "Progenitors",
                                                     "Differentiating_spermatogonia",
                                                     "Preleptotene_spermatocytes",
                                                     "Early_spermatocytes",
                                                     "Late_spermatocytes",
                                                     "Round_spermatids",
                                                     "Elongating_spermatids"
                                                     
                                          ))

list.of.genes <- c(
  
  "Sox9",
  #"Ctsl",
  #"Amhr2", 
  # "Clu",
  "Wt1", #sertoli
  # "Epcam", "Krt8", "Krt18","Cftr",
  "Cyp17a1",
  "Star",  #"Hsd3b1",#"Actg2",
  #"Col4a3", 
  "Etv1",
  
  "Acta2", #myoid
  "Myh11",
  #"Des",
  "Pdgfrb",
  "Mcam", "Col1a1", "Lama2",
  
  
  
  "Igf1",
  "Tcf21", #mesenchymal prog
  
  
  "Vwf","Tie1", "Tek", #endothelial
  
  "Ptprc",
  
  "Ccl4", # "Creb5",
  
  # "Cd74",
  # "H2-Ab1",
  "Csf1r", #"Slamf9", #macrophage
  # "Adgre1",
  
  # "Ccl8", 
  "F13a1",
  # "Itgax", #dendritic
  # "Flt3",
  # "Clec9a",
  # 
  # 
  # #Tcell
  # 
  # #"Cd3e",
  # "Il2ra",
  "Cd3d",
  
  "Cd8a", "Cd4",
  
  #"Tbx21",
  "Ifng",
  
  
  
  # "Ighm" #B
  
  
  
  
  
  # #"undiff"
  # "Etv5",
  # 
  "Id4",
  # 
  # "Ret",
  "Gfra1",
  # "Sdc4",
  # "Mxra7", 
  # "Mpped2", "Tubb5", "Zbtb16", 
  # "Fosl2", 
  # "Rims1",
  "Sohlh1",
  # "Crabp1",
  # #diff
  # 
  # "Tbl1x", "Hs6st2", "Ppil1",  "Marcksl1", "Cenpa", "Glis1",
  # "Igf1r","Itih3", "Pygo1","Tcea3", "Ccdc88c", "Usp26",  "Kmt2a",
  # "Gdi1", 
  # 
  # "Sh3rf3", "Jag2", "Colgalt1", "Snhg11", "Asph", "Rhob",
  
  # "Elavl2",
  # "Kit", 
  # 
  # #early meiotic
  # 
  # "Rpa2", 
  "Dazl",
  "Prdm9", #"Mei1", 
  #"Trank1",
  #"Dmc1", #"Prss50",
  # "Pbx3", 
  
  # "Figla",
  # "Stra8",
  # "Fbxo47","Ccnb3", "Peg3",
  # "Sycp1", 
  "Tex101",
  # "Hfm1",
  
  # "Sycp3",
  
  
  "Rad51ap2", 
  # "Ccdc152",
  "Meiob", 
  # "Gml", 
  # "Hormad1",
  
  #late meiotic
  
  #"Atr",  #"Tesmin",
  "Piwil2",
  
  # "Marcks",
  # "Msh3", 
  "Spo11",
  # "Ptchd3",
  # "Slc25a21",
  # "Dnm3", 
  "Tdrd5", 
  "Piwil1", 
  
  
  # "Abca15", "Abca17", "Gdpd4", "Fbxo43",
  
  # "Abca14", "Ccdc178", "Tbpl1",
  "Mns1", 
  
  #round
  # "Catsper3", 
  "Fam24a", 
  
  #"Spaca1"  , 
  #"Spag6" , 
  
  "Acrv1" , #"Saxo1" ,
  #elongating spermatids
  # "Sppl2c", 
  "Tnp1",
  
  "Prm1", 
  
  # "Tppp2", 
  "Spem1")


length(list.of.genes2)


DotPlot(human_normal_NOA_subset, 
        features = rev(toupper(list.of.genes)))+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=-0)) + 
  coord_flip()+ scale_colour_gradient2(low = "white", mid = "orange2", high = "red3")



DimPlot(human_normal_NOA_subset, reduction = "umap.rpca")

setwd("/Users/ewhelan/Documents/scRNAseq/NOA_GEO_ZHAO")
# save(human_normal_NOA_subset, file = "human_normal_NOA_subset.05.11.25.Robj")
load("human_normal_NOA_subset.05.11.25.Robj")
human_normal_NOA_subset_backup <- human_normal_NOA_subset

DimPlot(human_normal_NOA_subset_backup, group.by = "seurat_clusters", reduction = "umap.rpca")
human_normal_NOA_subset@reductions$umap <- human_normal_NOA_subset_backup@reductions$umap.rpca
DimPlot(human_normal_NOA_subset, reduction = "umap")

#make a bar chart of DEGs adjusted for spermatid contamination

table(Idents(human_normal_NOA_subset)) 
Idents(human_normal_NOA_subset) <- "treatment"
Idents(human_normal_NOA_subset) <- "cell.type"
table(Idents(human_normal_NOA_subset)) 
Sertoli.subset.human <- subset(human_normal_NOA_subset, idents = "Sertoli_cells")
smoothmuscle.subset.human <- subset(human_normal_NOA_subset, idents = "Perivascular_smooth_muscle")
Tcells.subset.human <- subset(human_normal_NOA_subset, idents = "T_cells")
endothelial.subset.human <- subset(human_normal_NOA_subset, idents = "Endothelial_cells")
Leydig.subset.human <- subset(human_normal_NOA_subset, idents = "Leydig_cells")
macrophage.subset.human <- subset(human_normal_NOA_subset, idents = "Macrophages")
Mesenchymal_progenitors.subset.human <- subset(human_normal_NOA_subset, idents = "Mesenchymal_progenitors")
myoid.subset.human <- subset(human_normal_NOA_subset, idents = "Peritubular_myoid")
mast.subset.human <- subset(human_normal_NOA_subset, idents = "Mast_cells")

Idents(human_normal_NOA_subset) <- "cell.type"
elongating.subset <- subset(human_normal_NOA_subset, idents = "Elongating_spermatids")
Idents(human_normal_NOA_subset) <- "fertility"
DimPlot(human_normal_NOA_subset, split.by = "fertility", reduction = "umap.rpca", group.by = "cell.type")
#make blacklist of spermatid genes
list.of.spermatid.genes.human <- FindMarkers(JoinLayers(human_normal_NOA_subset), ident.1 ="Elongating_spermatids" )



setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")
load("somatic.integrated.somatic6.Robj")

Idents(somatic.integrated.somatic6) <- "treatment"
table(Idents(somatic.integrated.somatic6))
somatic.integrated.somatic6 <- RenameIdents(somatic.integrated.somatic6, 
                                            'GFP transplanted' = "fertile",
                                            'control' = "infertile")
somatic.integrated.somatic6$fertility <- Idents(somatic.integrated.somatic6)
DimPlot(somatic.integrated.somatic6)
DefaultAssay(somatic.integrated.somatic6) <- "RNA"
Idents(somatic.integrated.somatic6) <- "treatment"
somatic.integrated.somatic6 <- RenameIdents(somatic.integrated.somatic6, 
                                            "GFP transplanted" = "fertile", 
                                            "control" = "infertile"
                                            
)
somatic.integrated.somatic6$fertility <- Idents(somatic.integrated.somatic6)

somatic.integrated.somatic6$fertility <- factor(somatic.integrated.somatic6$fertility,
                                                levels = c(
                                                  "fertile", 
                                                  "infertile"
                                                ))
Idents(somatic.integrated.somatic6) <- "cell.type3"
DimPlot(somatic.integrated.somatic6, label = T, split.by = "fertility")
table(Idents(somatic.integrated.somatic6))
cell.types.mouse <- c("Sertoli",
                      "Epithelial_cells",
                      "Leydig_cells",
                      "Peritubular_myoid",
                      "Perivascular_smooth_muscle",
                      "Mesenchymal_progenitors",
                      "Endothelial",
                      "Macrophages_peritubular",
                      "Macrophages_interstitial",
                      "Dendritic_cells",
                      "T_cells_immunomodulatory",
                      "T_cells_effector",
                      "B_cells"
)
mouse.sertoli.subset <- subset(somatic.integrated.somatic6, idents = "Sertoli")
mouse.mesenchymal.subset <- subset(somatic.integrated.somatic6, idents = "Mesenchymal_progenitors")
mouse.myoid.subset <- subset(somatic.integrated.somatic6, idents = "Peritubular_myoid")
mouse.epithelial.subset <- subset(somatic.integrated.somatic6, idents = "Epithelial_cells")
mouse.smooth.muscle.subset <- subset(somatic.integrated.somatic6, idents = "Perivascular_smooth_muscle")
mouse.Leydig.subset <- subset(somatic.integrated.somatic6, idents = "Leydig_cells")
mouse.macro_inter.subset <- subset(somatic.integrated.somatic6, idents = "Macrophages_interstitial")
mouse.macro_peri.subset <- subset(somatic.integrated.somatic6, idents = "Macrophages_peritubular")
mouse.B.subset <- subset(somatic.integrated.somatic6, idents = "B_cells")
mouse.Tmod.subset <- subset(somatic.integrated.somatic6, idents = "T_cells_immunomodulatory")
mouse.Teff.subset <- subset(somatic.integrated.somatic6, idents = "T_cells_effector")
mouse.dendritic.subset <- subset(somatic.integrated.somatic6, idents = "Dendritic_cells")
mouse.endothelial.subset <- subset(somatic.integrated.somatic6, idents = "Endothelial")

list.of.subsets.mouse <- c(mouse.sertoli.subset,
                           mouse.epithelial.subset,
                           mouse.Leydig.subset,
                           mouse.myoid.subset,
                           mouse.smooth.muscle.subset,
                           mouse.mesenchymal.subset,
                           mouse.endothelial.subset,
                           mouse.macro_peri.subset,
                           mouse.macro_inter.subset,
                           mouse.dendritic.subset,
                           mouse.Tmod.subset,
                           mouse.Teff.subset,
                           mouse.B.subset
)



list.of.subsets.human <- c(Sertoli.subset.human,
                           Leydig.subset.human,
                           myoid.subset.human,
                           smoothmuscle.subset.human,
                           Mesenchymal_progenitors.subset.human,
                           
                           endothelial.subset.human,
                           
                           macrophage.subset.human, 
                           
                           Tcells.subset.human,
                           mast.subset.human)

cell.types.human <- c("Sertoli_cells",
                      
                      "Leydig_cells",
                      "Peritubular_myoid",
                      "Perivascular_smooth_muscle",
                      "Mesenchymal_progenitors",
                      "Endothelial_cells",
                      "Macrophages",
                      
                      "T_cells",
                      "Mast_cells"
)

Idents(human_normal_NOA_subset) <- "cell.type"
# elongating.subset <- subset(human_normal_NOA_subset, idents = "Elongating_spermatids")
# Idents(human_normal_NOA_subset) <- "fertility"
#make blacklist of spermatid genes
list.of.spermatid.genes.human <- FindMarkers(JoinLayers(human_normal_NOA_subset), ident.1 ="Elongating_spermatids" )


setwd("/Users/ewhelan/Desktop/Somatic Project/somatic_samples")
load("somatic.integrated.germ4.Robj")
table(Idents(somatic.integrated.germ4))
list.of.spermatid.genes.mouse <- FindMarkers(somatic.integrated.germ4, ident.1 ="ElongatingSpermatids" )

setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")



filtered.genes.human <- list.of.spermatid.genes.human %>%
  filter(avg_log2FC > 1) %>%
  filter(p_val_adj < 0.05)

blacklist.human <- rownames(filtered.genes.human)

filtered.genes.mouse <- list.of.spermatid.genes.mouse %>%
  filter(avg_log2FC > 1) %>%
  filter(p_val_adj < 0.05)

blacklist.mouse <- rownames(filtered.genes.mouse)

cutoff <- 1

DEGs.list.mouse <- list.of.subsets.mouse
DEGs.list.human <- list.of.subsets.human

mouse.genes <- rownames(mouse.sertoli.subset)
mouse.clean.genes <- mouse.genes[!mouse.genes %in% blacklist.mouse]

summary.mouse <- data.frame(
  cell.type = cell.types.mouse,
  upregulated = integer(length(list.of.subsets.mouse)),
  downregulated = integer(length(list.of.subsets.mouse)),
  stringsAsFactors = FALSE
)

setwd("~/Documents/PROJECTS/Somatic cell analysis/finalDEGs")

for(current.loop in 1:length(list.of.subsets.mouse)){
  # current.loop <- 4
  current.subset <- list.of.subsets.mouse[[current.loop]]
  Idents(current.subset) <- "fertility"
  print(table(Idents(current.subset)))
  current.subset <- subset(current.subset, features = mouse.clean.genes)
  current.DEGs <- FindMarkers(current.subset, ident.1 = "infertile", ident.2 = "fertile")
  n_upregulated <- current.DEGs %>%
    filter(avg_log2FC > cutoff) %>%
    filter(p_val_adj < 0.05) %>%
    nrow()
  n_downregulated <- current.DEGs %>%
    filter(avg_log2FC < -cutoff) %>%
    filter(p_val_adj < 0.05) %>%
    nrow()
  print(cat("Upregulated genes:", n_upregulated, "\n", "Downregulated genes:", n_downregulated, "\n"))
  summary.mouse[current.loop, 2] <- n_upregulated
  summary.mouse[current.loop, 3] <- -n_downregulated
  DEGs.list.mouse[[current.loop]] <- current.DEGs
  # write.csv(current.DEGs, file = paste0(cell.types.mouse[current.loop], "_DEGs_mouse_1_infertile_2_fertile.csv"))
  
}
print(summary.mouse)
# write.csv(summary.mouse, file="summary.mouse.DEGs.csv")


#repeat for human #

human.genes <- rownames(Sertoli.subset.human)
human.clean.genes <- human.genes[!human.genes %in% blacklist.human]

summary.human <- data.frame(
  cell.type = cell.types.human,
  upregulated = integer(length(list.of.subsets.human)),
  downregulated = integer(length(list.of.subsets.human)),
  stringsAsFactors = FALSE
)

setwd("~/Documents/PROJECTS/Somatic cell analysis/finalDEGs")

for(current.loop in 1:length(list.of.subsets.human)){
  # current.loop <- 9
  current.subset <- list.of.subsets.human[[current.loop]]
  current.subset <- JoinLayers(current.subset)
  Idents(current.subset) <- "fertility"
  print(table(Idents(current.subset)))
  current.subset <- subset(current.subset, features = human.clean.genes)
  current.DEGs <- FindMarkers(current.subset, ident.1 = "infertile", ident.2 = "fertile")
  n_upregulated <- current.DEGs %>%
    filter(avg_log2FC > cutoff) %>%
    filter(p_val_adj < 0.05) %>%
    nrow()
  n_downregulated <- current.DEGs %>%
    filter(avg_log2FC < -cutoff) %>%
    filter(p_val_adj < 0.05) %>%
    nrow()
  print(cat("Upregulated genes:", n_upregulated, "\n", "Downregulated genes:", n_downregulated, "\n"))
  summary.human[current.loop, 2] <- n_upregulated
  summary.human[current.loop, 3] <- -n_downregulated
  DEGs.list.human[[current.loop]] <- current.DEGs
  # write.csv(current.DEGs, file = paste0(cell.types.human[current.loop], "_DEGs_human_1_infertile_2_fertile.csv"))
  
}
print(summary.human)
# write.csv(summary.human, file="summary.human.DEGs.csv")

###


### violin plots for figure 6 ###


Idents(Sertoli.subset.human) <- "fertility"
Idents(myoid.subset.human) <- "fertility"
Idents(immatureLeydig.subset.human) <- "fertility"


Idents(mouse.sertoli.subset) <- "fertility"
Idents(mouse.mesenchymal.subset) <- "fertility"
Idents(mouse.myoid.subset) <- "fertility"


p1 <- VlnPlot(Sertoli.subset.human, features = toupper(c("Sox9", "Sparc", "App", "Psap", "Dhh")), pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)
p2 <- VlnPlot(mouse.sertoli.subset, features = c("Sox9", "Sparc", "App", "Psap", "Dhh"), pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)
plot_grid(p1, p2)
# sertoli.subset <- JoinLayers(sertoli.subset)
# sertoli.DEGs <- FindMarkers(sertoli.subset, ident.1 = "infertile")

myoid.genes <- c("Irf8", "Cebpd", "Maf", "Rel", "Klf4", "Ar") #cebpd
myoid.genes <- c("Runx1", "Egr1", "Marb", "Irf1", "Jun", "Fozs") #egr1 is good
myoid.genes <- c("Nifb", "Frem1", "Fhl2") #fhl2 has the opposite trend
myoid.genes <- c("Klf12", "Stat1", "Runx1", "Cebpb", "Irf8", "Glis3")
myoid.genes <- c("Mafb", "Fhl2", "Fos", "Fosl2", "Fosb", "Hes1") #fosb is okay
myoid.genes <- c("Junb", "Jund", "Atf3", "Mef2a", "Nfib", "Frem1") #jund
myoid.genes <- c("Srebf2", "Nr4a1", "Klf9", "Rab18", "Jun", "Bhlhe40") #nr4a1 and klf9 

all.downstream <- c(
  "C4b", 
  "Hsp90aa1", 
  "H2-K1",
  "B2m", 
  "Ybx1",
  "Dcn",
  "Tcn2",
  "Penk",
  "Hgf",
  "Pros1",
  "P4hb",
  "Nampt",
  "Fstl1",
  "Nid1",
  "Ltbp1",
  "Ncam1",
  "Ntn1",
  "Hspg2",
  'Wnt5b',
  "Sorbs1",
  "Lamb1",
  "Trf",
  "Col5a1",
  "Spon1",
  "Igf1",
  "Lamc3",
  "Lamc1",
  "Tgm2",
  "Cx3cl1",
  "Adam9",
  "Lama2",
  "Timp3",
  "Fgf2",
  "Sertad1",
  "Hspa8",
  "Vcl",
  "Tln1",
  "Vcan",
  "Col4a2",
  "Col1a2",
  "Col1a1",
  "Col6a1",
  "Col4a4",
  "Col4a3",
  "Col4a1",
  "Col6a3",
  "Col8a1",
  "Col15a1",
  "Col13a1",
  "Vcam1",
  "Nlgn2",
  "Insl3",
  "Lcn2",
  "Sema5a"
  
)



for(a in 1:10){
  print(a)
  p3 <- VlnPlot(myoid.subset.human, features = toupper(all.downstream[((a*5)-4):(a*5)]), pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)
  p4 <- VlnPlot(mouse.myoid.subset, features = all.downstream[((a*5)-4):(a*5)], pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)
  
}

myoid.genes <- c("Egr1", "Cebpd", "Jund", "Hspa8", "Igf1")


p3 <- VlnPlot(myoid.subset.human, features = toupper(myoid.genes), pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)
p4 <- VlnPlot(mouse.myoid.subset, features = myoid.genes, pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)


print(plot_grid(p1, p2, p3, p4, ncol = 4))
# DimPlot(somatic.integrated.somatic6)
# mouse.myoid.subset <- subset(somatic.integrated.somatic6, idents = "Peritubular_myoid")
# human.myoid.subset <- subset(human_normal_NOA_subset, idents = "Sertoli_cells")




mesenchymal.genes <- c("Wt1", "Glis3", "Irf1", "Jund", "Creb5")
mesenchymal.genes <- c("Atf3", "Klf2", "Egr1", "Cebpd", "Klf4") #klf4 is good as is egr1
mesenchymal.genes <- c("Fosb", "Fos", "Junb", "Fhl2", "Pknox2") #fosb is okay 
mesenchymal.genes <- c("App", "Igfbp4", "Cxcl12", "Hspa1a", "C3") #hspa1a is okay, app and igfbp4 also

mesenchymal.genes <- c("Dcn", "Igf1", "Col4a2", "Col3a1", "Col4a4") #igf1 and col4a2 and col3a1 are okay
mesenchymal.genes <- c("Col1a2", "Lamc3", "Lama2")

mesenchymal.genes <- c("Klf4", "Fosb", "Hspa1a", "Igfbp4", "Col3a1")

p5 <- VlnPlot(immatureLeydig.subset.human, features = toupper(mesenchymal.genes), pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)
p6 <- VlnPlot(mouse.mesenchymal.subset, features = mesenchymal.genes, pt.size = 0, cols = c("#c05127", "#48b1e0"), split.by = "fertility", ncol = 1)


print(plot_grid(p1, p2, p3, p4,p5, p6, ncol = 6))




Idents(somatic.integrated.somatic6) <- "cell.type3"
DimPlot(somatic.integrated.somatic6, reduction = "umap")
table(Idents(somatic.integrated.somatic6))
somatic.integrated.somatic6 <- RenameIdents(somatic.integrated.somatic6, 
                                            "T_cells_effector" = "T_cells",
                                            "T_cells_immunomodulatory"= "T_cells",
                                            "Sertoli"= "Sertoli_cells",
                                            "Endothelial"= "Endothelial_cells",
                                            "Macrophages_interstitial"= "Macrophages",
                                            "Macrophages_peritubular"= "Macrophages")
somatic.integrated.somatic6$cell.type4human <- Idents(somatic.integrated.somatic6)
save(somatic.integrated.somatic6, file = "somatic.integrated.somatic6_15jun25.Robj")


VlnPlot(somatic.integrated.somatic6, features = c("Kcnma1"), pt.size = 0, split.by = "fertility")


somatic.integrated.somatic7_ATAConly$cell.type3 <- factor(somatic.integrated.somatic7_ATAConly$cell.type3, levels = c(
  "Sertoli", 
  "Leydig_cells",
  "Peritubular_myoid",
  "Perivascular_smooth_muscle",
  "Mesenchymal_progenitors",
  "Epithelial_cells",
  "Endothelial", 
  "Macrophages_peritubular",
  "Macrophages_interstitial",
  "Dendritic_cells",
  "T_cells_immunomodulatory",
  "T_cells_effector",
  "B_cells"
))

VlnPlot(somatic.integrated.somatic7_ATAConly, 
        group.by = "cell.type3",
        pt.size = 0, 
        cols = c("#b12625",
                 "#566b30",
                 "#d9a528",
                 "#cd8163",
                 "#88cdea",
                 "#7f8080",
                 "#d4a1ca",
                 "#7c57a4",
                 "#6ec6a8",
                 "#9bcb3c",
                 "#ee6463",
                 "#5570b6",
                 "#f06ca8",
                 "#f78c1e"
        ),
        features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal"))


somatic.integrated.somatic6$cell.type3 <- factor(somatic.integrated.somatic6$cell.type3, levels = c(
  "Sertoli", 
  "Leydig_cells",
  "Peritubular_myoid",
  "Perivascular_smooth_muscle",
  "Mesenchymal_progenitors",
  "Epithelial_cells",
  "Endothelial", 
  "Macrophages_peritubular",
  "Macrophages_interstitial",
  "Dendritic_cells",
  "T_cells_immunomodulatory",
  "T_cells_effector",
  "B_cells"
))
VlnPlot(somatic.integrated.somatic6, 
        group.by = "cell.type3",
        pt.size = 0, 
        cols = c("#b12625",
                 "#566b30",
                 "#d9a528",
                 "#cd8163",
                 "#88cdea",
                 "#7f8080",
                 "#d4a1ca",
                 "#7c57a4",
                 "#6ec6a8",
                 "#9bcb3c",
                 "#ee6463",
                 "#5570b6",
                 "#f06ca8",
                 "#f78c1e"
        ),
        features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))


VlnPlot(somatic.integrated.germ4, 
        group.by = "cell.type",
        pt.size = 0, 
        cols = c("chocolate4",
                 "cornsilk4",
                 "darkseagreen3",
                 "plum4",
                 "paleturquoise",
                 "darkgreen",
                 "blue4",
                 "gold2"
        ),
        features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))






### REVISION ###

library(dplyr) #version 1.1.4
library(patchwork) #version 1.2.0
library(gplots) #version 3.1.3.1
library(ggplot2) #version 3.5.1
library(monocle3) #version 1.3.1
library(cowplot) #version 1.1.3
library(viridis) #version 0.6.5
library(tidyverse) #version 2.0.0
library(scCustomize) #version 2.1.2
library(stringr) #version 1.5.1
library(ComplexHeatmap) #version 2.14.0
library(circlize) #version 0.4.16
library(lattice) #version 0.22-6
library(scales) #version 1.3.0
library(ggridges) #version 0.5.6
library(reshape2) #version 1.4.4
library(ggrastr) #version 1.0.2
library(tibble) #version 3.2.1
library(RColorBrewer) #version 1.1-3
library(DoubletFinder) #version 2.0.3
library(Seurat) #version 5.1.0
library(SeuratWrappers) #version 0.3.1
library(SeuratObject) #version 5.0.2
library(Matrix) #version 1.7-3
library(ggalluvial) #v0.12.5


#read in your data, change the paths to match wherever these files are on in your own computer

DHHE7Ae.data <- Read10X_h5('/Users/ewhelan/Documents/scRNAseq/DHH_matrices/DHHE7Ae_filtered_feature_bc_matrix.h5')
DHHE7Au.data <- Read10X_h5('/Users/ewhelan/Documents/scRNAseq/DHH_matrices/DHHE7Au_filtered_feature_bc_matrix.h5')
DHHE7Au2.data <- Read10X_h5('/Users/ewhelan/Documents/scRNAseq/DHH_matrices/DHHE7Au2_filtered_feature_bc_matrix.h5')
DHHKO15.data <- Read10X_h5('/Users/ewhelan/Documents/scRNAseq/DHH_matrices/DHHKO15_filtered_feature_bc_matrix.h5')
DHHKO18.data <- Read10X_h5('/Users/ewhelan/Documents/scRNAseq/DHH_matrices/DHHKO18_filtered_feature_bc_matrix.h5')
mouse0.data <- Read10X_h5('/Users/ewhelan/Documents/scRNAseq/10X_output/10X_mouse_output/WTouts/filtered_feature_bc_matrix.h5')
WWv_A.data.all <- Read10X_h5("/Users/ewhelan/Documents/scRNAseq/DHH_matrices/WWv_control_A9026_filtered_feature_bc_matrix.h5")
WWv_B.data.all <- Read10X_h5("/Users/ewhelan/Documents/scRNAseq/DHH_matrices/WWv_control_A9028_filtered_feature_bc_matrix.h5")
WWv_A.data <- WWv_A.data.all$`Gene Expression`
WWv_B.data <- WWv_B.data.all$`Gene Expression`
stra8ko.data <- Read10X("/Users/ewhelan/Desktop/Dhh_mouse/Stra8/Stra8_filtered_feature_bc_matrix/")

mouse.0.gene.data <- mouse0.data$`Gene Expression`


stra8ko.gene.data <- stra8ko.data$`Gene Expression`

WWv_2_control.data <- Read10X('/Users/ewhelan/Desktop/Somatic Project/somatic.datasets/WWv_2_control/filtered_feature_bc_matrix')
WWv_2_GFP1.data <- Read10X("/Users/ewhelan/Desktop/Somatic Project/somatic.datasets/WWv_2_GFP2/filtered_feature_bc_matrix")
WWv_2_GFP2.data <- Read10X("/Users/ewhelan/Desktop/Somatic Project/somatic.datasets/WWv_2_GFP1/filtered_feature_bc_matrix")


list.of.dataframes <- list(
  DHHE7Ae.data,
  DHHE7Au.data,
  DHHE7Au2.data,
  DHHKO15.data,
  DHHKO18.data,
  mouse.0.gene.data,
  
  WWv_A.data,
  WWv_B.data,
  stra8ko.gene.data,
  WWv_2_control.data,
  WWv_2_GFP1.data,
  WWv_2_GFP2.data
)

dim(WWv_A.data)
dim(WWv_B.data)

names.data2 <- c(  
  "DHHE7Ae",
  "DHHE7Au",
  "DHHE7Au2",
  "DHHKO15",
  "DHHKO18",
  "mouse.0",
  
  "WWv_A",
  "WWv_B",
  "STRA8KO",
  "WWv_2_control",
  "WWv_2_GFP1",
  "WWv_2_GFP2"
)




for(a in 1:length(list.of.dataframes)){
  print(length(Cells(list.of.dataframes[[a]])))
}

inflection.point2 <- c(1:length(list.of.dataframes))
log_lib_size_at_inflection2 <- c(1:length(list.of.dataframes))


par(mfrow = c(3,4))

for(y in 1:length(list.of.dataframes)){
  # y=1
  print(names.data2[y])
  expression.df <- as.data.frame(list.of.dataframes[[y]])
  
  umi_per_barcode <- colSums(expression.df)
  barcode_rank <- rank(-umi_per_barcode)
  # plot(barcode_rank, umi_per_barcode, 
  #      #ylim=c(1,2000),
  #      xlim=c(1,10000))
  
  log_lib_size <- log10(umi_per_barcode)
  #plot(barcode_rank, log_lib_size, xlim=c(1,10000))
  
  o <- order(barcode_rank)
  log_lib_size <- log_lib_size[o]
  barcode_rank <- barcode_rank[o]
  
  rawdiff <- diff(log_lib_size)/diff(barcode_rank)
  inflection <- which(rawdiff == min(rawdiff[500:4000], na.rm=TRUE))
  
  
  plot(barcode_rank, log_lib_size, xlim=c(1,25000),
       pch = 20,
       main = paste(y, names.data2[y], "cells", inflection, sep="_", "UMIs", round(10^log_lib_size[inflection]))
  )
  
  #in theory this code originally was supposed to figure out the inflection point but 
  #it never really worked well so I just set it manually
  
  # abline(v=inflection, col="blue", lwd=2)
  # abline(h=log_lib_size[inflection], col="green", lwd=2)
  # 
  
  # inflection.point2[y] <- inflection
  # log_lib_size_at_inflection2[y] <- log_lib_size[inflection]
}

min.UMIs2 <- 10^log_lib_size_at_inflection2


min.UMIs2[1] <- 10^3
min.UMIs2[2] <- 10^3
min.UMIs2[3] <- 10^3
min.UMIs2[4] <- 10^3
min.UMIs2[5] <- 10^3
min.UMIs2[6] <- 10^3

min.UMIs2[7] <- 10^2.5
min.UMIs2[8] <- 10^2.5
min.UMIs2[9] <- 10^3

min.UMIs2[10] <- 10^3
min.UMIs2[11] <- 10^3
min.UMIs2[12] <- 10^3

# min.UMIs <- 10^3.4
min.features <- 300
max.features <- 6000
max.mito <- 60

rm(list.of.dataframes)



DHHE7Ae <- CreateSeuratObject(DHHE7Ae.data, project = "DHHE7Ae")
DHHE7Au <- CreateSeuratObject(DHHE7Au.data, project = "DHHE7Au")
DHHE7Au2 <- CreateSeuratObject(DHHE7Au2.data, project = "DHHE7Au2")
DHHKO15 <- CreateSeuratObject(DHHKO15.data, project = "DHHKO15")
DHHKO18 <- CreateSeuratObject(DHHKO18.data, project = "DHHKO18")
mouse.0 <- CreateSeuratObject(mouse.0.gene.data, project = "mouse.0")

WWv_A <- CreateSeuratObject(WWv_A.data, project = "WWv_A")
WWv_B <- CreateSeuratObject(WWv_B.data, project = "WWv_B")
STRA8KO <- CreateSeuratObject(stra8ko.gene.data, project = "STRA8KO")

WWv_2_control <- CreateSeuratObject(WWv_2_control.data, project = "WWv_2_control")
WWv_2_GFP1 <- CreateSeuratObject(WWv_2_GFP1.data, project = "WWv_2_GFP1")
WWv_2_GFP2 <- CreateSeuratObject(WWv_2_GFP2.data, project = "WWv_2_GFP2")

list.of.seurats <- c(
  DHHE7Ae,
  DHHE7Au,
  DHHE7Au2,
  DHHKO15,
  DHHKO18,
  mouse.0,
  # mouse.1,
  # mouse.1e,
  # mouse.2,
  # mouse.2e,
  WWv_A, WWv_B,
  STRA8KO,
  WWv_2_control,
  WWv_2_GFP1,
  WWv_2_GFP2
)

genotype.data <-  c("DHHenhancer",
                    "DHHenhancer",
                    "DHHenhancer",
                    
                    
                    "DHHKO",
                    "DHHKO",
                    "WT",
                    # "WT",
                    # "WT",
                    # "WT",
                    # "WT",
                    "WWv",
                    "WWv",
                    "STRA8KO",
                    "WWv",
                    "WWv",
                    "WWv"
                    
)
selection.data <-  c("EpCAM+", "Unselected", "Unselected", "Unselected", 
                     "Unselected", "Unselected", 
                     # "Unselected", "EpCAM+", "Unselected", "EpCAM+", 
                     "Unselected", "Unselected", "Unselected",
                     "Unselected", "Unselected", "Unselected"
                     
                     
)

for(current.sample in 1:length(list.of.seurats)){
  # current.sample <- 1
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat$genotype <- genotype.data[current.sample]
  current.seurat$selection <- selection.data[current.sample]
  current.seurat[['percent.mito']] <- PercentageFeatureSet(current.seurat, pattern = "^mt")
  list.of.seurats[[current.sample]] <- current.seurat
}

Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)






for(current.sample in 1:(length(list.of.seurats))){
  print(current.sample)
  # current.sample <- 13
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat <- subset(x = current.seurat,
                           # subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mito < max.mito & nCount_RNA > min.UMIs)
                           # subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mito < max.mito & nCount_RNA > min.UMIs2[current.sample])
                           subset = percent.mito < max.mito & nCount_RNA > min.UMIs2[current.sample])
  # subset = percent.mito < max.mito & nCount_RNA > 10^3)
  
  # current.seurat <- NormalizeData(current.seurat)
  # current.seurat <- CellCycleScoring(current.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  list.of.seurats[[current.sample]] <- current.seurat
}

number.of.samples <- length(list.of.seurats)
number.of.cells.per.sample <-  data.frame(
  SampleName = c(1:number.of.samples),
  CellNumber = c(1:number.of.samples),
  GeneNumber = c(1:number.of.samples),
  SampleNumber = c(1:number.of.samples),
  SampleOrigin = c(1:number.of.samples),
  medianUMI = c(1:number.of.samples),
  medianGenes = c(1:number.of.samples),
  medianMito =c(1:number.of.samples),
  cutoff = c(1:number.of.samples)
)

for(a in 1:length(list.of.seurats)){
  number.of.cells.per.sample[a,4]<-a
  number.of.cells.per.sample[a,1]<-names.data2[a]
  number.of.cells.per.sample[a,2]<-length(Cells(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,3]<-length(rownames(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,5]<-genotype.data[a]
  number.of.cells.per.sample[a,6]<-median(list.of.seurats[[a]]$nCount_RNA)
  number.of.cells.per.sample[a,7]<-median(list.of.seurats[[a]]$nFeature_RNA)
  number.of.cells.per.sample[a,8]<-median(list.of.seurats[[a]]$percent.mito)
  number.of.cells.per.sample[a,9]<-10^3.3 #min.UMIs2[a]
}

print(number.of.cells.per.sample)




Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)



#Now loading all mouse samples from previous paper. Some of these are duplicated so let's make sure we don't have duplicates
# also all of the somatic paper samples.

setwd("/Users/ewhelan/Desktop/Somatic Project/somatic_samples")
load("somatic.integrated.germ4.Robj")
# load("somatic.integrated.somatic4.Robj")
setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis")
load("somatic.integrated.somatic6.Robj")




setwd("/Users/ewhelan/Documents/scRNAseq/multiome_project_files")
load("mouse.all.combined.2023.04.19.Robj")
DimPlot(mouse.all.combined, label = T)
DefaultAssay(mouse.all.combined) <- "RNA"
FeaturePlot(mouse.all.combined, features = c("Ddx4","Sox9", "Acta2", "Vwf", "Cyp17a1", "Tcf21"), order = T)
# FeaturePlot(mouse.all.combined, features = c("Piwil4"), order = T)
table(mouse.all.combined$orig.ident)

WT_mouse_somatic <- subset(mouse.all.combined, idents = c(20, 25, 22, 28, 26))
DimPlot(somatic.integrated.somatic6, group.by = "cell.type3")
somatic.integrated.somatic6$cell.type <- somatic.integrated.somatic6$cell.type3
DimPlot(somatic.integrated.germ4, group.by = "cell.type")

load("/Users/ewhelan/Documents/PROJECTS/Multiomics MS/Multiomics/Multiomics_MS/multiome_analysis_rebuttal/mouse.rat.merge.final2.Robj")
DimPlot(mouse.rat.merge, group.by = "cell.type.shared")
Idents(mouse.rat.merge) <- "species"
WT.mouse.germ <- subset(mouse.rat.merge, idents = "mouse")
WT.mouse.germ$cell.type <- WT.mouse.germ$cell.type.shared

DefaultAssay(somatic.integrated.somatic6) <- "RNA"
Idents(somatic.integrated.somatic6) <- "cell.type"
table(Idents(somatic.integrated.somatic6))

DefaultAssay(WT_mouse_somatic) <- "RNA"
Idents(WT_mouse_somatic) <- "cell.type"
table(Idents(WT_mouse_somatic))

DefaultAssay(somatic.integrated.germ4) <- "RNA"
Idents(somatic.integrated.germ4) <- "cell.type"
table(Idents(somatic.integrated.germ4))

DefaultAssay(WT.mouse.germ) <- "RNA"
Idents(WT.mouse.germ) <- "cell.type"
table(Idents(WT.mouse.germ))

WT.mouse.germ$genotype <- "WT"
WT_mouse_somatic$genotype <- "WT"

WT.mouse.germ$treatment <- "WT"
WT_mouse_somatic$treatment <- "WT"

# 
# DimPlot(somatic.integrated.somatic6)
# VlnPlot(somatic.integrated.somatic6, features = "Dhh", split.by = "treatment")
table(Cells.merged$genotype)
Cells.merged$treatment <- Cells.merged$genotype

#Add in cells from GUO object

#open GUO et al's data


setwd("~/Desktop/Dhh_mouse/GUO")
dge <- read.delim(
  "GSE112393_MergedAdultMouseST25_DGE.txt",
  check.names = FALSE,
  row.names = 1
)

dim(dge)
head(rownames(dge))
head(colnames(dge))

meta <- read.delim(
  "GSE112393_MergedAdultMouseST25_PerCellAttributes_edited.txt",
  check.names = FALSE
)

dim(meta)
head(meta)
colnames(meta)



rownames(meta) <- meta$CellBarcode

guo.obj <- CreateSeuratObject(
  counts = as.matrix(dge),
  meta.data = meta
)


guo.obj <- NormalizeData(guo.obj)
guo.obj <- FindVariableFeatures(guo.obj)

guo.obj <- ScaleData(guo.obj)
guo.obj <- RunPCA(guo.obj)

guo.obj$genotype <- "WT"
guo.obj$treatment <- "WT"

guo.obj$source <- "Guo"
Cells.merged$source <- "ThisPaper"
somatic.integrated.somatic6$source <- "ThisPaper"
somatic.integrated.germ4$source <- "ThisPaper"
WT.mouse.germ$source <- "MultiomicsPaper"
WT_mouse_somatic$source <- "MultiomicsPaper"

DefaultAssay(somatic.integrated.germ4) <-"RNA"
DefaultAssay(mouse.all.combined) <-"RNA"
DefaultAssay(WT.mouse.germ) <-"RNA"
gene.list.1 <- rownames(DHHE7Ae.data)
gene.list.2 <- rownames(mouse.0.gene.data)
gene.list.3 <- rownames(WWv_A.data)
gene.list.4 <- rownames(WWv_2_control.data)
gene.list.5 <- rownames(somatic.integrated.germ4)
gene.list.6 <- rownames(mouse.all.combined)
gene.list.7 <- rownames(WT.mouse.germ)
gene.list.8 <- rownames(guo.obj)
length(gene.list.1)
length(gene.list.2)
length(gene.list.3)
length(gene.list.4)
length(gene.list.5)
length(gene.list.6)
length(gene.list.7)
length(gene.list.8)
master.gene.list <- Reduce(intersect, list(gene.list.1,gene.list.2,
                                           gene.list.3,gene.list.4,
                                           gene.list.5,gene.list.6,
                                           gene.list.7,gene.list.8))
length(master.gene.list)




# obj <- merge(Cells.merged, somatic.integrated.somatic6)

obj <- merge(
  x = Cells.merged,
  y = c(somatic.integrated.somatic6, somatic.integrated.germ4, 
        WT.mouse.germ, WT_mouse_somatic, guo.obj))
table(obj$orig.ident)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)

obj <- ScaleData(obj)
obj <- RunPCA(obj)


options(future.globals.maxSize = 90000 * 1024^2)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)
obj <- FindNeighbors(obj, reduction = "integrated.harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "harmony_clusters")
obj <- RunUMAP(obj, reduction = "integrated.harmony", dims = 1:30)
FeaturePlot(obj, features = c("Ret", "Sdc4", "Ddx4", "Sohlh1", "Kit", "Stra8", "Meiob","Spo11", "Hormad1","Acrv1", "Prm1",
                              "Vwf", "Acta2", "Sox9", "Cyp17a1", "Ptprc"), order = T)  & DarkTheme() & scale_color_viridis(option = "C")
DimPlot(obj, label = T,  group.by = "orig.ident")
DimPlot(obj, label = T)
table(obj$orig.ident)
FeaturePlot(obj, #features = c("Sox9", "Dhh", "Amh", "Amhr2"), 
            features = c("Sox9", "Clu", "Amh", "nFeature_RNA"), 
            order = F,
            # split.by = "genotype"
) & DarkTheme() & scale_color_viridis(option = "C")

setwd("~/Desktop/Dhh_mouse")
# save(obj, file = "dhh_wt_experimental_wv_all_merged_v3.Robj")
# load("dhh_wt_experimental_wv_all_merged.Robj")

table(obj$orig.ident, obj$genotype)
Idents(obj) <- "harmony_clusters"
DimPlot(obj, label = T)


setwd("~/Desktop/Dhh_mouse")
save(obj, file = "dhh_wt_experimental_wv_all_merged_v5.Robj")
# load("dhh_wt_experimental_wv_all_merged.Robj")

table(obj$orig.ident)
DimPlot(obj, label = T, group.by = "cell.type")
DimPlot(obj, label = T)
all_mouse_merged_final <- subset(obj, idents = c(4, 31, 36, 38), invert = T)



all_mouse_merged_final <- FindNeighbors(all_mouse_merged_final, reduction = "integrated.harmony", dims = 1:30)
all_mouse_merged_final <- FindClusters(all_mouse_merged_final, resolution = 0.8, cluster.name = "harmony_clusters")
all_mouse_merged_final <- RunUMAP(all_mouse_merged_final, reduction = "integrated.harmony", dims = 1:30)
FeaturePlot(all_mouse_merged_final, features = c("Ret", "Sdc4", "Ddx4", "Tex101", "Kit", "Stra8", "Meiob","Spo11", "Hormad1","Acrv1", "Tcf21",
                                                 "Vwf", "Acta2", "Sox9", "Cyp17a1", "Ptprc"), order = T)  & DarkTheme() & scale_color_viridis(option = "C")
DimPlot(all_mouse_merged_final, label = T,  group.by = "cell.type")



Idents(all_mouse_merged_final) <- "harmony_clusters"
DimPlot(all_mouse_merged_final, label = T)

# all_mouse_merged_final <- subset(all_mouse_merged_final, idents = c(41, 38, 40, 16), invert = T)
# all_mouse_merged_final <- subset(all_mouse_merged_final, idents = c(37, 40), invert = T)
all_mouse_merged_final <- subset(all_mouse_merged_final, idents = c(40, 54, 50), invert = T)


all_mouse_merged_final <- FindNeighbors(all_mouse_merged_final, reduction = "integrated.harmony", dims = 1:30)
all_mouse_merged_final <- FindClusters(all_mouse_merged_final, resolution = 1)
all_mouse_merged_final <- RunUMAP(all_mouse_merged_final, reduction = "integrated.harmony", dims = 1:30)

p1 <- DimPlot(all_mouse_merged_final, label = T,  group.by = "cell.type") + NoLegend()
p2 <- DimPlot(all_mouse_merged_final, label = T) + NoLegend()
plot_grid(p1,p2)

FeaturePlot(all_mouse_merged_final, features = c("Ret", "Sdc4", "Mki67", "Tex101", "Kit", "Stra8", "Meiob","Spo11", "Hormad1","Acrv1", "Tcf21",
                                                 "Vwf", "Acta2", "Sox9", "Cyp17a1", "Ptprc"), order = T)  & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(all_mouse_merged_final, features = c("Kit", "Stra8", "Sohlh2", "Sdc4"), order = T)  & DarkTheme() & scale_color_viridis(option = "C")

# table(all_mouse_merged_final$treatment)
# 
Idents(all_mouse_merged_final) <- "seurat_clusters"
DimPlot(all_mouse_merged_final, label = T)

all_mouse_merged_final <- FindSubCluster(
  all_mouse_merged_final,
  cluster = 0,
  graph.name = "RNA_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.2,
  algorithm = 1
)

DimPlot(all_mouse_merged_final, group.by = "sub.cluster", label = T)
Idents(all_mouse_merged_final) <- "sub.cluster"

all_mouse_merged_final <- RenameIdents(all_mouse_merged_final, 
                                       '15' = "SSCs",
                                       '0_3' = "Progenitor_spermatogonia",
                                       '0_0' = "Progenitor_spermatogonia",
                                       # '36' = "Diff.spermatogonia",
                                       '0_2' = "Diff.spermatogonia",
                                       '0_1' = "Diff.spermatogonia",
                                       
                                       '30' = "Preleptotene",
                                       '4' = "Preleptotene",
                                       '32' = "Preleptotene",
                                       
                                       '18' = "Early_spermatocytes",
                                       '14' = "Early_spermatocytes",
                                       
                                       
                                       '23' = "Late_spermatocytes",
                                       '21' = "Late_spermatocytes",
                                       # '27' = "Late_spermatocytes",
                                       # '15' = "Late_spermatocytes",
                                       
                                       '24' = "Round_spermatids",
                                       '12' = "Round_spermatids",
                                       '19' = "Round_spermatids",
                                       # '35' = "Round_spermatids",
                                       
                                       # '33' = "Round_spermatids",
                                       # '22' = "Round_spermatids",
                                       
                                       # '29' = "Elongating_spermatids",
                                       '13' = "Elongating_spermatids",
                                       '39' = "Elongating_spermatids",
                                       '7' = "Elongating_spermatids",
                                       '35' = "Elongating_spermatids",
                                       '5' = "Elongating_spermatids",
                                       '25' = "Elongating_spermatids",
                                       # '6' = "Elongating_spermatids",
                                       '26' = "Elongating_spermatids",
                                       '16' = "Elongating_spermatids",
                                       
                                       '29' = "Sertoli_cells",
                                       '40' = "Sertoli_cells",
                                       '38' = "Sertoli_cells",
                                       '28' = "Sertoli_cells",
                                       '8' = "Sertoli_cells",
                                       '36' = "Sertoli_cells",
                                       # '26' = "Sertoli_cells",
                                       
                                       '34' = "Rete_epithelial",
                                       
                                       '2' = "Leydig_cells",
                                       # '4' = "Leydig_cells",
                                       '6' = "Leydig_cells",
                                       
                                       
                                       '11' = "Peritubular_myoid_cells",
                                       '37' = "Peritubular_myoid_cells",
                                       '1' = "Peritubular_myoid_cells",
                                       # '15' = "Peritubular_myoid_cells",
                                       # '12' = "Peritubular_myoid_cells",
                                       
                                       '20' = "Perivascular_smooth_muscle_cells",
                                       '17' = "Perivascular_smooth_muscle_cells",
                                       # '18' = "Perivascular_smooth_muscle_cells",
                                       
                                       '3' = "Mesenchymal_progenitors",
                                       '22' = "Mesenchymal_progenitors",
                                       # '37' = "Mesenchymal_progenitors",
                                       
                                       '10' = "Endothelial_cells",
                                       # '28' = "Endothelial_cells",
                                       
                                       '9' = "Macrophages",
                                       '31' = "Macrophages",
                                       # '20' = "Macrophages",
                                       '33' = "Dendritic",
                                       
                                       '27' = "T_cells")



DimPlot(all_mouse_merged_final, label = T)+NoLegend()

all_mouse_merged_final$final.cell.type <- Idents(all_mouse_merged_final)




col.vector <- c("#b12625",
                "#566b30",
                "#d9a528",
                "#cd8163",
                "#88cdea",
                "#7f8080",
                "#d4a1ca",
                "#7c57a4",
                "#6ec6a8",
                "#9bcb3c",
                "#ee6463",
                "#5570b6",
                "#f06ca8",
                "#f78c1e"
)


col.vector2 <- c( "chocolate4",
                  "cornsilk4",
                  "darkseagreen3",
                  "plum4",
                  "paleturquoise",
                  "darkgreen",
                  "blue4",
                  "gold2"
)

DimPlot(somatic.integrated.germ4, cols = col.vector2, split.by = "treatment")
DimPlot(somatic.integrated.germ4, cols = col.vector2)
DimPlot(somatic.integrated.somatic6, group.by = "cell.type3", split.by = "treatment", cols = col.vector, label = F)
DimPlot(somatic.integrated.somatic6,  cols = col.vector, label = T)


new.col.vector <- c("chocolate4",
                    "cornsilk4",
                    "darkseagreen3",
                    "plum4",
                    "paleturquoise",
                    "darkgreen",
                    "blue4",
                    "gold2",
                    
                    "#5570b6",
                    "#d4a1ca",
                    "#b12625",
                    "#9bcb3c",
                    
                    "#ee6463",
                    "#6ec6a8",
                    
                    
                    "#f06ca8",
                    "#cd8163",
                    # "#88cdea",
                    
                    "#7f8080",
                    "#d9a528"
                    # "#566b30"
                    
                    
                    
                    
                    
                    
                    
                    # "#7c57a4",
                    
                    
                    
                    
                    
                    # "#f78c1e"
)




DimPlot(all_mouse_merged_final, label = T, raster = F, cols = new.col.vector)


DimPlot(all_mouse_merged_final, 
        group.by = "final.cell.type", label = F, raster = F,
        split.by = "treatment", cols = new.col.vector)


DimPlot(all_mouse_merged_final, 
        group.by = "treatment")

table(all_mouse_merged_final$orig.ident, all_mouse_merged_final$treatment)

#there is a mistake in "treatment"!! 
table(all_mouse_merged_final$orig.ident)



unique(all_mouse_merged_final$orig.ident)
Idents(all_mouse_merged_final) <- "orig.ident"
all_mouse_merged_final <- RenameIdents(all_mouse_merged_final,
                                       # '15' = "SSCs",
                                       
                                       "DHHE7Ae" = "DHHenhancer",
                                       "DHHE7Au" = "DHHenhancer",
                                       "DHHE7Au2" = "DHHenhancer",                      
                                       "DHHKO15"  = "DHHKO",
                                       "DHHKO18"  = "DHHKO",
                                       "mouse.0" = "WT",            
                                       "WWv_A" = "WWv",
                                       "WWv_B" = "WWv",                            
                                       "STRA8KO"   = "STRA8KO",                       
                                       "WWv_2_control" = "Experimental_infertile",
                                       
                                       "WWv_2_GFP1"  = "Experimental_fertile",
                                       "WWv_2_GFP2"                        = "Experimental_fertile",
                                       
                                       "WWv3GFPpositive" = "Experimental_fertile",
                                       "WWv3GFPnegative" = "Experimental_fertile",                  
                                       "WWv3control" = "Experimental_infertile",                     
                                       "WWv5GFPpositive"   = "Experimental_fertile",                
                                       "WWv5GFPnegative"  = "Experimental_fertile",                 
                                       "WWv5control" = "Experimental_infertile",                     
                                       "WWv6GFPpositive"  = "Experimental_fertile",                 
                                       "WWv6GFPnegative"  = "Experimental_fertile",                 
                                       "WWv6control" = "Experimental_infertile",                     
                                       "WWv7GFPnegative"   = "Experimental_fertile",                
                                       "WWv7control" = "Experimental_infertile",                      
                                       "multiome1GFPpos"    = "Experimental_fertile",              
                                       "multiome1GFPnegA"  = "Experimental_fertile",                
                                       "multiome1GFPnegB"  = "Experimental_fertile",                
                                       "multiome1CTL" = "Experimental_infertile",                    
                                       "multiome2GFPpos"    = "Experimental_fertile",               
                                       "multiome2GFPnegA"   = "Experimental_fertile",               
                                       "multiome2GFPnegB"   = "Experimental_fertile",              
                                       "multiome2CTL" = "Experimental_infertile",                     
                                       "multiome3CTLAMHR2posFACS"  = "Experimental_infertile",        
                                       "multiome3transplantedAMHR2posFACS"= "Experimental_fertile",
                                       "multiome3CTLAMHR2negFACS" = "Experimental_infertile",         
                                       "multiome3transplantedAMHR2negFACS" = "Experimental_fertile",
                                       "multiome4CTLAMHR2posMACS"         = "Experimental_infertile",
                                       "multiome4transplantedAMHR2posMACS" = "Experimental_fertile",
                                       "WWv7GFPpositive"    = "Experimental_fertile",               
                                       "mouse1"     = "WT",                      
                                       "mouse2"   = "WT",                         
                                       "mouse3"   = "WT",                         
                                       "mouse4u"  = "WT",                        
                                       "mouse5"   = "WT",                         
                                       "mouse6"   = "WT",                         
                                       "mouse4e"   = "WT",                       
                                       "unselected.mouse.1" = "WT",               
                                       "EpCAM.mouse.1"   = "WT",                  
                                       "unselected.mouse.2"    = "WT",           
                                       "EpCAM.mouse.2"  = "WT",                   
                                       "ST1"      = "WT",                         
                                       "ST2"    = "WT",                          
                                       "ST3"   = "WT",                            
                                       "ST4"   = "WT",                            
                                       "ST5"   = "WT",                           
                                       "ST6"   = "WT",                            
                                       "ST7"   = "WT",                            
                                       "ST8"    = "WT",                          
                                       "SPG1"    = "WT",                          
                                       "SPG2"     = "WT",                         
                                       "SPG3"    = "WT",                         
                                       "INT1"    = "WT",                          
                                       "INT2"   = "WT",                           
                                       "INT3"   = "WT",                          
                                       "INT4"    = "WT",                          
                                       "INT5"    = "WT",                          
                                       "SER1"    = "WT",                         
                                       "SER2"    = "WT",                          
                                       "SER3"    = "WT",                          
                                       "SER4"     = "WT",                        
                                       "SER5"     = "WT",                         
                                       "SER6"     = "WT",                         
                                       "SER7"     = "WT",                        
                                       "SER8"     = "WT",                         
                                       "INT6" = "WT"
)

DimPlot(all_mouse_merged_final)

all_mouse_merged_final$treatment <- Idents(all_mouse_merged_final)


unique(all_mouse_merged_final$orig.ident)
Idents(all_mouse_merged_final) <- "orig.ident"
all_mouse_merged_final <- RenameIdents(all_mouse_merged_final,
                                       # '15' = "SSCs",
                                       
                                       "DHHE7Ae" = "EpCAM",
                                       "DHHE7Au" = "Unselected",
                                       "DHHE7Au2" = "Unselected",                      
                                       "DHHKO15"  = "Unselected",
                                       "DHHKO18"  = "Unselected",
                                       "mouse.0" = "Unselected",            
                                       "WWv_A" = "Unselected",
                                       "WWv_B" = "Unselected",                            
                                       "STRA8KO"   = "Unselected",                       
                                       "WWv_2_control" = "Unselected",
                                       
                                       "WWv_2_GFP1"  = "Unselected",
                                       "WWv_2_GFP2"                        = "Unselected",
                                       
                                       "WWv3GFPpositive" = "GFP_pos",
                                       "WWv3GFPnegative" = "GFP_neg",                  
                                       "WWv3control" = "Unselected",                     
                                       "WWv5GFPpositive"   = "GFP_pos",                
                                       "WWv5GFPnegative"  = "GFP_neg",                 
                                       "WWv5control" = "Unselected",                     
                                       "WWv6GFPpositive"  = "GFP_pos",                 
                                       "WWv6GFPnegative"  = "GFP_neg",                 
                                       "WWv6control" = "Unselected",                     
                                       "WWv7GFPnegative"   = "GFP_neg",                
                                       "WWv7control" = "Unselected",                      
                                       "multiome1GFPpos"    = "GFP_pos",              
                                       "multiome1GFPnegA"  = "GFP_neg",                
                                       "multiome1GFPnegB"  = "GFP_neg",                
                                       "multiome1CTL" = "Unselected",                    
                                       "multiome2GFPpos"    = "GFP_pos",               
                                       "multiome2GFPnegA"   = "GFP_neg",               
                                       "multiome2GFPnegB"   = "GFP_neg",              
                                       "multiome2CTL" = "Unselected",                     
                                       "multiome3CTLAMHR2posFACS"  = "AMHR2pos",        
                                       "multiome3transplantedAMHR2posFACS"= "AMHR2pos",
                                       "multiome3CTLAMHR2negFACS" = "AMHR2neg",         
                                       "multiome3transplantedAMHR2negFACS" = "AMHR2neg",
                                       "multiome4CTLAMHR2posMACS"         = "AMHR2pos",
                                       "multiome4transplantedAMHR2posMACS" = "AMHR2pos",
                                       "WWv7GFPpositive"    = "GFP_pos",               
                                       "mouse1"     = "Unselected",                      
                                       "mouse2"   = "Unselected",                         
                                       "mouse3"   = "Unselected",                         
                                       "mouse4u"  = "Unselected",                        
                                       "mouse5"   = "Unselected",                         
                                       "mouse6"   = "Unselected",                         
                                       "mouse4e"   = "EpCAM",                       
                                       "unselected.mouse.1" = "Unselected",               
                                       "EpCAM.mouse.1"   = "EpCAM",                  
                                       "unselected.mouse.2"    = "Unselected",           
                                       "EpCAM.mouse.2"  = "EpCAM",                   
                                       "ST1"      = "Unselected",                         
                                       "ST2"    = "Unselected",                          
                                       "ST3"   = "Unselected",                            
                                       "ST4"   = "Unselected",                            
                                       "ST5"   = "Unselected",                           
                                       "ST6"   = "Unselected",                            
                                       "ST7"   = "Unselected",                            
                                       "ST8"    = "Unselected",                          
                                       "SPG1"    = "SPG",                          
                                       "SPG2"     = "SPG",                         
                                       "SPG3"    = "SPG",                         
                                       "INT1"    = "INT",                          
                                       "INT2"   = "INT",                           
                                       "INT3"   = "INT",                          
                                       "INT4"    = "INT",                          
                                       "INT5"    = "INT",                          
                                       "SER1"    = "SER",                         
                                       "SER2"    = "SER",                          
                                       "SER3"    = "SER",                          
                                       "SER4"     = "SER",                        
                                       "SER5"     = "SER",                         
                                       "SER6"     = "SER",                         
                                       "SER7"     = "SER",                        
                                       "SER8"     = "SER",                         
                                       "INT6" = "INT"
)
DimPlot(all_mouse_merged_final)
all_mouse_merged_final$selection <- Idents(all_mouse_merged_final)

DimPlot(all_mouse_merged_final, 
        group.by = "final.cell.type", label = T, 
        raster = F,
        # split.by = "treatment", 
        cols = new.col.vector)





setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
# save(all_mouse_merged_final, file = "dhh_wt_experimental_wv_all_merged_v6.Robj")
# load("dhh_wt_experimental_wv_all_merged_v6.Robj")
load("dhh_wt_experimental_wv_all_merged_v6.Robj")
Idents(all_mouse_merged_final) <- "final.cell.type"
table(Idents(all_mouse_merged_final))
all_mouse_merged_final_sertoli <- subset(all_mouse_merged_final, idents = "Sertoli_cells")
VlnPlot(all_mouse_merged_final_sertoli, features = "Dhh", group.by = "treatment")

DimPlot(all_mouse_merged_final, group.by = "selection")

setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
# save(all_mouse_merged_final, file = "dhh_wt_experimental_wv_all_merged_v6.Robj")
# load("dhh_wt_experimental_wv_all_merged_v6.Robj")
load("dhh_wt_experimental_wv_all_merged_v6.Robj")


all_mouse_merged_final$treatment <- factor(all_mouse_merged_final$treatment, levels = c(
  "WT", "WWv", 
  "Experimental_fertile",
  "Experimental_infertile",
  "DHHKO",
  "DHHenhancer",
  "STRA8KO"
  
))


p <- DimPlot(
  all_mouse_merged_final,
  label = FALSE,
  raster = F,
  split.by = "treatment",
  group.by = "final.cell.type",
  cols = new.col.vector
)

ggsave(
  filename = "all_mouse_merged_final_split_by_treatment.pdf",
  plot = p,
  width = 24,
  height = 4,
  units = "in",
  limitsize = FALSE
)

table(all_mouse_merged_final$treatment, all_mouse_merged_final$final.cell.type)

table(all_mouse_merged_final$selection)


wwv <- c(0, 0.0243309)

experimental <- c(
  0.208707361,
  11.22988506,
  0,
  0,
  0.127144536
)

wilcox.test(
  wwv,
  experimental,
  exact = FALSE
)




all_mouse_merged_final










Idents(all_mouse_merged_final) <- "orig.ident"
all_mouse_merged_unselected <- subset(all_mouse_merged_final, idents = c(
  'SPG1',
  'SPG2',
  'SPG3',
  
  'INT1',
  'INT2',
  'INT3',
  'INT4',
  'INT5',
  'INT6',
  
  'SER1',
  'SER2',
  'SER3',
  'SER4',
  'SER5',
  'SER6',
  'SER7',
  'SER8'), invert = T)

Idents(all_mouse_merged_unselected) <- "selection"
all_mouse_merged_unselected <- subset(all_mouse_merged_unselected, idents = c("EpCAM",        "EpCAM+"), invert = T)

#note that other selection methods still exist but since they're not in the WT or WWv, ignore for now

table(all_mouse_merged_unselected$orig.ident, all_mouse_merged_unselected$final.cell.type)



table(all_mouse_merged_final$orig.ident)
FeaturePlot(all_mouse_merged_final, #features = c("Sox9", "Dhh", "Amh", "Amhr2"), 
            features = c("Sox9", "Clu", "Amh", "nFeature_RNA"), 
            order = F,
            # split.by = "genotype"
) & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(all_mouse_merged_final, #features = c("Sox9", "Dhh", "Amh", "Amhr2"), 
            features = c("Etv5", "Gfra1"), 
            order = F,
            # split.by = "genotype"
) & DarkTheme() & scale_color_viridis(option = "C")

setwd("~/Desktop/Dhh_mouse")
save(all_mouse_merged_final, file = "dhh_wt_experimental_wv_all_merged_v4.Robj")


spermatogonia.subset <- subset(all_mouse_merged_final, idents = 16)


setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
# save(all_mouse_merged_final, file = "dhh_wt_experimental_wv_all_merged_v6.Robj")
# load("dhh_wt_experimental_wv_all_merged_v6.Robj")
load("dhh_wt_experimental_wv_all_merged_v6.Robj")
Idents(all_mouse_merged_final) <- "final.cell.type"
table(Idents(all_mouse_merged_final))
all_mouse_merged_final_sertoli <- subset(all_mouse_merged_final, idents = "Sertoli_cells")
VlnPlot(all_mouse_merged_final_sertoli, features = "Dhh", group.by = "treatment")

DimPlot(all_mouse_merged_final, group.by = "selection")

# 
# setwd("~/Desktop/Dhh_mouse")
# # save(all_mouse_merged_final, file = "all_mouse_merged_final_v1.Robj")
# # load("dhh_wt_experimental_wv_all_merged.Robj")
#   
# Idents(all_mouse_merged_final) <- "harmony_clusters"
# DimPlot(all_mouse_merged_final, label = T)
# table(Idents(all_mouse_merged_final))
# sertoli.subset <- subset(all_mouse_merged_final, idents = c(27, 42, 11, 34, 38, 33))
sertoli.subset <- subset(all_mouse_merged_final_sertoli)
sertoli.subset <- FindNeighbors(sertoli.subset, reduction = "integrated.harmony", dims = 1:30)
sertoli.subset <- FindClusters(sertoli.subset, resolution = 0.8, cluster.name = "harmony_clusters")
sertoli.subset <- RunUMAP(sertoli.subset, reduction = "integrated.harmony", dims = 1:30)
DimPlot(sertoli.subset, 
        # group.by = c("cell.type"),
        # split.by = "genotype", 
        label = T,
        raster = F)
FeaturePlot(sertoli.subset, 
            # features = c("Amh", "Sox9", "Wt1", "Gata1"), 
            features = c("Ddx4", "Vwf", "Cyp17a1", "Acta2", "Krt8"), 
            # split.by = "genotype",
            order =T
            
)



sertoli.subset2 <- sertoli.subset




DefaultAssay(sertoli.subset2) <- "RNA"

# keep only genes in your master list
genes.keep <- master.gene.list

sertoli.subset2 <- subset(
  sertoli.subset2,
  features = genes.keep
)


# Set identities to treatment
Idents(sertoli.subset2) <- "treatment"
# 
# table(Idents(sertoli.subset2))
sertoli.subset.subset <- subset(sertoli.subset2, idents = c("WT", "WWv", "Experimental_fertile", "Experimental_infertile"))
# sertoli.subset.subset <- sertoli.subset2
# sertoli.subset.subset <- RenameIdents(sertoli.subset.subset, 
#                                        'GFP transplanted' = "Transplanted_plus_germ_cells",
#                                        'WT' = "WT_plus_germ_cells",
#                                        'WWv' = "WWv_no_germ_cells",
#                                        'control' = "Transplanted_no_germ_cells")
# 
# sertoli.subset.subset$treatment <- Idents(sertoli.subset.subset)

DimPlot(sertoli.subset2)
Idents(sertoli.subset2) <- "treatment"
sertoli.subset2 <- RenameIdents(sertoli.subset2, 'WT' = 'fertile', "WWv" = "infertile",
                                "Experimental_fertile" = "fertile", "Experimental_infertile" = "infertile",
                                "DHHKO" = "infertile", "DHHenhancer" = "fertile", "STRA8KO" = "infertile")

sertoli.subset2$class <- Idents(sertoli.subset2)
DimPlot(sertoli.subset2, split.by = "class", group.by = "treatment",
        cols = c(
          
          "#937700", #WT
          "#545249", #WWv
          "#48b1e0", #"fertile",
          "#c05127", #infertile
          "#f26b68", #DHHKO
          "#426bb4", #DHHEnh
          "#58aa44" #stra8ko
          
        )
        
)

table(sertoli.subset2$treatment)


#make a volcano plot.



###############################
## Use RNA assay
###############################

DefaultAssay(sertoli.subset2) <- "RNA"

Idents(sertoli.subset2) <- sertoli.subset2$treatment

###############################
## Differential expression
###############################
sertoli.subset2 <- JoinLayers(sertoli.subset2)
markers <- FindMarkers(
  sertoli.subset2,
  ident.1 = "DHHKO",
  ident.2 = "WT",
  assay = "RNA"
)

markers$gene <- rownames(markers)

###############################
## Parameters
###############################

fc_cutoff <- 1
padj_cutoff <- 0.05

genes_to_label <- c(
  "Cldn3","Cldn10","Dab1","Chst2","Aif1","App","Jam2",
  "Sox9","Nfe2l2","Itga9","Pdgfa","Gdnf","Fgf11",
  "Wnt5a","Jag1","Csf1","Cxcl12"
)

###############################
## Prepare dataframe
###############################

markers <- markers %>%
  mutate(
    gene = rownames(.),
    
    # replace NA or 0 adjusted p-values
    p_val_adj_clean = ifelse(is.na(p_val_adj) | p_val_adj == 0,
                             .Machine$double.xmin,
                             p_val_adj),
    
    negLog10Padj = -log10(p_val_adj_clean),
    
    # cap extreme values (important)
    negLog10Padj_plot = pmin(negLog10Padj, 300),
    
    avg_log2FC_plot = pmax(pmin(avg_log2FC, 15), -15),
    
    # group = case_when(
    #   p_val_adj >= padj_cutoff ~ "NS",
    #   avg_log2FC >= fc_cutoff ~ "WT",
    #   avg_log2FC <= -fc_cutoff ~ "DHHKO",
    #   TRUE ~ "NS"
    group = case_when(
      p_val_adj >= padj_cutoff ~ "NS",
      avg_log2FC >= fc_cutoff ~ "DHHKO",
      avg_log2FC <= -fc_cutoff ~ "WT",
      TRUE ~ "NS"
    )
  )

label_df <- markers %>%
  dplyr::filter(
    gene %in% genes_to_label,
    group != "NS"
  )

###############################
## Volcano plot
###############################

ggplot(
  markers,
  aes(avg_log2FC_plot, negLog10Padj_plot)
) +
  
  geom_point(
    aes(fill = group),
    shape = 21,
    colour = "grey25",
    stroke = 0.1,
    size = 1.6,
    alpha = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      NS = "grey80",
      WT = "#937700",
      DHHKO = "#f26b68"
    )
  ) +
  
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    colour = "grey75",
    linetype = "dotted"
  ) +
  
  geom_hline(
    yintercept = -log10(padj_cutoff),
    colour = "grey75",
    linetype = "dotted"
  ) +
  
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    seed = 123
  ) +
  
  labs(
    x = "Log2 Fold Change (WT / DHHKO)",
    y = expression(-log[10]("adjusted "*italic(p))),
    title = "WT vs DHHKO Differential Expression"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none"
  )



# now for Stra8 KO


#make a volcano plot.



###############################
## Use RNA assay
###############################

DefaultAssay(sertoli.subset2) <- "RNA"

Idents(sertoli.subset2) <- sertoli.subset2$treatment

###############################
## Differential expression
###############################
sertoli.subset2 <- JoinLayers(sertoli.subset2)
markers <- FindMarkers(
  sertoli.subset2,
  ident.1 = "STRA8KO",
  ident.2 = "WT",
  assay = "RNA"
)

markers$gene <- rownames(markers)

###############################
## Parameters
###############################

fc_cutoff <- 1
padj_cutoff <- 0.05

genes_to_label <- c(
  "Cldn3","Cldn10","Dab1","Chst2","Aif1","App","Jam2",
  "Sox9","Nfe2l2","Itga9","Pdgfa","Gdnf","Fgf11",
  "Wnt5a","Jag1","Csf1","Cxcl12"
)

###############################
## Prepare dataframe
###############################

markers <- markers %>%
  mutate(
    gene = rownames(.),
    
    # replace NA or 0 adjusted p-values
    p_val_adj_clean = ifelse(is.na(p_val_adj) | p_val_adj == 0,
                             .Machine$double.xmin,
                             p_val_adj),
    
    negLog10Padj = -log10(p_val_adj_clean),
    
    # cap extreme values (important)
    negLog10Padj_plot = pmin(negLog10Padj, 300),
    
    avg_log2FC_plot = pmax(pmin(avg_log2FC, 8), -8),
    
    # group = case_when(
    #   p_val_adj >= padj_cutoff ~ "NS",
    #   avg_log2FC >= fc_cutoff ~ "WT",
    #   avg_log2FC <= -fc_cutoff ~ "DHHKO",
    #   TRUE ~ "NS"
    group = case_when(
      p_val_adj >= padj_cutoff ~ "NS",
      avg_log2FC >= fc_cutoff ~ "STRA8KO",
      avg_log2FC <= -fc_cutoff ~ "WT",
      TRUE ~ "NS"
    )
  )

label_df <- markers %>%
  dplyr::filter(
    gene %in% genes_to_label,
    group != "NS"
  )

###############################
## Volcano plot
###############################

ggplot(
  markers,
  aes(avg_log2FC_plot, negLog10Padj_plot)
) +
  
  geom_point(
    aes(fill = group),
    shape = 21,
    colour = "grey25",
    stroke = 0.1,
    size = 1.6,
    alpha = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      NS = "grey80",
      WT = "#937700",
      STRA8KO = "#58aa44"
    )
  ) +
  
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    colour = "grey75",
    linetype = "dotted"
  ) +
  
  geom_hline(
    yintercept = -log10(padj_cutoff),
    colour = "grey75",
    linetype = "dotted"
  ) +
  
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    seed = 123
  ) +
  
  labs(
    x = "Log2 Fold Change",
    y = "-log10(adjusted p)",
    title = "WT vs STRA8 Differential Expression"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none"
  )








table(sertoli.subset.subset$source)
Idents(sertoli.subset.subset) <- "source"
sertoli.subset.subset <- subset(sertoli.subset.subset, idents = "ThisPaper", invert = F)
Idents(sertoli.subset.subset) <- "treatment"
sertoli.subset.subset <- FindVariableFeatures(sertoli.subset.subset, nfeatures = 800)

# Average expression per treatment (pseudo-bulk)
avg.exp <- AverageExpression(
  sertoli.subset.subset,
  assays = "RNA",
  slot = "data",
  group.by = "treatment"
)$RNA

# Get variable genes
hvg <- VariableFeatures(sertoli.subset.subset)

# Keep only HVGs that are present in avg.exp
hvg <- intersect(hvg, rownames(avg.exp))

# Subset to variable genes
avg.exp.hvg <- avg.exp[hvg, ]

# Calculate Pearson correlation between treatments
cor.mat <- cor(
  as.matrix(avg.exp.hvg),
  method = "pearson"
)

# Heatmap
pheatmap(
  cor.mat,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  border_color = NA,
  display_numbers = TRUE,
  number_format = "%.3f",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = paste0(
    "Transcriptome Correlation Between Treatments\n(",
    length(hvg),
    " Variable Genes)"
  )
)

pheatmap(
  cor.mat,
  color = colorRampPalette(c("cyan", "yellow", "orange", "red"))(100),
  border_color = NA,
  display_numbers = TRUE,
  number_format = "%.3f",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = paste0(
    "Transcriptome Correlation Between Treatments\n(",
    length(hvg),
    " Variable Genes)"
  )
)
pheatmap(
  cor.mat,
  color = hcl.colors(100, "Inferno"),
  border_color = NA,
  display_numbers = TRUE,
  number_format = "%.3f",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = paste0(
    "Transcriptome Correlation Between Treatments\n(",
    length(hvg),
    " Variable Genes)"
  )
)




table(sertoli.subset2$treatment)

# from: https://journals.biologists.com/dev/article/139/23/4347/45398/Retinoic-acid-signaling-in-Sertoli-cells-regulates?guestAccessKey=
sertoli.table <- read_csv("/Users/ewhelan/Desktop/Dhh_mouse/dev080119-sup-supplemental_table1_modified_brief.csv")
sertoli.table <- sertoli.table[, 1:2]

oneThree <- sertoli.table$Gene[
  sertoli.table$Type == "oneThree"
]

fourSix <- sertoli.table$Gene[
  sertoli.table$Type == "fourSix"
]

sevenEight <- sertoli.table$Gene[
  sertoli.table$Type == "sevenEight"
]

nineTwelve <- sertoli.table$Gene[
  sertoli.table$Type == "nineTwelve"
]

length(oneThree)
length(fourSix)
length(sevenEight)
length(nineTwelve)


oneThree <- intersect(oneThree, rownames(sertoli.subset2))
fourSix <- intersect(fourSix, rownames(sertoli.subset2))
sevenEight <- intersect(sevenEight, rownames(sertoli.subset2))
nineTwelve <- intersect(nineTwelve, rownames(sertoli.subset2))

sapply(
  list(
    oneThree = oneThree,
    fourSix = fourSix,
    sevenEight = sevenEight,
    nineTwelve = nineTwelve
  ),
  length
)

sertoli.subset2 <- AddModuleScore(
  object = sertoli.subset2,
  features = list(
    oneThree,
    fourSix,
    sevenEight,
    nineTwelve
  ),
  name = c(
    "oneThree",
    "fourSix2",
    "sevenEight3",
    "nineTwelve4"
  )
)

FeaturePlot(
  sertoli.subset2,
  features = c(
    "oneThree1",
    "fourSix22",
    "sevenEight33",
    "nineTwelve44",
    "Sox9", "Dhh"
  ),
  order = T
) & DarkTheme() & scale_color_viridis(option = "C")

VlnPlot(sertoli.subset2, features = c(
  "oneThree1",
  "fourSix2",
  "sevenEight3",
  "nineTwelve4"), group.by = "harmony_clusters" )

# save(sertoli.subset2, file = "sertoli.subset2.Robj")







Idents(sertoli.subset2) <- "treatment"
table(Idents(sertoli.subset2))
sertoli.subset.wt <- subset(sertoli.subset2, idents = "WT")

Idents(sertoli.subset.wt) <- "source"
sertoli.subset.wt <- subset(sertoli.subset.wt, idents = "Guo", invert = T)

sertoli.subset.wt <- FindNeighbors(sertoli.subset.wt, reduction = "integrated.harmony", dims = 1:30)
sertoli.subset.wt <- FindClusters(sertoli.subset.wt, resolution = 0.6, cluster.name = "harmony_clusters")
sertoli.subset.wt <- RunUMAP(sertoli.subset.wt, reduction = "integrated.harmony", dims = 1:30)
DimPlot(sertoli.subset.wt, 
        # group.by = c("cell.type"),
        # split.by = "genotype", 
        label = T,
        raster = F)
FeaturePlot(sertoli.subset.wt, 
            # features = c("Amh", "Sox9", "Wt1", "Gata1"), 
            features = c("Ddx4", "Vwf", "Cyp17a1", "Acta2", "Krt8"), 
            # split.by = "genotype",
            order =T
            
)

sertoli.subset.wt <- AddModuleScore(
  object = sertoli.subset.wt,
  features = list(
    oneThree,
    fourSix,
    sevenEight,
    nineTwelve
  ),
  name = c(
    "oneThree",
    "fourSix",
    "sevenEight",
    "nineTwelve"
  )
)

FeaturePlot(
  sertoli.subset.wt,
  features = c(
    "oneThree1",
    "fourSix2",
    "sevenEight3",
    "nineTwelve4"
  ),
  order = T
) & DarkTheme() & scale_color_viridis(option = "C")


VlnPlot(sertoli.subset.wt)


#TRY TAKING THE TOP 10 FROM EACH LIST



DimPlot(sertoli.subset2, label = T)
FeaturePlot(sertoli.subset2, features = "Dhh", split.by = "treatment")

VlnPlot(sertoli.subset2, features = "Sox9", group.by = "treatment")






Idents(obj) <- "orig.ident"
table(Idents(obj))
obj2 <- RenameIdents(obj, 'DHHE7Ae' = "DHHE_rep1",
                     'DHHE7Au' = "DHHE_rep1",
                     'DHHE7Au2' = "DHHE_rep1",
                     'DHHKO15' = "DHHKO_rep1",
                     'DHHKO18' = "DHHKO_rep2",
                     'mouse.0' = "WT_rep1",
                     'mouse1' = "WT_rep2",
                     
                     'mouse2' = "WT_rep3",
                     'mouse3' = "WT_rep4",
                     'mouse4u' = "WT_rep5",
                     'mouse5' = "WT_rep6",
                     'mouse6' = "WT_rep7",
                     'mouse4e' = "WT_rep5",
                     'unselected.mouse.1' = "WT_rep8",
                     'EpCAM.mouse.1' = "WT_rep8",
                     'unselected.mouse.2' = "WT_rep9",
                     'EpCAM.mouse.2' = "WT_rep9",
                     
                     'ST1' = "WT_rep10",
                     'ST2' = "WT_rep11",
                     'ST3' = "WT_rep12",
                     'ST4' = "WT_rep13",
                     'ST5' = "WT_rep14",
                     'ST6' = "WT_rep15",
                     'ST7' = "WT_rep16",
                     'ST8' = "WT_rep17",
                     
                     'SPG1' = "WT_rep18",
                     'SPG2' = "WT_rep19",
                     'SPG3' = "WT_rep20",
                     
                     'INT1' = "WT_rep21",
                     'INT2' = "WT_rep22",
                     'INT3' = "WT_rep23",
                     'INT4' = "WT_rep24",
                     'INT5' = "WT_rep25",
                     'INT6' = "WT_rep26",
                     
                     'SER1' = "WT_rep27",
                     'SER2' = "WT_rep28",
                     'SER3' = "WT_rep29",
                     'SER4' = "WT_rep30",
                     'SER5' = "WT_rep31",
                     'SER6' = "WT_rep32",
                     'SER7' = "WT_rep33",
                     'SER8' = "WT_rep34",  
                     
                     
                     
                     
                     
                     "STRA8KO" = "STRA8KO",
                     
                     'WWv3GFPpositive' = "fertile_rep1",
                     'WWv3GFPnegative' = "fertile_rep1",
                     
                     'WWv3control' = "infertile_rep1",
                     'WWv5GFPpositive' = "fertile_rep2",
                     'WWv5GFPnegative'  = "fertile_rep2",
                     
                     'WWv5control' = "infertile_rep2",
                     'WWv6GFPpositive' = "fertile_rep3",
                     'WWv6GFPnegative' = "fertile_rep3", 
                     
                     'WWv6control' = "infertile_rep3",
                     'WWv7GFPnegative' = "fertile_rep4",
                     'WWv7control' = "infertile_rep4",
                     
                     'multiome1GFPpos' = "fertile_rep5",
                     'multiome1GFPnegA' = "fertile_rep5",
                     'multiome1GFPnegB' = "fertile_rep5" ,
                     
                     'multiome1CTL' = "infertile_rep5",
                     'multiome2GFPpos' = "fertile_rep6",
                     'multiome2GFPnegA' = "fertile_rep6",
                     'multiome2GFPnegB' = "fertile_rep6",
                     'multiome2CTL' = "infertile_rep6",
                     'multiome3CTLAMHR2posFACS' = "infertile_rep6",
                     
                     'multiome3transplantedAMHR2posFACS' = "fertile_rep7",
                     'multiome3CTLAMHR2negFACS' = "infertile_rep7",
                     'multiome3transplantedAMHR2negFACS'  = "fertile_rep7",
                     
                     'multiome4CTLAMHR2posMACS' = "infertile_rep8",
                     'multiome4transplantedAMHR2posMACS'  = "fertile_rep8",
                     'WWv_A' = "WWv_rep1",
                     'WWv_B' = "WWv_rep2"
                     
)
obj2$replicate <- Idents(obj2)
Idents(obj2) <- "orig.ident"
table(Idents(obj2))
obj2 <- RenameIdents(obj2, 'DHHE7Ae' = "Whole_cell",
                     'DHHE7Au' = "Whole_cell",
                     'DHHE7Au2' = "Whole_cell",
                     'DHHKO15' = "Whole_cell",
                     'DHHKO18' = "Whole_cell",
                     'mouse.0' = "Whole_cell",
                     "STRA8KO" = "Whole_cell",
                     'mouse1' = "Nuclear",
                     
                     'mouse2' = "Nuclear",
                     'mouse3' = "Nuclear",
                     'mouse4u' = "Nuclear",
                     'mouse5' = "Nuclear",
                     'mouse6' = "Nuclear",
                     'mouse4e' = "Nuclear",
                     'unselected.mouse.1' = "Whole_cell",
                     'EpCAM.mouse.1' = "Whole_cell",
                     'unselected.mouse.2' = "Whole_cell",
                     'EpCAM.mouse.2' = "Whole_cell",
                     
                     'ST1' = "Whole_cell",
                     'ST2' = "Whole_cell",
                     'ST3' = "Whole_cell",
                     'ST4' = "Whole_cell",
                     'ST5' = "Whole_cell",
                     'ST6' = "Whole_cell",
                     'ST7' = "Whole_cell",
                     'ST8' = "Whole_cell",
                     
                     'SPG1' = "Whole_cell",
                     'SPG2' = "Whole_cell",
                     'SPG3' = "Whole_cell",
                     
                     'INT1' = "Whole_cell",
                     'INT2' = "Whole_cell",
                     'INT3' = "Whole_cell",
                     'INT4' = "Whole_cell",
                     'INT5' = "Whole_cell",
                     'INT6' = "Whole_cell",
                     
                     'SER1' = "Whole_cell",
                     'SER2' = "Whole_cell",
                     'SER3' = "Whole_cell",
                     'SER4' = "Whole_cell",
                     'SER5' = "Whole_cell",
                     'SER6' = "Whole_cell",
                     'SER7' = "Whole_cell",
                     'SER8' = "Whole_cell",                     
                     
                     
                     'WWv3GFPpositive' = "Whole_cell",
                     'WWv3GFPnegative' = "Whole_cell",
                     
                     'WWv3control' = "Whole_cell",
                     'WWv5GFPpositive' = "Whole_cell",
                     'WWv5GFPnegative'  = "Whole_cell",
                     
                     'WWv5control' = "Whole_cell",
                     'WWv6GFPpositive' = "Whole_cell",
                     'WWv6GFPnegative' = "Whole_cell", 
                     
                     'WWv6control' = "Whole_cell",
                     'WWv7GFPnegative' = "Whole_cell",
                     'WWv7control' = "Whole_cell",
                     
                     'multiome1GFPpos' = "Nuclear",
                     'multiome1GFPnegA' = "Nuclear",
                     'multiome1GFPnegB' = "Nuclear" ,
                     
                     'multiome1CTL' = "Nuclear",
                     'multiome2GFPpos' = "Nuclear",
                     'multiome2GFPnegA' = "Nuclear",
                     'multiome2GFPnegB' = "Nuclear",
                     'multiome2CTL' = "Nuclear",
                     'multiome3CTLAMHR2posFACS' = "Nuclear",
                     
                     'multiome3transplantedAMHR2posFACS' = "Nuclear",
                     'multiome3CTLAMHR2negFACS' = "Nuclear",
                     'multiome3transplantedAMHR2negFACS'  = "Nuclear",
                     
                     'multiome4CTLAMHR2posMACS' = "Nuclear",
                     'multiome4transplantedAMHR2posMACS'  = "Nuclear",
                     'WWv_A' = "Nuclear",
                     'WWv_B' = "Nuclear"
                     
)

obj2$encapsulation <- Idents(obj2)
table(obj2$encapsulation)
table(obj2$replicate)

# --- subset Sertoli cells ---
sertoli.cells <-Cells(sertoli.subset2)


obj2.sertoli <- subset(obj2, cells = sertoli.cells)

table(obj2.sertoli$orig.ident)


# --- remove unwanted replicates ---
Idents(obj2.sertoli) <- "source"
table(obj2.sertoli$source)
obj2.sertoli <- subset(obj2.sertoli, idents = "Guo", invert = T)

Idents(obj2.sertoli) <- "orig.ident"
table(obj2.sertoli$orig.ident)
#remove singles
obj2.sertoli <- subset(obj2.sertoli, idents = c(
  "EpCAM.mouse.1",
  # "INT1", "INT2", "INT3",
  # "INT4",
  #                               "INT6",
  # "mouse2", "mouse3", "mouse4u","mouse5",
  "multiome3CTLAMHR2posFACS",
  "multiome2GFPpos"
  
  # "multiome3transplantedAMHR2negFACS",
  # "multiome3transplantedAMHR2posFACS",
  # "multiome4CTLAMHR2posMACS",
  # "multiome4transplantedAMHR2posMACS",
  # "multiome1GFPnegA",
  # "multiome1GFPnegB",
  
  
  # "SPG1",
  # "SPG2",
  # "SPG3",
  #                               "ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "ST7", "ST8"
), invert = TRUE)



# --- set identities ---
Idents(obj2.sertoli) <- "treatment"
DefaultAssay(obj2.sertoli) <- "RNA"
table(obj2.sertoli$treatment, obj2.sertoli$replicate)

# optional checks
VlnPlot(obj2.sertoli, features = "Dhh")



obj2.sertoli <- JoinLayers(obj2.sertoli)


avg_rep <- AverageExpression(
  obj2.sertoli,
  assays = "RNA",
  group.by = "replicate"
)

rep_vals <- avg_rep$RNA["Sox9", ]

df_rep <- data.frame(
  replicate = names(rep_vals),
  value = as.numeric(rep_vals)
)

rep_info <- obj2.sertoli@meta.data %>%
  group_by(replicate) %>%
  summarise(
    orig.ident = names(which.max(table(orig.ident))),
    .groups = "drop"
  )

df_rep <- df_rep %>%
  mutate(replicate = gsub("-rep", "_rep", replicate)) %>%
  left_join(rep_info, by = "replicate") %>%
  mutate(label = orig.ident)

# --- map replicate -> treatment ---
df_rep <- df_rep %>%
  mutate(treatment = case_when(
    grepl("^infertile", replicate) ~ "control",
    grepl("^fertile", replicate) ~ "GFP transplanted",
    grepl("^WT", replicate) ~ "WT",
    grepl("^DHHKO", replicate) ~ "DHHKO",
    grepl("^DHHE", replicate) ~ "DHHenhancer",
    grepl("^WWv", replicate) ~ "WWv",
    grepl("^STRA8KO", replicate) ~ "STRA8KO",
    TRUE ~ NA_character_
  )) %>%
  mutate(treatment = recode(
    treatment,
    "WT" = "WT_fertile",
    "GFP transplanted" = "Regenerated_fertile",
    "control" = "infertile",
    "DHHenhancer" = "DHHenhancer",
    "DHHKO" = "DHHKO",
    "WWv" = "WWv",
    "STRA8KO" = "STRA8KO"
  )) %>%
  mutate(treatment = factor(
    treatment,
    levels = c(
      "WT_fertile",
      "Regenerated_fertile",
      "infertile",
      "DHHenhancer",
      "DHHKO",
      "STRA8KO",
      "WWv"
    )
  ))

# --- compute mean + SEM ---
df_mean <- df_rep %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )



ggplot() +
  geom_col(
    data = df_mean,
    aes(x = treatment, y = mean),
    fill = "skyblue",
    color = "black"
  ) +
  geom_errorbar(
    data = df_mean,
    aes(
      x = treatment,
      ymin = mean - sem,
      ymax = mean + sem
    ),
    width = 0.2,
    size = 0.7
  ) +
  geom_jitter(
    data = df_rep,
    aes(x = treatment, y = value),
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  geom_text_repel(
    data = df_rep,
    aes(
      x = treatment,
      y = value,
      label = label
    ),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_classic() +
  labs(
    y = "Sox9 expression",
    x = "Treatment"
  )


### calculate proportions ###


DimPlot(somatic.integrated.germ4, group.by = "cell.type")
DimPlot(somatic.integrated.somatic6, group.by = "cell.type3")
somatic.integrated.somatic6$cell.type <- somatic.integrated.somatic6$cell.type3
table(somatic.integrated.germ4$replicate, somatic.integrated.germ4$cell.type, somatic.integrated.germ4$treatment, somatic.integrated.germ4$selection)
somatic.integrated.germ4$germline <- "germline"
somatic.integrated.somatic6$germline <- "somatic"

merge.somatic.germ <- merge(somatic.integrated.germ4, somatic.integrated.somatic6)
table(merge.somatic.germ$selection)
Idents(merge.somatic.germ) <- "selection"
merge.somatic.germ.GFPnegpos <- subset(merge.somatic.germ, idents = c("GFPnegative", "GFPpositive"))

table(merge.somatic.germ.GFPnegpos$cell.type, merge.somatic.germ.GFPnegpos$replicate, merge.somatic.germ.GFPnegpos$germline, merge.somatic.germ.GFPnegpos$treatment,  merge.somatic.germ.GFPnegpos$selection)
# write.csv(table.counts, file = "table.counts.csv")



# Desired biological ordering
cell_order <- c(
  # Germline
  "SSCs",
  "Progenitors",
  "DiffSpermatogonia",
  "PrelepSpermatocytes",
  "EarlySpermatocytes",
  "LateSpermatocytes",
  "RoundSpermatids",
  "ElongatingSpermatids",
  
  # Somatic
  "Sertoli",
  "Epithelial_cells",
  "Leydig_cells",
  "Peritubular_myoid",
  "Perivascular_smooth_muscle",
  "Mesenchymal_progenitors",
  "Endothelial",
  "Macrophages_peritubular",
  "Macrophages_interstitial",
  "Dendritic_cells",
  "T_cells_immunomodulatory",
  "T_cells_effector",
  "B_cells"
)

meta <- merge.somatic.germ.GFPnegpos@meta.data %>%
  as.data.frame() %>%
  dplyr::select(
    cell.type,
    replicate,
    germline,
    treatment,
    selection
  )

counts <- meta %>%
  group_by(
    cell.type,
    replicate,
    germline,
    treatment,
    selection
  ) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(!(selection == "GFPpositive" & treatment == "control")) %>%
  mutate(
    condition = case_when(
      selection == "GFPnegative" & treatment == "control" ~ "GFPneg_Control",
      selection == "GFPnegative" & treatment == "GFP transplanted" ~ "GFPneg_Transplant",
      selection == "GFPpositive" & treatment == "GFP transplanted" ~ "GFPpos_Transplant"
    ),
    condition = factor(
      condition,
      levels = c(
        "GFPneg_Control",
        "GFPneg_Transplant",
        "GFPpos_Transplant"
      )
    )
  ) %>%
  group_by(germline) %>%
  complete(
    cell.type,
    replicate,
    condition,
    fill = list(n_cells = 0)
  ) %>%
  ungroup()

means <- counts %>%
  group_by(cell.type, germline, condition) %>%
  summarise(
    mean_cells = mean(n_cells),
    .groups = "drop"
  )

means$cell.type <- factor(means$cell.type, levels = cell_order)
counts$cell.type <- factor(counts$cell.type, levels = cell_order)

ggplot(
  means,
  aes(
    x = cell.type,
    y = mean_cells,
    fill = condition
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_point(
    data = counts,
    aes(
      y = n_cells,
      fill = condition,
      shape = replicate
    ),
    colour = "black",
    stroke = 0.7,
    size = 3.2,
    position = position_jitterdodge(
      dodge.width = 0.8,
      jitter.width = 0.12
    )
  ) +
  facet_wrap(
    ~ germline,
    scales = "free_x"
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24, 25, 21)
  ) +
  scale_fill_manual(
    values = c(
      "GFPneg_Control" = "#999999",
      "GFPneg_Transplant" = "#56B4E9",
      "GFPpos_Transplant" = "#E69F00"
    )
  ) +
  labs(
    x = NULL,
    y = "Mean cells per replicate",
    fill = "Condition",
    shape = "Replicate"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )



### NOW WITH PROPORTIONS ###


cell_order <- c(
  "SSCs",
  "Progenitors",
  "DiffSpermatogonia",
  "PrelepSpermatocytes",
  "EarlySpermatocytes",
  "LateSpermatocytes",
  "RoundSpermatids",
  "ElongatingSpermatids",
  "Sertoli",
  "Epithelial_cells",
  "Leydig_cells",
  "Peritubular_myoid",
  "Perivascular_smooth_muscle",
  "Mesenchymal_progenitors",
  "Endothelial",
  "Macrophages_peritubular",
  "Macrophages_interstitial",
  "Dendritic_cells",
  "T_cells_immunomodulatory",
  "T_cells_effector",
  "B_cells"
)

meta <- merge.somatic.germ.GFPnegpos@meta.data %>%
  as.data.frame() %>%
  dplyr::select(
    cell.type,
    replicate,
    germline,
    treatment,
    selection
  )

counts <- meta %>%
  group_by(
    cell.type,
    replicate,
    germline,
    treatment,
    selection
  ) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(!(selection == "GFPpositive" & treatment == "control")) %>%
  mutate(
    condition = case_when(
      selection == "GFPnegative" & treatment == "control" ~ "GFPneg_Control",
      selection == "GFPnegative" & treatment == "GFP transplanted" ~ "GFPneg_Transplant",
      selection == "GFPpositive" & treatment == "GFP transplanted" ~ "GFPpos_Transplant"
    ),
    condition = factor(
      condition,
      levels = c(
        "GFPneg_Control",
        "GFPneg_Transplant",
        "GFPpos_Transplant"
      )
    )
  ) %>%
  group_by(germline) %>%
  complete(
    cell.type,
    replicate,
    condition,
    fill = list(n_cells = 0)
  ) %>%
  ungroup()

props <- counts %>%
  group_by(
    replicate,
    germline,
    condition
  ) %>%
  mutate(
    prop_cells = n_cells / sum(n_cells)
  ) %>%
  ungroup()

means <- props %>%
  group_by(
    cell.type,
    germline,
    condition
  ) %>%
  summarise(
    mean_prop = mean(prop_cells),
    .groups = "drop"
  )

means$cell.type <- factor(means$cell.type, levels = cell_order)
props$cell.type <- factor(props$cell.type, levels = cell_order)

ggplot(
  means,
  aes(
    x = cell.type,
    y = mean_prop,
    fill = condition
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_point(
    data = props,
    aes(
      y = prop_cells,
      fill = condition,
      shape = replicate
    ),
    colour = "black",
    stroke = 0.7,
    size = 3.2,
    position = position_jitterdodge(
      dodge.width = 0.8,
      jitter.width = 0.12
    )
  ) +
  facet_wrap(
    ~ germline,
    scales = "free_x"
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24, 25, 21)
  ) +
  scale_fill_manual(
    values = c(
      "GFPneg_Control" = "#999999",
      "GFPneg_Transplant" = "#56B4E9",
      "GFPpos_Transplant" = "#E69F00"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "Mean proportion of cells",
    fill = "Condition",
    shape = "Replicate"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )

### STACKED BAR PLOTS ###





meta <- merge.somatic.germ.GFPnegpos@meta.data %>%
  as.data.frame()

# ----------------------------
# Orders
# ----------------------------

germ_order <- c(
  "SSCs","Progenitors","DiffSpermatogonia","PrelepSpermatocytes",
  "EarlySpermatocytes","LateSpermatocytes","RoundSpermatids","ElongatingSpermatids"
)

somatic_order <- c(
  "Sertoli","Epithelial_cells","Leydig_cells","Peritubular_myoid",
  "Perivascular_smooth_muscle","Mesenchymal_progenitors","Endothelial",
  "Macrophages_peritubular","Macrophages_interstitial","Dendritic_cells",
  "T_cells_immunomodulatory","T_cells_effector","B_cells"
)

cell_order <- c(germ_order, somatic_order)

# ----------------------------
# Colours
# ----------------------------

germ_colors <- c(
  "SSCs"="chocolate4","Progenitors"="cornsilk4","DiffSpermatogonia"="darkseagreen3",
  "PrelepSpermatocytes"="plum4","EarlySpermatocytes"="paleturquoise",
  "LateSpermatocytes"="darkgreen","RoundSpermatids"="blue4","ElongatingSpermatids"="gold2"
)

som_colors <- c(
  "Sertoli"="#b12625","Epithelial_cells"="#88cdea","Leydig_cells"="#566b30",
  "Peritubular_myoid"="#d9a528","Perivascular_smooth_muscle"="#cd8163",
  "Mesenchymal_progenitors"="#88cdea","Endothelial"="#7f8080",
  "Macrophages_peritubular"="#d4a1ca","Macrophages_interstitial"="#7c57a4",
  "Dendritic_cells"="#6ec6a8","T_cells_immunomodulatory"="#9bcb3c",
  "T_cells_effector"="#ee6463","B_cells"="#5570b6"
)

cell_cols <- c(germ_colors, som_colors)

# ----------------------------
# Prepare data
# ----------------------------

df <- meta %>%
  filter(!(selection == "GFPpositive" & treatment == "control")) %>%
  mutate(
    condition = case_when(
      selection == "GFPnegative" & treatment == "control" ~ "GFPneg_Control",
      selection == "GFPnegative" & treatment == "GFP transplanted" ~ "GFPneg_Transplant",
      selection == "GFPpositive" & treatment == "GFP transplanted" ~ "GFPpos_Transplant"
    ),
    condition = factor(
      condition,
      levels = c("GFPneg_Control","GFPneg_Transplant","GFPpos_Transplant")
    ),
    cell.type = factor(cell.type, levels = cell_order)
  ) %>%
  dplyr::count(condition, cell.type) %>%
  complete(condition, cell.type, fill = list(n = 0)) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ----------------------------
# Alluvial plot (THIS is the correct geometry)
# ----------------------------

ggplot(df,
       aes(x = condition,
           stratum = cell.type,
           alluvium = cell.type,
           y = prop,
           fill = cell.type)) +
  
  geom_flow(alpha = 0.35, width = 0.6) +
  
  geom_stratum(
    width = 0.6,
    colour = "black",
    linewidth = 0.2
  ) +
  
  scale_fill_manual(values = cell_cols, drop = FALSE) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  
  labs(
    x = NULL,
    y = "Cell composition (%)",
    fill = "Cell type"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold")
  )




# now let's look at the replicates.

meta <- merge.somatic.germ.GFPnegpos@meta.data %>%
  as.data.frame()

# ----------------------------
# Define conditions
# ----------------------------

df <- meta %>%
  filter(!(selection == "GFPpositive" & treatment == "control")) %>%
  mutate(
    condition = case_when(
      selection == "GFPnegative" & treatment == "control" ~ "GFPneg_Control",
      selection == "GFPnegative" & treatment == "GFP transplanted" ~ "GFPneg_Transplant",
      selection == "GFPpositive" & treatment == "GFP transplanted" ~ "GFPpos_Transplant"
    ),
    condition = factor(
      condition,
      levels = c("GFPneg_Control","GFPneg_Transplant","GFPpos_Transplant")
    )
  )

# ----------------------------
# Germline subset
# ----------------------------

germ_cells <- c(
  "SSCs",
  "Progenitors",
  "DiffSpermatogonia",
  "PrelepSpermatocytes",
  "EarlySpermatocytes",
  "LateSpermatocytes",
  "RoundSpermatids",
  "ElongatingSpermatids"
)

df_germ <- df %>%
  filter(cell.type %in% germ_cells)

# ----------------------------
# Summarise per replicate
# ----------------------------

plot_df <- df_germ %>%
  group_by(condition, replicate) %>%
  summarise(n_germ = n(), .groups = "drop")

# ----------------------------
# Plot
# ----------------------------

ggplot(plot_df, aes(x = condition, y = n_germ)) +
  
  geom_col(
    aes(fill = condition),
    colour = "black",
    width = 0.7,
    alpha = 0.8
  ) +
  
  geom_point(
    size = 3,
    colour = "black",
    position = position_jitter(width = 0.12, height = 0)
  ) +
  
  scale_fill_manual(values = c(
    "GFPneg_Control" = "#4C78A8",
    "GFPneg_Transplant" = "#F58518",
    "GFPpos_Transplant" = "#54A24B"
  )) +
  
  labs(
    x = NULL,
    y = "Total germ cells per replicate",
    fill = "Condition"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold")
  )


table(somatic.integrated.germ4$selection, somatic.integrated.germ4$treatment)



meta <- merge.somatic.germ.GFPnegpos@meta.data %>%
  as.data.frame()

# ----------------------------
# Define conditions
# ----------------------------

df <- meta %>%
  filter(!(selection == "GFPpositive" & treatment == "control")) %>%
  mutate(
    condition = case_when(
      selection == "GFPnegative" & treatment == "control" ~ "GFPneg_Control",
      selection == "GFPnegative" & treatment == "GFP transplanted" ~ "GFPneg_Transplant",
      selection == "GFPpositive" & treatment == "GFP transplanted" ~ "GFPpos_Transplant"
    ),
    condition = factor(
      condition,
      levels = c("GFPneg_Control","GFPneg_Transplant","GFPpos_Transplant")
    )
  )

# ----------------------------
# Germline only
# ----------------------------

germ_cells <- c(
  "SSCs","Progenitors","DiffSpermatogonia","PrelepSpermatocytes",
  "EarlySpermatocytes","LateSpermatocytes","RoundSpermatids","ElongatingSpermatids"
)

df_germ <- df %>%
  filter(cell.type %in% germ_cells)

# ----------------------------
# Summarise per replicate
# ----------------------------

plot_df <- df_germ %>%
  group_by(condition, replicate) %>%
  summarise(n_germ = n(), .groups = "drop")

# ----------------------------
# Plot (simple faceted bars)
# ----------------------------

condition_cols <- c(
  "GFPneg_Control"    = "#c05127",  # blue
  "GFPneg_Transplant" = "#48b1e0",  # orange
  "GFPpos_Transplant" = "#54A24B"   # green
)

ggplot(plot_df, aes(x = replicate, y = n_germ, fill = condition)) +
  
  geom_col(
    colour = "black",
    width = 0.7
  ) +
  
  facet_wrap(~condition, nrow = 1) +
  
  scale_fill_manual(values = condition_cols) +
  
  labs(
    x = "Replicate",
    y = "Total germ cells"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"   # hide legend since facets already label conditions
  )



Idents(somatic.integrated.somatic6) <- "treatment"
table(Idents(somatic.integrated.somatic6))
somatic.integrated.somatic6.control <- subset(somatic.integrated.somatic6, idents = "control")
table(somatic.integrated.somatic6.control$selection)
Idents(somatic.integrated.somatic6.control) <- "selection"
somatic.integrated.somatic6.control <- subset(somatic.integrated.somatic6.control, idents = "GFPnegative")

Idents(somatic.integrated.somatic6.control) <- "replicate"
table(somatic.integrated.somatic6.control$replicate)
somatic.integrated.somatic6.control <- RenameIdents(somatic.integrated.somatic6.control, 
                                                    "rep3" = "high",
                                                    "rep7" = "high",
                                                    "rep8" = "high",
                                                    
                                                    "rep5" = "low",
                                                    "rep6" = "low",
                                                    "rep9" = "low"
)
somatic.integrated.somatic6.control$germ_cell_high_or_low <- Idents(somatic.integrated.somatic6.control)

DimPlot(somatic.integrated.somatic6.control)

Idents(somatic.integrated.somatic6.control) <- "cell.type"
VlnPlot(somatic.integrated.somatic6.control, features = "Sox9", split.by = "germ_cell_high_or_low")

table(Idents(somatic.integrated.somatic6.control))
DefaultAssay(somatic.integrated.somatic6.control) <- "RNA"

somatic.integrated.somatic6.control.sertoli <- subset(somatic.integrated.somatic6.control, idents = "Sertoli")
somatic.integrated.somatic6.control.myoid <- subset(somatic.integrated.somatic6.control, idents = "Peritubular_myoid")
somatic.integrated.somatic6.control.mesenchymal <- subset(somatic.integrated.somatic6.control, idents = "Mesenchymal_progenitors")


p1 <- VlnPlot(somatic.integrated.somatic6.control.sertoli, pt.size = 0,ncol = 5,
              features = c("Sox9", "Sparc", #"App", "Psap", 
                           "Dhh"), cols = c("#c05127", '#d8795b'),
              group.by = "germ_cell_high_or_low")
p2 <- VlnPlot(somatic.integrated.somatic6.control.myoid, pt.size = 0,ncol = 5,
              features = c("Egr1", "Cebpd", #"Jund", "Hspa8", 
                           "Igf1"),cols = c("#c05127", '#d8795b'),
              group.by = "germ_cell_high_or_low")
p3 <- VlnPlot(somatic.integrated.somatic6.control.mesenchymal, pt.size = 0,ncol = 5,
              features = c("Klf4", "Fosb", #"Hspa1a", "Igfbp4", 
                           "Col3a1"),cols = c("#c05127", '#d8795b'),
              group.by = "germ_cell_high_or_low")

plot_grid(p1,p2,p3, ncol = 1)

Idents(somatic.integrated.somatic6.control.sertoli) <- "germ_cell_high_or_low"
Idents(somatic.integrated.somatic6.control.myoid) <- "germ_cell_high_or_low"
Idents(somatic.integrated.somatic6.control.mesenchymal) <- "germ_cell_high_or_low"
specific.markers.sertoli <- FindMarkers(somatic.integrated.somatic6.control.sertoli, 
                                        ident.1 = "high", ident.2 = "low",
                                        features = c("Sox9", "Sparc", "App", "Psap", "Dhh"), 
                                        logfc.threshold = 0, min.pct = 0)
specific.markers.myoid <- FindMarkers(somatic.integrated.somatic6.control.myoid, 
                                      ident.1 = "high", ident.2 = "low", 
                                      features = c("Egr1", "Cebpd", "Jund", "Hspa8", "Igf1"), 
                                      logfc.threshold = 0, min.pct = 0)
specific.markers.mesenchymal <- FindMarkers(somatic.integrated.somatic6.control.mesenchymal, 
                                            ident.1 = "high", ident.2 = "low", 
                                            features = c("Klf4", "Fosb", "Hspa1a", "Igfbp4", "Col3a1"), 
                                            logfc.threshold = 0, min.pct = 0)
specific.markers.sertoli
specific.markers.myoid
specific.markers.mesenchymal









### okay now proportions with the all-merged dataset.


setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
load("dhh_wt_experimental_wv_all_merged_v6.Robj")

DimPlot(all_mouse_merged_final, group.by = "final.cell.type")
table(all_mouse_merged_final$final.cell.type, all_mouse_merged_final$treatment)


all_mouse_merged_final$treatment <- factor(all_mouse_merged_final$treatment, levels = c(
  "WT", "WWv", 
  "Experimental_fertile",
  "Experimental_infertile",
  "DHHKO",
  "DHHenhancer",
  "STRA8KO"
  
))

table(all_mouse_merged_final$final.cell.type, all_mouse_merged_final$treatment)


### ### ###



meta <- all_mouse_merged_final@meta.data %>%
  as.data.frame()

# ----------------------------
# Cell type order
# ----------------------------

germ_order <- c(
  "SSCs",
  "Progenitor_spermatogonia",
  "Diff.spermatogonia",
  "Preleptotene",
  "Early_spermatocytes",
  "Late_spermatocytes",
  "Round_spermatids",
  "Elongating_spermatids"
)

somatic_order <- c(
  "Sertoli_cells",
  "Rete_epithelial",
  "Leydig_cells",
  "Peritubular_myoid_cells",
  "Perivascular_smooth_muscle_cells",
  "Mesenchymal_progenitors",
  "Endothelial_cells",
  "Macrophages",
  "Dendritic",
  "T_cells"
)

cell_order <- c(germ_order, somatic_order)

# ----------------------------
# Colours
# ----------------------------

cell_cols <- c(
  
  # Germline
  "SSCs"                     = "chocolate4",
  "Progenitor_spermatogonia" = "cornsilk4",
  "Diff.spermatogonia"       = "darkseagreen3",
  "Preleptotene"             = "plum4",
  "Early_spermatocytes"      = "paleturquoise",
  "Late_spermatocytes"       = "darkgreen",
  "Round_spermatids"         = "blue4",
  "Elongating_spermatids"    = "gold2",
  
  # Somatic
  "Sertoli_cells"                    = "#b12625",
  "Rete_epithelial"                  = "#88cdea",
  "Leydig_cells"                     = "#566b30",
  "Peritubular_myoid_cells"          = "#d9a528",
  "Perivascular_smooth_muscle_cells" = "#cd8163",
  "Mesenchymal_progenitors"          = "#88cdea",
  "Endothelial_cells"                = "#7f8080",
  "Macrophages"                      = "#7c57a4",
  "Dendritic"                        = "#6ec6a8",
  "T_cells"                          = "#9bcb3c"
)

# ----------------------------
# Prepare data
# ----------------------------

df <- meta %>%
  mutate(
    final.cell.type = factor(
      final.cell.type,
      levels = cell_order
    )
  ) %>%
  dplyr::count(treatment, final.cell.type) %>%
  complete(
    treatment,
    final.cell.type,
    fill = list(n = 0)
  ) %>%
  group_by(treatment) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

# ----------------------------
# Plot
# ----------------------------

ggplot(
  df,
  aes(
    x = treatment,
    stratum = final.cell.type,
    alluvium = final.cell.type,
    y = prop,
    fill = final.cell.type
  )
) +
  
  geom_flow(
    alpha = 0.35,
    width = 0.6
  ) +
  
  geom_stratum(
    width = 0.6,
    colour = "black",
    linewidth = 0.2
  ) +
  
  scale_fill_manual(
    values = cell_cols,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  
  labs(
    x = NULL,
    y = "Cell composition (%)",
    fill = "Cell type"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    )
  )

#SUBSET ONLY THE undiff spermatogonia and tubular somatic cells (sertoli + myoid)

Idents(all_mouse_merged_final) <- "selection"
table(Idents(all_mouse_merged_final))
table(all_mouse_merged_final$selection, all_mouse_merged_final$treatment)
all_mouse_merged_unselected <- subset(all_mouse_merged_final, idents = "Unselected")
meta <- all_mouse_merged_unselected@meta.data %>%
  as.data.frame()

# ----------------------------
# Treatments to include
# ----------------------------

treatment_order <- c(
  "WT",
  "WWv",  "Experimental_infertile",
  "Experimental_fertile"
  
)

# ----------------------------
# Cell types to include
# ----------------------------

cell_order <- c(
  "SSCs",
  "Progenitor_spermatogonia",
  "Peritubular_myoid_cells",
  "Sertoli_cells"
)

# ----------------------------
# Colours
# ----------------------------

cell_cols <- c(
  "SSCs"                     = "chocolate4",
  "Progenitor_spermatogonia" = "cornsilk4",
  "Peritubular_myoid_cells"  = "#d9a528",
  "Sertoli_cells"            = "#b12625"
)

# ----------------------------
# Prepare data
# ----------------------------

df <- meta %>%
  dplyr::filter(
    treatment %in% treatment_order,
    final.cell.type %in% cell_order
  ) %>%
  dplyr::mutate(
    treatment = factor(
      treatment,
      levels = treatment_order
    ),
    final.cell.type = factor(
      final.cell.type,
      levels = cell_order
    )
  ) %>%
  dplyr::count(treatment, final.cell.type) %>%
  complete(
    treatment,
    final.cell.type,
    fill = list(n = 0)
  ) %>%
  group_by(treatment) %>%
  dplyr::mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

# ----------------------------
# Plot
# ----------------------------

ggplot(
  df,
  aes(
    x = treatment,
    stratum = final.cell.type,
    alluvium = final.cell.type,
    y = prop,
    fill = final.cell.type
  )
) +
  
  geom_flow(
    alpha = 0.4,
    width = 0.6
  ) +
  
  geom_stratum(
    width = 0.6,
    colour = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_manual(
    values = cell_cols,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  
  labs(
    x = NULL,
    y = "Cell composition (%)",
    fill = "Cell type"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    )
  )

### now statistics:



ratio_df <- all_mouse_merged_unselected@meta.data %>%
  as.data.frame() %>%
  filter(
    treatment %in% c(
      "WT",
      "WWv",
      "Experimental_fertile",
      "Experimental_infertile"
    )
  ) %>%
  mutate(
    group = case_when(
      final.cell.type %in% c(
        "SSCs",
        "Progenitor_spermatogonia"
      ) ~ "undiff",
      
      final.cell.type %in% c(
        "Sertoli_cells",
        "Peritubular_myoid_cells"
      ) ~ "tubular",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  dplyr::count(treatment, orig.ident, group) %>%
  pivot_wider(
    names_from = group,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    ratio = undiff / tubular
  )

ratio_df

ggplot(
  ratio_df,
  aes(treatment, ratio, colour = treatment)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  theme_classic() +
  ylab("Undifferentiated / Tubular ratio")

#WT vs WWV
wilcox.test(
  ratio ~ treatment,
  data = ratio_df %>%
    filter(treatment %in% c("WT", "WWv"))
)

# WWV vs experimental infertile
wilcox.test(
  ratio ~ treatment,
  data = ratio_df %>%
    filter(treatment %in% c("WWv", "Experimental_infertile"))
)





### NOW ANALYSIS OF DHHKO AND STRA8KO ###

setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
# save(all_mouse_merged_final, file = "dhh_wt_experimental_wv_all_merged_v6.Robj")
load("dhh_wt_experimental_wv_all_merged_v6.Robj")

# need to make sure the gene list is the same across the samples
# made "master.gene.list" earlier

DefaultAssay(all_mouse_merged_final) <- "RNA"

all_mouse_merged_final_subset <- subset(
  all_mouse_merged_final,
  features = intersect(master.gene.list, rownames(all_mouse_merged_final))
)

# save(all_mouse_merged_final_subset, file="all_mouse_merged_final_subset_common_genes.Robj")
load("all_mouse_merged_final_subset_common_genes.Robj")
master.gene.list<-rownames(all_mouse_merged_final_subset)

DimPlot(all_mouse_merged_final_subset, group.by = "treatment")
Idents(all_mouse_merged_final_subset) <- "treatment"
stra8ko.analysis <- subset(all_mouse_merged_final_subset, idents = c("WT", "STRA8KO"))
dhhko.analysis <- subset(all_mouse_merged_final_subset, idents = c("WT", "DHHKO"))

DefaultAssay(stra8ko.analysis) <- "RNA"
DefaultAssay(dhhko.analysis) <- "RNA"
Idents(stra8ko.analysis) <- "final.cell.type" 
Idents(dhhko.analysis) <- "final.cell.type" 
table(Idents(stra8ko.analysis))

dhhko.analysis.sertoli <- subset(dhhko.analysis, idents = "Sertoli_cells")
dhhko.analysis.myoid <- subset(dhhko.analysis, idents = "Peritubular_myoid_cells")
dhhko.analysis.mesenchymal <- subset(dhhko.analysis, idents = "Mesenchymal_progenitors")

stra8ko.analysis.sertoli <- subset(stra8ko.analysis, idents = "Sertoli_cells")
stra8ko.analysis.myoid <- subset(stra8ko.analysis, idents = "Peritubular_myoid_cells")
stra8ko.analysis.mesenchymal <- subset(stra8ko.analysis, idents = "Mesenchymal_progenitors")

sertoli.genes <- intersect(c("Sox9", "Sparc", "App", "Psap", 
                             "Dhh"), master.gene.list)
myoid.genes <- intersect(c("Egr1", "Cebpd", "Jund", "Hspa8", 
                           "Igf1"), master.gene.list)
mesench.genes <- intersect(c("Klf4", "Fosb", "Hspa1a", "Igfbp4", 
                             "Col3a1"), master.gene.list)

p1 <- VlnPlot(stra8ko.analysis.sertoli, pt.size = 0,ncol = 5,
              features = sertoli.genes, cols = c("#c05127", "#48b1e0"),
              group.by = "treatment")
p2 <- VlnPlot(stra8ko.analysis.myoid, pt.size = 0,ncol = 5,
              features = myoid.genes ,cols = c("#c05127", "#48b1e0"),
              group.by = "treatment")
p3 <- VlnPlot(stra8ko.analysis.mesenchymal, pt.size = 0,ncol = 5,
              features = mesench.genes,cols = c("#c05127", "#48b1e0"),
              group.by = "treatment")

plot_grid(p1,p2,p3, ncol = 1)

p1 <- VlnPlot(dhhko.analysis.sertoli, pt.size = 0,ncol = 5,
              features = sertoli.genes, cols = c("#c05127", "#48b1e0"),
              group.by = "treatment")
p2 <- VlnPlot(dhhko.analysis.myoid, pt.size = 0,ncol = 5,
              features = myoid.genes,cols = c("#c05127", "#48b1e0"),
              group.by = "treatment")
p3 <- VlnPlot(dhhko.analysis.mesenchymal, pt.size = 0,ncol = 5,
              features = mesench.genes,cols = c("#c05127", "#48b1e0"),
              group.by = "treatment")

plot_grid(p1,p2,p3, ncol = 1)

#may need to create this object further down:
spermatid.blacklist <- spermatid.genes

# IPA of just the Sertoli cells
dhhko.analysis.sertoli <- JoinLayers(dhhko.analysis.sertoli)
stra8ko.analysis.sertoli <- JoinLayers(stra8ko.analysis.sertoli)

Idents(dhhko.analysis.sertoli) <- "treatment"
Idents(stra8ko.analysis.sertoli) <- "treatment"

markers.dhhko.analysis.sertoli <- FindMarkers(
  dhhko.analysis.sertoli,
  ident.1 = "DHHKO",
  ident.2 = "WT"
)

markers.stra8.analysis.sertoli <- FindMarkers(
  stra8ko.analysis.sertoli,
  ident.1 = "STRA8KO",
  ident.2 = "WT"
)

# Remove spermatid genes
markers.dhhko.analysis.sertoli <- markers.dhhko.analysis.sertoli[
  !rownames(markers.dhhko.analysis.sertoli) %in% spermatid.blacklist,
]

markers.stra8.analysis.sertoli <- markers.stra8.analysis.sertoli[
  !rownames(markers.stra8.analysis.sertoli) %in% spermatid.blacklist,
]

setwd("~/Desktop/Dhh_mouse/DHHKO_filtered_markers")

write.csv(
  markers.dhhko.analysis.sertoli,
  file = "markers.dhhko.analysis.sertoli_filtered.csv"
)

write.csv(
  markers.stra8.analysis.sertoli,
  file = "markers.stra8.analysis.sertoli_filtered.csv"
)



table(Idents(dhhko.analysis))

dhhko.analysis.joined <-JoinLayers(dhhko.analysis)

dhhko.analysis.SSCs <- subset(dhhko.analysis.joined, idents = "SSCs")
dhhko.analysis.Progenitor_spermatogonia <- subset(dhhko.analysis.joined, idents = "Progenitor_spermatogonia")
dhhko.analysis.Diff.spermatogonia <- subset(dhhko.analysis.joined, idents = "Diff.spermatogonia")
dhhko.analysis.Preleptotene <- subset(dhhko.analysis.joined, idents = "Preleptotene")
dhhko.analysis.Early_spermatocytes <- subset(dhhko.analysis.joined, idents = "Early_spermatocytes")
dhhko.analysis.Late_spermatocytes <- subset(dhhko.analysis.joined, idents = "Late_spermatocytes")
dhhko.analysis.Round_spermatids <- subset(dhhko.analysis.joined, idents = "Round_spermatids")
dhhko.analysis.Elongating_spermatids <- subset(dhhko.analysis.joined, idents = "Elongating_spermatids")


Idents(dhhko.analysis.SSCs) <- "treatment"
Idents(dhhko.analysis.Progenitor_spermatogonia) <- "treatment"
Idents(dhhko.analysis.Diff.spermatogonia) <- "treatment"
Idents(dhhko.analysis.Preleptotene) <- "treatment"
Idents(dhhko.analysis.Early_spermatocytes) <- "treatment"
Idents(dhhko.analysis.Late_spermatocytes) <- "treatment"
Idents(dhhko.analysis.Round_spermatids) <- "treatment"
Idents(dhhko.analysis.Elongating_spermatids) <- "treatment"


markers.dhhko.analysis.SSCs <- FindMarkers(dhhko.analysis.SSCs, 
                                           ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Progenitor_spermatogonia <- FindMarkers(dhhko.analysis.Progenitor_spermatogonia, 
                                                               ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Diff.spermatogonia <- FindMarkers(dhhko.analysis.Diff.spermatogonia, 
                                                         ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Preleptotene <- FindMarkers(dhhko.analysis.Preleptotene, 
                                                   ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Early_spermatocytes <- FindMarkers(dhhko.analysis.Early_spermatocytes, 
                                                          ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Late_spermatocytes <- FindMarkers(dhhko.analysis.Late_spermatocytes, 
                                                         ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Round_spermatids <- FindMarkers(dhhko.analysis.Round_spermatids, 
                                                       ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Elongating_spermatids <- FindMarkers(dhhko.analysis.Elongating_spermatids, 
                                                            ident.1 = "DHHKO", ident.2 = "WT")


all.markers.dhhko.analysis <- FindAllMarkers(dhhko.analysis.joined)
all.markers.dhhko.analysis$pct.diff <- all.markers.dhhko.analysis$pct.1 - all.markers.dhhko.analysis$pct.2

# spermatid populations
spermatid.clusters <- c("Round_spermatids", "Elongating_spermatids")

# all other populations
other.clusters <- setdiff(
  unique(all.markers.dhhko.analysis$cluster),
  spermatid.clusters
)

markers.filt <- subset(
  all.markers.dhhko.analysis,
  avg_log2FC > 2 &
    pct.diff > 0.20 &
    p_val_adj < 0.05
)

# genes appearing as markers in spermatids
spermatid.genes <- unique(
  markers.filt$gene[
    markers.filt$cluster %in% spermatid.clusters
  ]
)
length(spermatid.genes)





markers.dhhko.analysis.SSCs$pct.diff <- markers.dhhko.analysis.SSCs$pct.1 - markers.dhhko.analysis.SSCs$pct.2
markers.dhhko.analysis.Progenitor_spermatogonia$pct.diff <- markers.dhhko.analysis.Progenitor_spermatogonia$pct.1 - markers.dhhko.analysis.Progenitor_spermatogonia$pct.2
markers.dhhko.analysis.Diff.spermatogonia$pct.diff <- markers.dhhko.analysis.Diff.spermatogonia$pct.1 - markers.dhhko.analysis.Diff.spermatogonia$pct.2
markers.dhhko.analysis.Preleptotene$pct.diff <- markers.dhhko.analysis.Preleptotene$pct.1 - markers.dhhko.analysis.Preleptotene$pct.2
markers.dhhko.analysis.Early_spermatocytes$pct.diff <- markers.dhhko.analysis.Early_spermatocytes$pct.1 - markers.dhhko.analysis.Early_spermatocytes$pct.2
markers.dhhko.analysis.Late_spermatocytes$pct.diff <- markers.dhhko.analysis.Late_spermatocytes$pct.1 - markers.dhhko.analysis.Late_spermatocytes$pct.2
markers.dhhko.analysis.Round_spermatids$pct.diff <- markers.dhhko.analysis.Round_spermatids$pct.1 - markers.dhhko.analysis.Round_spermatids$pct.2
markers.dhhko.analysis.Elongating_spermatids$pct.diff <- markers.dhhko.analysis.Elongating_spermatids$pct.1 - markers.dhhko.analysis.Elongating_spermatids$pct.2

setwd("~/Desktop/Dhh_mouse")
spermatid.blacklist <- spermatid.genes

# List of marker tables
marker.list <- list(
  SSCs = markers.dhhko.analysis.SSCs,
  Progenitor_spermatogonia = markers.dhhko.analysis.Progenitor_spermatogonia,
  Diff_spermatogonia = markers.dhhko.analysis.Diff.spermatogonia,
  Preleptotene = markers.dhhko.analysis.Preleptotene,
  Early_spermatocytes = markers.dhhko.analysis.Early_spermatocytes,
  Late_spermatocytes = markers.dhhko.analysis.Late_spermatocytes,
  Round_spermatids = markers.dhhko.analysis.Round_spermatids,
  Elongating_spermatids = markers.dhhko.analysis.Elongating_spermatids
)


filtered.markers <- list()

for (nm in names(marker.list)) {
  
  df <- marker.list[[nm]]
  
  # remove spermatid genes only from non-spermatid stages
  if (!nm %in% c("Round_spermatids", "Elongating_spermatids")) {
    df <- df[!(rownames(df) %in% spermatid.blacklist), ]
  }
  
  filtered.markers[[nm]] <- df
}

outdir <- "DHHKO_filtered_markers"
dir.create(outdir, showWarnings = FALSE)

for (nm in names(filtered.markers)) {
  
  write.csv(
    filtered.markers[[nm]],
    file = file.path(outdir, paste0(nm, "_filtered.csv")),
    row.names = TRUE
  )
}

# # filtered versions of spermatids too
# round.filtered <- markers.dhhko.analysis.Round_spermatids[
#   !(rownames(markers.dhhko.analysis.Round_spermatids) %in% spermatid.blacklist),
# ]
# 
# elong.filtered <- markers.dhhko.analysis.Elongating_spermatids[
#   !(rownames(markers.dhhko.analysis.Elongating_spermatids) %in% spermatid.blacklist),
# ]
# 
# write.csv(round.filtered,
#           "DHHKO_filtered_markers/Round_spermatids_filtered.csv")
# 
# write.csv(elong.filtered,
#           "DHHKO_filtered_markers/Elongating_spermatids_filtered.csv")


FeaturePlot(dhhko.analysis, features = "Tonsl", split.by = "treatment")
VlnPlot(dhhko.analysis, features = "Prps1l1", split.by = "treatment")



marker_list <- list(
  SSCs = markers.dhhko.analysis.SSCs,
  Progenitor_spermatogonia = markers.dhhko.analysis.Progenitor_spermatogonia,
  Diff_spermatogonia = markers.dhhko.analysis.Diff.spermatogonia,
  Preleptotene = markers.dhhko.analysis.Preleptotene,
  Early_spermatocytes = markers.dhhko.analysis.Early_spermatocytes,
  Late_spermatocytes = markers.dhhko.analysis.Late_spermatocytes,
  Round_spermatids = markers.dhhko.analysis.Round_spermatids,
  Elongating_spermatids = markers.dhhko.analysis.Elongating_spermatids
)

gene_counts <- lapply(names(marker_list), function(ct){
  
  df <- marker_list[[ct]]
  
  data.frame(
    CellType = ct,
    Up = sum(df$avg_log2FC >= 1 & df$p_val_adj < 0.05, na.rm = TRUE),
    Down = sum(df$avg_log2FC <= -1 & df$p_val_adj < 0.05, na.rm = TRUE)
  )
  
}) %>% bind_rows()

plot_df <- gene_counts %>%
  mutate(Down = -Down) %>%
  pivot_longer(
    cols = c(Up, Down),
    names_to = "Direction",
    values_to = "Genes"
  )

ggplot(plot_df,
       aes(x = Genes,
           y = factor(CellType, levels = rev(gene_counts$CellType)),
           fill = Direction)) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = c(
    "Up" = "#D55E00",
    "Down" = "#0072B2"
  )) +
  labs(
    x = "Number of DE genes",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14)

#make a dotplot of all treatments.

genes <- c("Dhh", "Gdnf", "Sox9", "Cldn10", 
           "Jag1", "Pdgfa", "Itga9")
DotPlot(
  object = all_mouse_merged_final_subset, 
  # group.by = "orig.ident", 
  features = rev((genes)), 
  cols = c("gray", "#890600")
) + 
  theme(
    axis.text.x = element_text(
      angle = -45,          # rotate labels 30 degrees
      vjust = 1, 
      hjust = 1, 
      size = 12
    ),
    axis.text.y = element_text(
      size = 12, 
      face = "italic"
    ),                  
    axis.title = element_text(size = 12)
  ) +   
  # scale_color_scico(palette = "vikO") + 
  scale_size(range = c(-5, 10), limits = c(0, 100)) +
  scale_y_discrete(position = "right") + 
  coord_flip() +
  ggtitle("in vivo")


## that didn't work so let's just do IPA now with just Sertoli cells


dhhko.analysis.SSCs <- subset(dhhko.analysis.joined, idents = "SSCs")
dhhko.analysis.Progenitor_spermatogonia <- subset(dhhko.analysis.joined, idents = "Progenitor_spermatogonia")
dhhko.analysis.Diff.spermatogonia <- subset(dhhko.analysis.joined, idents = "Diff.spermatogonia")
dhhko.analysis.Preleptotene <- subset(dhhko.analysis.joined, idents = "Preleptotene")
dhhko.anal

Idents(dhhko.analysis.SSCs) <- "treatment"
Idents(dhhko.analysis.Progenitor_spermatogonia) <- "treatment"
Idents(dhhko.analysis.Diff.spermatogonia) <- "treatment"
Idents(dhhko.analysis.Preleptotene) <- "treatment"
Idents(dhhko.analysis.Early_spermatocytes) <- "treatment"
Idents(dhhko.analysis.Late_spermatocytes) <- "treatment"
Idents(dhhko.analysis.Round_spermatids) <- "treatment"
Idents(dhhko.analysis.Elongating_spermatids) <- "treatment"


markers.dhhko.analysis.SSCs <- FindMarkers(dhhko.analysis.SSCs, 
                                           ident.1 = "DHHKO", ident.2 = "WT")
markers.dhhko.analysis.Progenitor_spermatogonia <- FindMarkers(dhhko.analysis.Progenitor_spermatogonia, 
                                                               ident.1 = "DHHKO", ident.2 = "WT")


## look if Wv are mature sertoli ##

Idents(all_mouse_merged_final) <- "final.cell.type"
DimPlot(all_mouse_merged_final)
sertoli.subset <- subset(all_mouse_merged_final, idents = "Sertoli_cells")
Idents(sertoli.subset) <- "treatment"
sertoli.subset <- subset(sertoli.subset, idents = c("DHHKO", "DHHenhancer", "STRA8KO"), invert = T)
VlnPlot(sertoli.subset, cols = c(
  "#937700", "#c9c9c9",
  "#426bb4", "#D55E00"
),
group.by = "treatment",pt.size = 0,
features = c("Amh", "Cldn11"))



## additional spermatogonia analysis ##



setwd("/Users/ewhelan/Desktop/Somatic Project/somatic_samples")
load("somatic.integrated.germ4.Robj")

Idents(somatic.integrated.germ4) <- "cell.type"
DimPlot(somatic.integrated.germ4)

spermatogonia.subset <- subset(somatic.integrated.germ4, idents = c("SSCs", "Progenitors" #, "DiffSpermatogonia"
))

# spermatogonia.subset <- FindVariableFeatures(spermatogonia.subset)
# spermatogonia.subset <- ScaleData(spermatogonia.subset)
# spermatogonia.subset <- RunPCA(spermatogonia.subset, dims = 1:30)
DefaultAssay(spermatogonia.subset) <- "integrated"
spermatogonia.subset <- FindNeighbors(spermatogonia.subset, dims = 1:30)
spermatogonia.subset <- FindClusters(spermatogonia.subset, resolution = 0.2)
spermatogonia.subset <- RunUMAP(spermatogonia.subset, dims = 1:30)
DimPlot(spermatogonia.subset, group.by = "cell.type", cols = c("chocolate4",
                                                               "cornsilk4"))
DefaultAssay(spermatogonia.subset) <- "RNA"
FeaturePlot(spermatogonia.subset, 
            features = c("Gfra1", "Etv5", 
                         "Id4", "Ret",# "Bcl6b", #"Kit",
                         "Utf1", #"Sall4",  #"Zbtb16", 
                         #"Foxo1",
                         "Lin28a", 
                         "Neurog3", 
                         "Rarg", #"Sohlh2"
                         "Mki67"
            ), ncol = 2,
            
            order = T) & DarkTheme() & scale_color_viridis(option = "C")




#banana
Idents(all_mouse_merged_final_subset) <- "final.cell.type"
DimPlot(all_mouse_merged_final_subset)
table(Idents(all_mouse_merged_final_subset))
all_mouse_merged_GCs <- subset(all_mouse_merged_final_subset, idents = c(
  "SSCs",
  "Progenitor_spermatogonia",
  "Diff.spermatogonia",
  "Preleptotene",
  "Early_spermatocytes",
  "Late_spermatocytes",
  "Round_spermatids",
  "Elongating_spermatids"
))
Idents(all_mouse_merged_GCs) <- "treatment"
table(Idents(all_mouse_merged_GCs))
all_mouse_merged_GCs_DHHKO <- subset(all_mouse_merged_GCs, idents = c("WT", "DHHKO"))
Idents(all_mouse_merged_GCs_DHHKO) <- "final.cell.type"


genes <- c(
  "Dcn", "Pou5f1", "Apoe", "Igf1", "Gfra1","Morc1",
  "Id4","Sparc", "Stra8", "Sohlh2", "Hormad1", "Wt1",
  "Pou3f2", "Fgfr4", "Rhox8", "Ybx2", "Mnd1", "Syce3",
  "Fox", "Jund", "Pax8", "Runx1", "Eomes", "Tex101",
  "Rhox10", "Dnd1", "Meig1","Nanog", "Ldhc",
  "Tex15", "Sycp1", "Sycp2", "Rad51ap2", "Spata16",
  "Topbp1", "Cdk1",
  "Dazl", "Bcl6b","Utf1", "Etv5", "Aurkb", "Aurkc",
  "Spata33", "Odf1", "Tekt1", "Prm1"
)

obj <- JoinLayers(all_mouse_merged_GCs_DHHKO)

# metadata
meta <- obj@meta.data

# cell types
celltypes <- unique(meta$final.cell.type)

# calculate WT vs DHHKO log2FC per cell type
fc_list <- lapply(celltypes, function(ct){
  
  cells_ct <- rownames(meta[meta$final.cell.type == ct, ])
  
  sub <- subset(obj, cells = cells_ct)
  
  Idents(sub) <- "treatment"
  
  markers <- FindMarkers(
    sub,
    ident.1 = "DHHKO",
    ident.2 = "WT",
    features = genes,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  markers$gene <- rownames(markers)
  markers$cell.type <- ct
  
  markers[,c("gene","cell.type","avg_log2FC")]
})

fc_df <- bind_rows(fc_list)


# add expression detection (% expressed) for dot size
expr <- FetchData(
  obj,
  vars = c(genes,"final.cell.type")
)

pct_df <- expr %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(final.cell.type, gene) %>%
  summarise(
    pct.exp = mean(expression > 0)*100,
    .groups="drop"
  )

plot_df <- left_join(
  fc_df,
  pct_df,
  by=c("cell.type"="final.cell.type","gene")
)


genes <- c(
  "Dcn", "Pou5f1", "Apoe", "Igf1", "Gfra1","Morc1",
  "Id4","Sparc", "Stra8", "Sohlh2", "Hormad1", "Wt1",
  "Pou3f2", "Fgfr4", "Rhox8", "Ybx2", "Mnd1", "Syce3",
  "Fox", "Jund", "Pax8", "Runx1", "Eomes", "Tex101",
  "Rhox10", "Dnd1", "Meig1","Nanog", "Ldhc",
  "Tex15", "Sycp1", "Sycp2", "Rad51ap2", "Spata16",
  "Topbp1", "Cdk1",
  "Dazl", "Bcl6b","Utf1", "Etv5", "Aurkb", "Aurkc",
  "Spata33", "Odf1", "Tekt1"
)

obj <- JoinLayers(all_mouse_merged_GCs_DHHKO)

# metadata
meta <- obj@meta.data

# cell types
celltypes <- unique(meta$final.cell.type)

genes_use <- genes[genes %in% rownames(obj)]


fc_list <- lapply(unique(obj$final.cell.type), function(ct){
  
  message("Processing ", ct)
  
  sub <- subset(
    obj,
    subset = final.cell.type == ct
  )
  
  Idents(sub) <- "treatment"
  
  markers <- FindMarkers(
    sub,
    ident.1 = "DHHKO",
    ident.2 = "WT",
    features = genes_use,
    logfc.threshold = 0,
    min.pct = 0,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox"
  )
  
  markers$gene <- rownames(markers)
  markers$cell.type <- ct
  
  markers[, c("gene","cell.type","avg_log2FC")]
})

fc_df <- bind_rows(fc_list)




# add expression detection (% expressed) for dot size
expr <- FetchData(
  obj,
  vars = c(genes_use, "final.cell.type")
)

pct_df <- expr %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  group_by(final.cell.type, gene) %>%
  summarise(
    pct.exp = mean(expression > 0) * 100,
    .groups = "drop"
  )

plot_df <- left_join(
  fc_df,
  pct_df,
  by = c("cell.type" = "final.cell.type",
         "gene")
)

# preserve gene ordering (only genes that exist)
plot_df$gene <- factor(
  plot_df$gene,
  levels = genes_use
)

# plot
ggplot(plot_df,
       aes(
         x=gene,
         y=cell.type,
         size=pct.exp,
         color=avg_log2FC
       )) +
  geom_point() +
  scale_color_gradient2(
    low="blue",
    mid="white",
    high="red",
    midpoint=0,
    name="DHHKO vs WT\nlog2FC"
  ) +
  scale_size(
    range=c(1,8),
    name="% expressed"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle=90,
      hjust=1,
      vjust=0.5
    )
  ) +
  labs(
    x=NULL,
    y=NULL,
    title="WT vs DHHKO differential expression by germ cell type"
  )
# order genes
plot_df$gene <- factor(plot_df$gene, levels = genes)

# Calculate mean log2FC across all cell types for gene ordering
gene_order <- plot_df %>%
  group_by(gene) %>%
  summarise(
    mean_log2FC = mean(avg_log2FC, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_log2FC)) %>%
  pull(gene)

# apply ordering
plot_df$gene <- factor(
  plot_df$gene,
  levels = gene_order
)

# plot
ggplot(plot_df,
       aes(
         x=gene,
         y=cell.type,
         size=pct.exp,
         color=avg_log2FC
       )) +
  geom_point() +
  scale_color_gradient2(
    low="#604e00",
    mid="white",
    high="#ba5656",
    midpoint=0,
    name="DHHKO vs WT\nlog2FC"
  ) +
  scale_size(
    range=c(1,8),
    name="% expressed"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle=90,
      hjust=1,
      vjust=0.5
    )
  ) +
  labs(
    x=NULL,
    y=NULL,
    title="WT vs DHHKO differential expression by germ cell type"
  )




###########################
## Revision 2 #############
###########################


# ============================================================
# WT vs WWv: NUMBER OF DEGs BY CELL TYPE
# ============================================================



# ------------------------------------------------------------
# 1. Keep only WT and WWv
# ------------------------------------------------------------

obj <- subset(
  all_mouse_merged_final,
  subset = treatment %in% c("WT", "WWv")
)
obj <- JoinLayers(obj)
Idents(obj) <- "treatment"
table(all_mouse_merged_final$final.cell.type)

# ------------------------------------------------------------
# 2. Run DE within each cell type
# ------------------------------------------------------------

cell.types <- unique(obj$final.cell.type)

deg.list <- list()

for (ct in cell.types) {
  
  message("Processing: ", ct)
  
  x <- subset(
    obj,
    subset = final.cell.type == ct
  )
  
  n.WT  <- sum(x$treatment == "WT")
  n.WWv <- sum(x$treatment == "WWv")
  
  # Skip cell types without enough cells in both groups
  if (n.WT < 10 | n.WWv < 10) {
    message("  Skipping: WT = ", n.WT, ", WWv = ", n.WWv)
    next
  }
  
  Idents(x) <- "treatment"
  
  deg <- FindMarkers(
    x,
    ident.1 = "WWv",
    ident.2 = "WT",
    test.use = "LR",
    latent.vars = "nCount_RNA",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  
  deg$gene <- rownames(deg)
  deg$cell.type <- ct
  
  deg.list[[ct]] <- deg
}


# ------------------------------------------------------------
# 3. Combine DEGs
# ------------------------------------------------------------

deg.all <- bind_rows(deg.list)

# Significant DEGs
deg.sig <- deg.all %>%
  filter(p_val_adj < 0.05)


# ------------------------------------------------------------
# 4. Count UP and DOWN genes
#
# Positive log2FC = higher in WWv
# Negative log2FC = lower in WWv
# ------------------------------------------------------------

deg.counts <- deg.sig %>%
  mutate(
    direction = ifelse(avg_log2FC > 0, "Up in WWv", "Down in WWv")
  ) %>%
  group_by(cell.type, direction) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    n.plot = ifelse(direction == "Down in WWv", -n, n)
  )

# ------------------------------------------------------------
# 4b. Keep SOMATIC cell types only and set plotting order
# ------------------------------------------------------------

somatic.order <- c(
  "Sertoli_cells",
  "Rete_epithelial",
  "Leydig_cells",
  "Peritubular_myoid_cells",
  "Perivascular_smooth_muscle_cells",
  "Mesenchymal_progenitors",
  "Endothelial_cells",
  "Macrophages",
  "Dendritic",
  "T_cells"
)

deg.counts <- deg.counts %>%
  filter(cell.type %in% somatic.order) %>%
  mutate(
    cell.type = factor(
      cell.type,
      levels = rev(somatic.order)
    )
  )

# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------

ggplot(
  deg.counts,
  aes(
    x = cell.type,
    y = n.plot,
    fill = direction
  )
) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  coord_flip() +
  scale_y_continuous(
    labels = abs
  ) +
  scale_fill_manual(
    values = c(
      "Down in WWv" = "#937700",  # higher in WT
      "Up in WWv"   = "#c9c9c9"   # higher in WWv
    ),
    labels = c(
      "Down in WWv" = "Higher in WT",
      "Up in WWv"   = "Higher in WWv"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of DEGs",
    fill = NULL,
    title = "WT vs WWv differential expression"
  ) +
  theme_classic(base_size = 12)






sertoli.genes <- c(
  "Sox9",
  "Sparc",
  "App",
  "Psap",
  "Dhh"
)

myoid.genes <- c(
  "Egr1",
  "Cebpd",
  "Jund",
  "Hspa8",
  "Igf1"
)

mesenchymal.genes <- c(
  "Klf4",
  "Fosb",
  "Hspa1a",
  "Igfbp4",
  "Col3a1"
)



table(all_mouse_merged_final$source)
Idents(all_mouse_merged_final) <- "source"
all_mouse_merged_our_data <- subset(all_mouse_merged_final, idents = "Guo", invert = T)
all_mouse_merged_our_data_wtwwv

Idents(all_mouse_merged_our_data) <- "final.cell.type"
obj_Sertoli_all <- subset(all_mouse_merged_our_data, idents = "Sertoli_cells")
Idents(obj_Sertoli_all) <- "treatment"
p1 <- VlnPlot(obj_Sertoli_all, features = sertoli.genes, pt.size = 0, ncol = 5)

Idents(all_mouse_merged_our_data) <- "final.cell.type"
table(Idents(all_mouse_merged_our_data))
obj_myoid_all <- subset(all_mouse_merged_our_data, idents = "Peritubular_myoid_cells")
Idents(obj_myoid_all) <- "treatment"
p2 <- VlnPlot(obj_myoid_all, features = myoid.genes, pt.size = 0, ncol = 5)

Idents(all_mouse_merged_our_data) <- "final.cell.type"
obj_mesenchymal_all <- subset(all_mouse_merged_our_data, idents = "Mesenchymal_progenitors")
Idents(obj_mesenchymal_all) <- "treatment"
p3 <- VlnPlot(obj_mesenchymal_all, features = mesenchymal.genes, pt.size = 0, ncol = 5)
plot_grid(p1, p2, p3, ncol = 1)

Idents(all_mouse_merged_our_data) <- "final.cell.type"
obj_Sertoli <- subset(all_mouse_merged_our_data, idents = "Sertoli_cells")
Idents(obj_Sertoli) <- "treatment"
p1 <- VlnPlot(obj_Sertoli, features = sertoli.genes, pt.size = 0, cols = c("#937700","#c9c9c9"), ncol = 5)

Idents(all_mouse_merged_our_data) <- "final.cell.type"
table(Idents(all_mouse_merged_our_data))
obj_myoid <- subset(all_mouse_merged_our_data, idents = "Peritubular_myoid_cells")
Idents(obj_myoid) <- "treatment"
p2 <- VlnPlot(obj_myoid, features = myoid.genes, pt.size = 0, cols = c("#937700","#c9c9c9"), ncol = 5)

Idents(all_mouse_merged_our_data) <- "final.cell.type"
obj_mesenchymal <- subset(all_mouse_merged_our_data, idents = "Mesenchymal_progenitors")
Idents(obj_mesenchymal) <- "treatment"
p3 <- VlnPlot(obj_mesenchymal, features = mesenchymal.genes, pt.size = 0, cols = c("#937700","#c9c9c9"), ncol = 5)

plot_grid(p1, p2, p3, ncol = 1)

# ============================================================
# SERTOLI CELL DEG OVERLAP:
# WT vs WWv
# Experimental_fertile vs Experimental_infertile
# ============================================================




# ------------------------------------------------------------
# 1. Subset Sertoli cells from the FULL object
# ------------------------------------------------------------

obj_Sertoli_all <- subset(
  all_mouse_merged_final,
  subset = final.cell.type == "Sertoli_cells"
)

DefaultAssay(obj_Sertoli_all) <- "RNA"

# Important for merged Seurat v5 object
obj_Sertoli_all <- JoinLayers(
  obj_Sertoli_all,
  assay = "RNA"
)

table(obj_Sertoli_all$treatment)


log2fc_threshold <- 0.25
min.pct <- 0.1

# ------------------------------------------------------------
# 2. WT vs WWv
#
# Positive FC = higher in WWv
# Negative FC = higher in WT
# ------------------------------------------------------------

sertoli_wv <- subset(
  obj_Sertoli_all,
  subset = treatment %in% c("WT", "WWv")
)

Idents(sertoli_wv) <- "treatment"

deg_WWv_vs_WT <- FindMarkers(
  sertoli_wv,
  ident.1 = "WWv",
  ident.2 = "WT",
  test.use = "LR",
  latent.vars = "nCount_RNA",
  logfc.threshold = log2fc_threshold,
  min.pct = min.pct
)

deg_WWv_vs_WT$gene <- rownames(deg_WWv_vs_WT)


# ------------------------------------------------------------
# 3. Experimental infertile vs fertile
#
# Positive FC = higher in Experimental_infertile
# Negative FC = higher in Experimental_fertile
# ------------------------------------------------------------

sertoli_exp <- subset(
  obj_Sertoli_all,
  subset = treatment %in% c(
    "Experimental_fertile",
    "Experimental_infertile"
  )
)

Idents(sertoli_exp) <- "treatment"

deg_infertile_vs_fertile <- FindMarkers(
  sertoli_exp,
  ident.1 = "Experimental_infertile",
  ident.2 = "Experimental_fertile",
  test.use = "LR",
  latent.vars = "nCount_RNA",
  logfc.threshold = log2fc_threshold,
  min.pct = min.pct
)

deg_infertile_vs_fertile$gene <- rownames(deg_infertile_vs_fertile)


# ------------------------------------------------------------
# 4. Keep significant DEGs
# ------------------------------------------------------------

sig_WWv_vs_WT <- deg_WWv_vs_WT %>%
  filter(p_val_adj < 0.05)

sig_infertile_vs_fertile <- deg_infertile_vs_fertile %>%
  filter(p_val_adj < 0.05)


# ------------------------------------------------------------
# 5. Make directional gene sets
# ------------------------------------------------------------

# Higher in abnormal / infertile state
WWv_up <- sig_WWv_vs_WT %>%
  filter(avg_log2FC > 0) %>%
  pull(gene) %>%
  unique()

infertile_up <- sig_infertile_vs_fertile %>%
  filter(avg_log2FC > 0) %>%
  pull(gene) %>%
  unique()


# Higher in control / fertile state
WT_up <- sig_WWv_vs_WT %>%
  filter(avg_log2FC < 0) %>%
  pull(gene) %>%
  unique()

fertile_up <- sig_infertile_vs_fertile %>%
  filter(avg_log2FC < 0) %>%
  pull(gene) %>%
  unique()


# ------------------------------------------------------------
# 6. Print overlap numbers
# ------------------------------------------------------------

cat("\n====================================\n")
cat("HIGHER IN WWv / INFERTILE\n")
cat("====================================\n")

cat("WWv:", length(WWv_up), "\n")
cat("Experimental infertile:", length(infertile_up), "\n")
cat(
  "Overlap:",
  length(intersect(WWv_up, infertile_up)),
  "\n"
)


cat("\n====================================\n")
cat("HIGHER IN WT / FERTILE\n")
cat("====================================\n")

cat("WT:", length(WT_up), "\n")
cat("Experimental fertile:", length(fertile_up), "\n")
cat(
  "Overlap:",
  length(intersect(WT_up, fertile_up)),
  "\n"
)


# ------------------------------------------------------------
# 7. Euler diagram:
# Higher in WWv / Experimental infertile
# ------------------------------------------------------------

fit_up <- euler(
  list(
    "WWv" = WWv_up,
    "Experimental infertile" = infertile_up
  )
)

p1 <- plot(
  fit_up,
  fills = c("#c9c9c9", "#c05127"),
  edges = "black",
  quantities = TRUE,
  labels = list(font = 2),
  main = "Sertoli cells: genes higher in infertile conditions"
)


# ------------------------------------------------------------
# 8. Euler diagram:
# Higher in WT / Experimental fertile
# ------------------------------------------------------------

fit_down <- euler(
  list(
    "WT" = WT_up,
    "Experimental fertile" = fertile_up
  )
)

p2 <- plot(
  fit_down,
  fills = c("#937700", "#48b1e0"),
  edges = "black",
  quantities = TRUE,
  labels = list(font = 2),
  main = "Sertoli cells: genes higher in fertile conditions"
)


plot_grid(p1, p2, ncol = 1)

# ============================================================
# SERTOLI CELL CONCORDANCE ANALYSIS
#
# Compare:
#
#   1. WWv vs WT
#   2. Experimental infertile vs Experimental fertile
#
# Goal:
# Determine whether the transcriptional response associated
# with infertility is concordant between the two comparisons.
# ============================================================


# ============================================================
# LIBRARIES
# ============================================================



# ============================================================
# 1. SETTINGS
# ============================================================

# Used to DEFINE DEGs after FindMarkers
deg_log2fc_cutoff <- 0.25
deg_padj_cutoff <- 0.05

# FindMarkers expression filter
min_pct_cutoff <- 0.1


# ============================================================
# 2. START FROM FULL OBJECT
# Keep only groups required for this analysis
# ============================================================

obj_overlap <- subset(
  all_mouse_merged_final,
  subset = treatment %in% c(
    "WT",
    "WWv",
    "Experimental_fertile",
    "Experimental_infertile"
  )
)

DefaultAssay(obj_overlap) <- "RNA"


# ============================================================
# 3. REMOVE GUO DATA
#
# We want this analysis to use only internally generated data.
# ============================================================

if ("source" %in% colnames(obj_overlap@meta.data)) {
  
  cat("\nSources before removal:\n")
  print(table(obj_overlap$source))
  
  obj_overlap <- subset(
    obj_overlap,
    subset = source != "Guo"
  )
  
  cat("\nSources after removal:\n")
  print(table(obj_overlap$source))
}


# ============================================================
# 4. KEEP RNA ASSAY ONLY
#
# Removes peaks / ADT so they cannot cause subsetting problems.
# ============================================================

DefaultAssay(obj_overlap) <- "RNA"

for (a in c("peaks", "ADT")) {
  
  if (a %in% names(obj_overlap@assays)) {
    obj_overlap[[a]] <- NULL
  }
  
}

cat("\nRemaining assays:\n")
print(names(obj_overlap@assays))


# ============================================================
# 5. JOIN RNA LAYERS
# ============================================================

obj_overlap <- JoinLayers(
  obj_overlap,
  assay = "RNA"
)

DefaultAssay(obj_overlap) <- "RNA"

cat("\nRNA layers:\n")
print(
  Layers(obj_overlap[["RNA"]])
)


# ============================================================
# 6. CHECK CELL NUMBERS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CELL NUMBERS\n")
cat("============================================================\n")

print(
  table(
    obj_overlap$final.cell.type,
    obj_overlap$treatment
  )
)


# ============================================================
# 7. SUBSET SERTOLI CELLS
# ============================================================

sertoli <- subset(
  obj_overlap,
  subset = final.cell.type == "Sertoli_cells"
)

DefaultAssay(sertoli) <- "RNA"

cat("\n")
cat("============================================================\n")
cat("SERTOLI CELL NUMBERS\n")
cat("============================================================\n")

print(
  table(sertoli$treatment)
)


# ============================================================
# 8. CREATE THE TWO COMPARISON OBJECTS
# ============================================================


# ------------------------------------------------------------
# WWv versus WT
# ------------------------------------------------------------

sertoli_wv <- subset(
  sertoli,
  subset = treatment %in% c(
    "WT",
    "WWv"
  )
)

DefaultAssay(sertoli_wv) <- "RNA"

Idents(sertoli_wv) <- "treatment"


# ------------------------------------------------------------
# Experimental infertile versus fertile
# ------------------------------------------------------------

sertoli_exp <- subset(
  sertoli,
  subset = treatment %in% c(
    "Experimental_fertile",
    "Experimental_infertile"
  )
)

DefaultAssay(sertoli_exp) <- "RNA"

Idents(sertoli_exp) <- "treatment"


# ============================================================
# 9. FINDMARKERS:
# WWv vs WT
#
# IMPORTANT:
# logfc.threshold = 0
#
# We want the full fold-change vector rather than selecting
# genes based on logFC before calculating concordance.
#
# Positive FC = higher in WWv / infertility
# Negative FC = higher in WT / fertility
# ============================================================

deg_wv_full <- FindMarkers(
  sertoli_wv,
  ident.1 = "WWv",
  ident.2 = "WT",
  test.use = "LR",
  latent.vars = "nCount_RNA",
  logfc.threshold = 0,
  min.pct = min_pct_cutoff
)

deg_wv_full$gene <- rownames(
  deg_wv_full
)


# ============================================================
# 10. FINDMARKERS:
# Experimental infertile vs fertile
#
# Positive FC = higher in Experimental infertile
# Negative FC = higher in Experimental fertile
# ============================================================

deg_exp_full <- FindMarkers(
  sertoli_exp,
  ident.1 = "Experimental_infertile",
  ident.2 = "Experimental_fertile",
  test.use = "LR",
  latent.vars = "nCount_RNA",
  logfc.threshold = 0,
  min.pct = min_pct_cutoff
)

deg_exp_full$gene <- rownames(
  deg_exp_full
)


# ============================================================
# 11. COMBINE THE TWO CONTRASTS BY GENE
# ============================================================

comparison <- inner_join(
  
  deg_wv_full %>%
    select(
      gene,
      FC_WWv_WT = avg_log2FC,
      padj_WWv_WT = p_val_adj
    ),
  
  deg_exp_full %>%
    select(
      gene,
      FC_ExpInf_ExpFert = avg_log2FC,
      padj_Exp = p_val_adj
    ),
  
  by = "gene"
)


cat(
  "\nGenes tested in both comparisons:",
  nrow(comparison),
  "\n"
)


# ============================================================
# 12. DEFINE DEGs
#
# DEG requires:
#
# adjusted P < 0.05
# AND
# |log2FC| >= 0.25
#
# This is applied AFTER FindMarkers.
# ============================================================

comparison <- comparison %>%
  
  mutate(
    
    significant_WWv =
      padj_WWv_WT < deg_padj_cutoff &
      abs(FC_WWv_WT) >= deg_log2fc_cutoff,
    
    significant_exp =
      padj_Exp < deg_padj_cutoff &
      abs(FC_ExpInf_ExpFert) >= deg_log2fc_cutoff,
    
    # significant in BOTH comparisons
    significant_both =
      significant_WWv &
      significant_exp,
    
    # Same direction in both comparisons
    concordant_sign =
      sign(FC_WWv_WT) ==
      sign(FC_ExpInf_ExpFert)
    
  )


# ============================================================
# 13. DEGs SIGNIFICANT IN BOTH COMPARISONS
# ============================================================

both_sig <- comparison %>%
  filter(
    significant_both
  )


n_shared_deg <- nrow(
  both_sig
)


n_concordant_deg <- sum(
  both_sig$concordant_sign,
  na.rm = TRUE
)


percent_concordant_sig <-
  
  100 *
  n_concordant_deg /
  n_shared_deg


cat("\n")
cat("============================================================\n")
cat("CONCORDANCE RESULTS\n")
cat("============================================================\n")

cat(
  "DEGs significant in both comparisons:",
  n_shared_deg,
  "\n"
)

cat(
  "Same direction in both comparisons:",
  n_concordant_deg,
  "\n"
)

cat(
  "Directional concordance:",
  round(percent_concordant_sig, 1),
  "%\n"
)


# ============================================================
# 14. OPTIONAL: SPEARMAN CORRELATION
#
# We calculate it and keep it available,
# but do NOT display it on the figure.
# ============================================================

spearman_result <- cor.test(
  comparison$FC_WWv_WT,
  comparison$FC_ExpInf_ExpFert,
  method = "spearman",
  exact = FALSE
)

spearman_rho <- unname(
  spearman_result$estimate
)

cat(
  "Spearman rho across all tested genes:",
  round(spearman_rho, 3),
  "\n"
)


# ============================================================
# 15. CLASSIFY GENES FOR PLOTTING
# ============================================================

comparison <- comparison %>%
  
  mutate(
    
    plot_group = case_when(
      
      # Significant in both and increased in infertility
      significant_both &
        FC_WWv_WT > 0 &
        FC_ExpInf_ExpFert > 0
      ~ "Higher in infertile",
      
      # Significant in both and increased in fertility
      significant_both &
        FC_WWv_WT < 0 &
        FC_ExpInf_ExpFert < 0
      ~ "Higher in fertile",
      
      # Significant in both but opposite directions
      significant_both &
        !concordant_sign
      ~ "Discordant",
      
      # Everything else
      TRUE
      ~ "Other genes"
    )
  )


# ============================================================
# 16. PRINT CATEGORY NUMBERS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PLOT GROUP COUNTS\n")
cat("============================================================\n")

print(
  table(comparison$plot_group)
)


# ============================================================
# 17. MANUSCRIPT-FOCUSED SERTOLI GENES
#
# These genes are particularly relevant to the Sertoli GRN
# and biological model in the manuscript.
# ============================================================



genes_to_label <- c(
  "Sox9",
  "Dhh",
  "App",
  "Jam2",
  "Inha",
  "Sparc",
  "Gstm7"
)


# Only label genes if:
#
# 1. they are significant in BOTH comparisons
# 2. they change in the same direction
#
# This prevents us from highlighting a manuscript gene that
# does not actually support the concordance analysis.
# ============================================================

# label_data <- comparison %>%
#   
#   filter(
#     gene %in% genes_to_label,
#     significant_both,
#     concordant_sign
#   )
final_label_genes <- c(
  # Higher in fertile
  "Dab1",
  "Cldn3",
  "Rhox5",
  "Lgals1",
  
  # Higher in infertile
  "App",
  "Jam2",
  "Nfe2l2",
  "Itga9"
)

label_data <- comparison.clean %>%
  filter(
    gene %in% final_label_genes,
    significant_both,
    concordant_sign
  )

cat("\n")
cat("============================================================\n")
cat("MANUSCRIPT GENES INCLUDED IN PLOT\n")
cat("============================================================\n")

print(
  label_data %>%
    select(
      gene,
      FC_WWv_WT,
      FC_ExpInf_ExpFert,
      padj_WWv_WT,
      padj_Exp
    )
)


# ============================================================
# 18. FIGURE ANNOTATION
# ============================================================

annotation_text <- paste0(
  
  "Concordant direction among DEGs\n",
  "significant in both = ",
  
  round(
    percent_concordant_sig,
    1
  ),
  
  "% (",
  
  n_concordant_deg,
  
  "/",
  
  n_shared_deg,
  
  ")"
)


# ============================================================
# 19. CONCORDANCE PLOT
#
# Plot order is deliberate:
#
# 1. grey genes FIRST
# 2. discordant DEGs
# 3. concordant DEGs
# 4. labelled manuscript genes
#
# Therefore informative DEGs remain visible on top.
# ============================================================

p_concordance <- ggplot(
  comparison,
  aes(
    x = FC_WWv_WT,
    y = FC_ExpInf_ExpFert
  )
) +
  
  
  # ----------------------------------------------------------
# ZERO REFERENCE LINES
# ----------------------------------------------------------

geom_hline(
  yintercept = 0,
  linewidth = 0.35,
  color = "grey60"
) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey60"
  ) +
  
  
  # ----------------------------------------------------------
# PERFECT 1:1 RELATIONSHIP
# ----------------------------------------------------------

geom_abline(
  slope = 1,
  intercept = 0,
  linetype = "dashed",
  linewidth = 0.5,
  color = "grey40"
) +
  
  
  # ----------------------------------------------------------
# BACKGROUND:
# ALL NON-SHARED DEGs / NON-DEGs
# ----------------------------------------------------------

geom_point(
  
  data = comparison %>%
    filter(
      plot_group == "Other genes"
    ),
  
  color = "#c9c9c9",
  
  size = 1,
  
  alpha = 0.30
) +
  
  
  # ----------------------------------------------------------
# DISCORDANT DEGs
# ----------------------------------------------------------

geom_point(
  
  data = comparison %>%
    filter(
      plot_group == "Discordant"
    ),
  
  aes(
    color = plot_group
  ),
  
  size = 1.5,
  
  alpha = 0.75
) +
  
  
  # ----------------------------------------------------------
# CONCORDANT DEGs
#
# Draw these AFTER grey genes so they sit on top.
# ----------------------------------------------------------

geom_point(
  
  data = comparison %>%
    filter(
      plot_group %in% c(
        "Higher in infertile",
        "Higher in fertile"
      )
    ),
  
  aes(
    color = plot_group
  ),
  
  size = 1.6,
  
  alpha = 0.85
) +
  
  
  # ----------------------------------------------------------
# COLORS
# ----------------------------------------------------------

scale_color_manual(
  
  values = c(
    
    "Higher in infertile" =
      "#c05127",
    
    "Higher in fertile" =
      "#937700",
    
    "Discordant" =
      "#555555"
  ),
  
  breaks = c(
    "Higher in infertile",
    "Higher in fertile",
    "Discordant"
  )
) +
  
  
  # ----------------------------------------------------------
# HIGHLIGHT MANUSCRIPT GENES
# ----------------------------------------------------------

geom_point(
  
  data = label_data,
  
  color = "black",
  
  size = 3
) +
  
  
  ggrepel::geom_text_repel(
    
    data = label_data,
    
    aes(
      label = gene
    ),
    
    color = "black",
    
    size = 4,
    
    fontface = "italic",
    
    box.padding = 0.5,
    
    point.padding = 0.25,
    
    max.overlaps = Inf,
    
    min.segment.length = 0
  ) +
  
  
  # ----------------------------------------------------------
# CONCORDANCE STATISTIC
# ----------------------------------------------------------

annotate(
  
  "text",
  
  x = -Inf,
  
  y = Inf,
  
  label = annotation_text,
  
  hjust = -0.05,
  
  vjust = 1.1,
  
  size = 4.5
) +
  
  
  # ----------------------------------------------------------
# AXIS LABELS
# ----------------------------------------------------------

labs(
  
  x = expression(
    log[2] *
      " fold change: W/W"^v *
      " vs WT"
  ),
  
  y = expression(
    log[2] *
      " fold change: Experimental infertile vs fertile"
  ),
  
  color = NULL
) +
  
  
  # ----------------------------------------------------------
# THEME
# ----------------------------------------------------------

theme_classic(
  base_size = 13
) +
  
  theme(
    
    legend.position = "bottom",
    
    legend.text = element_text(
      size = 11
    ),
    
    axis.title = element_text(
      size = 13
    ),
    
    axis.text = element_text(
      size = 11,
      color = "black"
    )
  )


# ============================================================
# 20. DISPLAY
# ============================================================

p_concordance






# ============================================================
# INVESTIGATE STRONGEST CONCORDANT GENES
# ============================================================

candidate_table <- comparison %>%
  filter(
    significant_both,
    concordant_sign
  ) %>%
  mutate(
    # Conservative score:
    # gene must have a strong effect in BOTH comparisons
    min_abs_FC = pmin(
      abs(FC_WWv_WT),
      abs(FC_ExpInf_ExpFert)
    ),
    
    mean_abs_FC = (
      abs(FC_WWv_WT) +
        abs(FC_ExpInf_ExpFert)
    ) / 2
  )


# ============================================================
# TOP RIGHT
# Higher in BOTH infertile conditions
# ============================================================

top_infertile <- candidate_table %>%
  filter(
    FC_WWv_WT > 0,
    FC_ExpInf_ExpFert > 0
  ) %>%
  arrange(desc(min_abs_FC)) %>%
  select(
    gene,
    FC_WWv_WT,
    FC_ExpInf_ExpFert,
    padj_WWv_WT,
    padj_Exp,
    min_abs_FC
  ) %>%
  head(30)

cat("\n====================================\n")
cat("TOP CONCORDANT: HIGHER IN INFERTILE\n")
cat("====================================\n")

print(top_infertile)


# ============================================================
# BOTTOM LEFT
# Higher in BOTH fertile conditions
# ============================================================

top_fertile <- candidate_table %>%
  filter(
    FC_WWv_WT < 0,
    FC_ExpInf_ExpFert < 0
  ) %>%
  arrange(desc(min_abs_FC)) %>%
  select(
    gene,
    FC_WWv_WT,
    FC_ExpInf_ExpFert,
    padj_WWv_WT,
    padj_Exp,
    min_abs_FC
  ) %>%
  head(50)

cat("\n====================================\n")
cat("TOP CONCORDANT: HIGHER IN FERTILE\n")
cat("====================================\n")

print(top_fertile)




# ============================================================
# REDO WITH THE GERM CELL GENES REMOVED
# ============================================================


# ============================================================
# 1. RECREATE THE ORIGINAL SPERMATID GENE BLACKLIST
# ============================================================

# Original analysis used WT + DHHKO
DefaultAssay(all_mouse_merged_final) <- "RNA"
Idents(all_mouse_merged_final) <- "treatment"

blacklist.obj <- subset(
  all_mouse_merged_final,
  idents = c("WT", "DHHKO")
)

# RNA only avoids problems from other assays
for (a in c("peaks", "ADT")) {
  if (a %in% names(blacklist.obj@assays)) blacklist.obj[[a]] <- NULL
}

DefaultAssay(blacklist.obj) <- "RNA"
blacklist.obj <- JoinLayers(blacklist.obj, assay = "RNA")
Idents(blacklist.obj) <- "final.cell.type"

# Same strategy as original script
all.markers.blacklist <- FindAllMarkers(
  blacklist.obj,
  assay = "RNA",
  only.pos = TRUE
)

all.markers.blacklist$pct.diff <-
  all.markers.blacklist$pct.1 - all.markers.blacklist$pct.2

spermatid.clusters <- c(
  "Round_spermatids",
  "Elongating_spermatids"
)

spermatid.genes <- all.markers.blacklist %>%
  filter(
    cluster %in% spermatid.clusters,
    avg_log2FC > 2,
    pct.diff > 0.20,
    p_val_adj < 0.05
  ) %>%
  pull(gene) %>%
  unique()

spermatid.blacklist <- spermatid.genes

cat("Number of spermatid genes in blacklist:", length(spermatid.genes), "\n")

# Sanity check: these should mostly be TRUE
check.genes <- c("Prm1","Prm2","Tnp1","Tnp2","Ldhc","Fabp9","Smcp")
print(data.frame(
  gene = check.genes,
  blacklisted = check.genes %in% spermatid.genes
))


# ============================================================
# 2. APPLY BLACKLIST TO CURRENT CONCORDANCE ANALYSIS
#
# 'comparison' should already contain:
# gene
# FC_WWv_WT
# FC_ExpInf_ExpFert
# padj_WWv_WT
# padj_Exp
# ============================================================

cat("Genes before blacklist:", nrow(comparison), "\n")

comparison.clean <- comparison %>%
  filter(!gene %in% spermatid.blacklist)

cat("Genes after blacklist:", nrow(comparison.clean), "\n")
cat("Genes removed:", nrow(comparison) - nrow(comparison.clean), "\n")


# ============================================================
# 3. REDEFINE SIGNIFICANCE AFTER REMOVING BLACKLIST
# ============================================================

deg_log2fc_cutoff <- 0.25
deg_padj_cutoff <- 0.05

comparison.clean <- comparison.clean %>%
  mutate(
    significant_WWv =
      padj_WWv_WT < deg_padj_cutoff &
      abs(FC_WWv_WT) >= deg_log2fc_cutoff,
    
    significant_exp =
      padj_Exp < deg_padj_cutoff &
      abs(FC_ExpInf_ExpFert) >= deg_log2fc_cutoff,
    
    significant_both =
      significant_WWv & significant_exp,
    
    concordant_sign =
      sign(FC_WWv_WT) == sign(FC_ExpInf_ExpFert)
  )


# ============================================================
# 4. RECALCULATE CONCORDANCE
# ============================================================

both_sig.clean <- comparison.clean %>%
  filter(significant_both)

n_shared_deg <- nrow(both_sig.clean)
n_concordant_deg <- sum(both_sig.clean$concordant_sign, na.rm = TRUE)
percent_concordant_sig <- 100 * n_concordant_deg / n_shared_deg

cat("\nCLEANED RESULTS\n")
cat("DEGs significant in both:", n_shared_deg, "\n")
cat("Concordant DEGs:", n_concordant_deg, "\n")
cat("Percent concordant:", round(percent_concordant_sig, 1), "%\n")


# ============================================================
# 5. CLASSIFY FOR PLOT
# ============================================================

comparison.clean <- comparison.clean %>%
  mutate(
    plot_group = case_when(
      significant_both &
        FC_WWv_WT > 0 &
        FC_ExpInf_ExpFert > 0 ~ "Higher in infertile",
      
      significant_both &
        FC_WWv_WT < 0 &
        FC_ExpInf_ExpFert < 0 ~ "Higher in fertile",
      
      significant_both &
        !concordant_sign ~ "Discordant",
      
      TRUE ~ "Other genes"
    )
  )

print(table(comparison.clean$plot_group))


# ============================================================
# 6. LOOK AGAIN AT STRONGEST CLEANED CANDIDATES
# ============================================================

candidate.clean <- comparison.clean %>%
  filter(significant_both, concordant_sign) %>%
  mutate(
    min_abs_FC = pmin(
      abs(FC_WWv_WT),
      abs(FC_ExpInf_ExpFert)
    )
  )

top.infertile.clean <- candidate.clean %>%
  filter(FC_WWv_WT > 0, FC_ExpInf_ExpFert > 0) %>%
  arrange(desc(min_abs_FC)) %>%
  select(gene, FC_WWv_WT, FC_ExpInf_ExpFert, min_abs_FC) %>%
  head(30)

top.fertile.clean <- candidate.clean %>%
  filter(FC_WWv_WT < 0, FC_ExpInf_ExpFert < 0) %>%
  arrange(desc(min_abs_FC)) %>%
  select(gene, FC_WWv_WT, FC_ExpInf_ExpFert, min_abs_FC) %>%
  head(30)

print(top.infertile.clean)
print(top.fertile.clean)


# ============================================================
# 7. GENES TO LABEL
# ============================================================

# Higher in fertile conditions / lower in infertility
fertile_label_genes <- c(
  "Aif1",
  "Dab1",
  "Cldn3",
  "Chst2",
  "Cldn10"
)

# Higher in infertile conditions
infertile_label_genes <- c(
  "App",
  "Jam2",
  "Nfe2l2",
  "Itga9",
  "Pdgfa",
  "Gdnf",
  "Fgf11",
  "Wnt5a",
  "Jag1",
  "Csf1",
  "Cxcl12"
)


# ------------------------------------------------------------
# Only label genes that:
# 1. are significant in BOTH comparisons
# 2. move concordantly
# 3. move in the expected direction
# ------------------------------------------------------------

label_fertile <- comparison.clean %>%
  filter(
    gene %in% fertile_label_genes,
    significant_both,
    FC_WWv_WT < 0,
    FC_ExpInf_ExpFert < 0
  )

label_infertile <- comparison.clean %>%
  filter(
    gene %in% infertile_label_genes,
    significant_both,
    FC_WWv_WT > 0,
    FC_ExpInf_ExpFert > 0
  )

label_data <- bind_rows(
  label_fertile,
  label_infertile
)


# ============================================================
# 8. CHECK WHICH CANDIDATES ACTUALLY QUALIFY
# ============================================================

label_check <- comparison.clean %>%
  filter(
    gene %in% c(
      fertile_label_genes,
      infertile_label_genes
    )
  ) %>%
  mutate(
    expected_group = case_when(
      gene %in% fertile_label_genes ~ "Higher in fertile",
      gene %in% infertile_label_genes ~ "Higher in infertile"
    )
  ) %>%
  select(
    gene,
    expected_group,
    FC_WWv_WT,
    FC_ExpInf_ExpFert,
    padj_WWv_WT,
    padj_Exp,
    significant_WWv,
    significant_exp,
    significant_both,
    concordant_sign
  ) %>%
  arrange(expected_group, gene)

print(label_check)


# ============================================================
# 8. REDRAW CLEANED CONCORDANCE PLOT
# ============================================================

annotation_text <- paste0(
  "Concordant direction among DEGs\n",
  "significant in both = ",
  round(percent_concordant_sig, 1),
  "% (", n_concordant_deg, "/", n_shared_deg, ")"
)

p_concordance_clean <- ggplot(
  comparison.clean,
  aes(x = FC_WWv_WT, y = FC_ExpInf_ExpFert)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.5, color = "grey40") +
  
  # Grey genes first
  geom_point(
    data = comparison.clean %>% filter(plot_group == "Other genes"),
    color = "#c9c9c9", size = 1, alpha = 0.30
  ) +
  
  # Discordant DEGs
  geom_point(
    data = comparison.clean %>% filter(plot_group == "Discordant"),
    aes(color = plot_group), size = 1.5, alpha = 0.75
  ) +
  
  # Concordant DEGs on top
  geom_point(
    data = comparison.clean %>%
      filter(plot_group %in% c("Higher in infertile", "Higher in fertile")),
    aes(color = plot_group), size = 1.6, alpha = 0.85
  ) +
  
  scale_color_manual(
    values = c(
      "Higher in infertile" = "#c05127",
      "Higher in fertile" = "#48b1e0",
      "Discordant" = "#555555"
    ),
    breaks = c(
      "Higher in infertile",
      "Higher in fertile",
      "Discordant"
    )
  ) +
  
  
  geom_point(
    data = label_data,
    color = "black",
    size = 3
  ) +
  
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    color = "black",
    size = 4,
    fontface = "italic",
    box.padding = 0.5,
    point.padding = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = annotation_text,
    hjust = -0.05, vjust = 1.1,
    size = 4.5
  ) +
  
  labs(
    x = expression(log[2] * " fold change: W/W"^v * " vs WT"),
    y = expression(log[2] * " fold change: Experimental infertile vs fertile"),
    color = NULL
  ) +
  
  coord_cartesian(
    xlim = c(-8, 7),
    ylim = c(-7, 5),
    clip = "off"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom"
  )

p_concordance_clean

Idents(sertoli) <- "source"
sertoli_thispaper <- subset(sertoli, idents = "ThisPaper")
VlnPlot(sertoli_thispaper, features = c("Dhh", "Sox9"), group.by = "treatment")





# how about a regular volcano plot?


# ============================================================
# VOLCANO PLOT: W/Wv vs WT SERTOLI CELLS
# Uses cleaned ThisPaper-only DEG results
# ============================================================



# Recreate cleaned W/Wv vs WT DEG table
deg_wv_clean <- deg_wv_full %>%
  mutate(gene = rownames(deg_wv_full)) %>%
  filter(!gene %in% spermatid.blacklist)

cat("Genes before blacklist:", nrow(deg_wv_full), "\n")
cat("Genes after blacklist:", nrow(deg_wv_clean), "\n")

# Same thresholds as your previous manuscript volcano plots
fc_cutoff <- 1
padj_cutoff <- 0.05

# Genes to label
genes_to_label <- c(
  "Dab1",
  "Cldn3",
  "Rhox5",
  "Lgals1",
  "App",
  "Jam2",
  "Nfe2l2",
  "Itga9",
  "Sox9"
)

# ============================================================
# PREPARE DATA
#
# deg_wv_clean was generated as:
# ident.1 = WWv
# ident.2 = WT
#
# Therefore:
# positive log2FC = higher in W/Wv
# negative log2FC = higher in WT
# ============================================================

volcano_wv <- deg_wv_clean %>%
  mutate(
    p_val_adj_clean = ifelse(
      is.na(p_val_adj) | p_val_adj == 0,
      .Machine$double.xmin,
      p_val_adj
    ),
    
    negLog10Padj = -log10(p_val_adj_clean),
    
    # Cap only for plotting
    negLog10Padj_plot = pmin(negLog10Padj, 300),
    avg_log2FC_plot = pmax(pmin(avg_log2FC, 8), -8),
    
    group = case_when(
      p_val_adj < padj_cutoff &
        avg_log2FC >= fc_cutoff ~ "W/Wv",
      
      p_val_adj < padj_cutoff &
        avg_log2FC <= -fc_cutoff ~ "WT",
      
      TRUE ~ "NS"
    )
  )

# Only label genes that pass the volcano thresholds
label_df <- volcano_wv %>%
  filter(
    gene %in% genes_to_label,
    group != "NS"
  )

# See which labels qualify
print(
  volcano_wv %>%
    filter(gene %in% genes_to_label) %>%
    select(
      gene,
      avg_log2FC,
      p_val_adj,
      group
    )
)

# ============================================================
# PLOT
# ============================================================

p_volcano_wwv_wt <- ggplot(
  volcano_wv,
  aes(
    x = avg_log2FC,
    y = negLog10Padj_plot
  )
) +
  geom_point(
    aes(fill = group),
    shape = 21,
    colour = "grey25",
    stroke = 0.1,
    size = 1.6,
    alpha = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      "NS" = "grey80",
      "WT" = "#937700",
      "W/Wv" = "#545249"
    )
  ) +
  
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    colour = "grey75",
    linetype = "dotted"
  ) +
  
  geom_hline(
    yintercept = -log10(padj_cutoff),
    colour = "grey75",
    linetype = "dotted"
  ) +
  
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 4,
    fontface = "italic",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    seed = 123
  ) +
  
  coord_cartesian(
    xlim = c(-7.25, 7),
    ylim = c(0, 300)
  ) +
  
  labs(
    x = expression(log[2] * " fold change: W/W"^v * " vs WT"),
    y = expression(-log[10]("adjusted " * italic(p))),
    title = "W/Wv vs WT Sertoli cells",
    fill = NULL
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black")
  )

p_volcano_wwv_wt







# ============================================================
# EXPORT WT vs W/Wv SERTOLI DEGs FOR INGENUITY
# ============================================================

ipa_wwv_wt <- deg_wv_clean %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) >= 1
  ) %>%
  transmute(
    GeneSymbol = gene,
    Log2FC = avg_log2FC,
    AdjP = p_val_adj
  ) %>%
  arrange(AdjP)

write.table(
  ipa_wwv_wt,
  file = "WT_vs_WWv_Sertoli_DEGs_IPA.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = ""
)

cat("Exported", nrow(ipa_wwv_wt), "DEGs\n")



# ============================================================
# WT vs W/Wv: NUMBER OF DEGs BY SOMATIC CELL TYPE
#
# Reviewer analysis:
# - ThisPaper only
# - WT vs W/Wv
# - LR test adjusting for nCount_RNA
# - min.pct = 0.1
# - DEG = adjusted P < 0.05 AND |log2FC| >= 1
# - spermatid.blacklist removed after DE testing
# ============================================================



# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

fc_cutoff <- 1
padj_cutoff <- 0.05
min_pct_cutoff <- 0.1

stopifnot(exists("spermatid.blacklist"))

somatic.order <- c(
  "Sertoli_cells",
  "Rete_epithelial",
  "Leydig_cells",
  "Peritubular_myoid_cells",
  "Perivascular_smooth_muscle_cells",
  "Mesenchymal_progenitors",
  "Endothelial_cells",
  "Macrophages",
  "Dendritic",
  "T_cells"
)

# ------------------------------------------------------------
# 2. Keep ThisPaper WT + W/Wv only
# ------------------------------------------------------------

obj_bar <- subset(
  all_mouse_merged_final,
  subset =
    source == "ThisPaper" &
    treatment %in% c("WT", "WWv")
)

DefaultAssay(obj_bar) <- "RNA"

# Remove non-RNA assays to avoid subsetting/layer problems
for (a in c("peaks", "ADT")) {
  if (a %in% names(obj_bar@assays)) {
    obj_bar[[a]] <- NULL
  }
}

obj_bar <- JoinLayers(
  obj_bar,
  assay = "RNA"
)

cat("\nCell numbers:\n")
print(
  table(
    obj_bar$final.cell.type,
    obj_bar$treatment
  )
)

# ------------------------------------------------------------
# 3. Run WT vs W/Wv DE within each somatic cell type
#
# Positive FC = higher in W/Wv
# Negative FC = higher in WT
# ------------------------------------------------------------

deg_bar_list <- list()

for (ct in somatic.order) {
  
  message("Processing: ", ct)
  
  x <- subset(
    obj_bar,
    subset = final.cell.type == ct
  )
  
  n_WT <- sum(x$treatment == "WT")
  n_WWv <- sum(x$treatment == "WWv")
  
  # Same general safeguard as before
  if (n_WT < 10 | n_WWv < 10) {
    message(
      "  Skipping: WT = ", n_WT,
      ", W/Wv = ", n_WWv
    )
    next
  }
  
  Idents(x) <- "treatment"
  
  deg <- FindMarkers(
    x,
    ident.1 = "WWv",
    ident.2 = "WT",
    test.use = "LR",
    latent.vars = "nCount_RNA",
    logfc.threshold = 0,
    min.pct = min_pct_cutoff
  )
  
  deg$gene <- rownames(deg)
  deg$cell.type <- ct
  
  # Apply the same spermatid blacklist used elsewhere
  deg <- deg %>%
    filter(!gene %in% spermatid.blacklist)
  
  deg_bar_list[[ct]] <- deg
}

# ------------------------------------------------------------
# 4. Combine results
# ------------------------------------------------------------

deg_bar_all <- bind_rows(deg_bar_list)

# Paper-style DEG cutoff:
# adjusted P < 0.05 AND |log2FC| >= 1
deg_bar_sig <- deg_bar_all %>%
  filter(
    p_val_adj < padj_cutoff,
    abs(avg_log2FC) >= fc_cutoff
  )

cat("\nTotal significant DEGs after blacklist:",
    nrow(deg_bar_sig), "\n")

# ------------------------------------------------------------
# 5. Count genes by direction
# ------------------------------------------------------------

deg_bar_counts <- deg_bar_sig %>%
  mutate(
    direction = ifelse(
      avg_log2FC > 0,
      "Higher in W/Wv",
      "Higher in WT"
    )
  ) %>%
  group_by(cell.type, direction) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

# Add zero-count combinations so every bar is represented
deg_bar_counts <- expand.grid(
  cell.type = somatic.order,
  direction = c("Higher in WT", "Higher in W/Wv"),
  stringsAsFactors = FALSE
) %>%
  left_join(
    deg_bar_counts,
    by = c("cell.type", "direction")
  ) %>%
  mutate(
    n = ifelse(is.na(n), 0, n),
    
    # WT goes left, W/Wv goes right
    n.plot = ifelse(
      direction == "Higher in WT",
      -n,
      n
    ),
    
    cell.type = factor(
      cell.type,
      levels = rev(somatic.order)
    )
  )

# ------------------------------------------------------------
# 6. Print counts
# ------------------------------------------------------------

print(
  deg_bar_counts %>%
    arrange(cell.type, direction)
)

# ------------------------------------------------------------
# 7. Plot
#
# WT = #937700
# W/Wv = #545249
# ------------------------------------------------------------

p_deg_bar <- ggplot(
  deg_bar_counts,
  aes(
    x = cell.type,
    y = n.plot,
    fill = direction
  )
) +
  
  geom_col(
    width = 0.75
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    labels = abs
  ) +
  
  scale_fill_manual(
    values = c(
      "Higher in WT" = "#937700",
      "Higher in W/Wv" = "#545249"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Number of DEGs",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black")
  )

p_deg_bar


# ============================================================
# LEYDIG ABUNDANCE AS % OF ALL SOMATIC CELLS
#
# Fresh-session script
#
# WT / WWv / DHHKO:
#   Unselected samples
#
# Experimental infertile:
#   Unselected control samples
#
# Experimental fertile:
#   GFP-negative fraction
#   + WWv_2_GFP1 / WWv_2_GFP2
#
# Technical replicates from the same biological mouse are
# pooled BEFORE calculating Leydig %.
#
# Each dot = one biological mouse
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================



# ============================================================
# 2. LOAD MERGED OBJECT
# ============================================================

object_file <- paste0(
  "/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles/",
  "dhh_wt_experimental_wv_all_merged_v6.Robj"
)

load(object_file)

stopifnot(exists("all_mouse_merged_final"))


# ============================================================
# 3. TREATMENT ORDER / COLOURS
# ============================================================

treatment_order <- c(
  "WT",
  "WWv",
  "DHHKO",
  "Experimental_fertile",
  "Experimental_infertile"
)

treatment_cols <- c(
  "WT"                     = "#937700",
  "WWv"                    = "#545249",
  "DHHKO"                  = "#f26b68",
  "Experimental_fertile"   = "#48b1e0",
  "Experimental_infertile" = "#c05127"
)


# ============================================================
# 4. DEFINE SOMATIC CELL TYPES
# ============================================================

somatic_celltypes <- c(
  "Sertoli_cells",
  "Rete_epithelial",
  "Leydig_cells",
  "Peritubular_myoid_cells",
  "Perivascular_smooth_muscle_cells",
  "Mesenchymal_progenitors",
  "Endothelial_cells",
  "Macrophages",
  "Dendritic",
  "T_cells"
)


# ============================================================
# 5. START FROM METADATA AND REMOVE GUO
# ============================================================

meta <- all_mouse_merged_final@meta.data %>%
  as.data.frame() %>%
  dplyr::filter(
    is.na(source) | source != "Guo"
  )


# ============================================================
# 6. SELECT APPROPRIATE FRACTION FOR EACH TREATMENT
#
# WT / WWv / DHHKO:
#   Unselected
#
# Experimental infertile:
#   Unselected
#
# Experimental fertile:
#   GFP-negative
#
# WWv_2_GFP1 / GFP2 are older fertile samples which are
# annotated as Unselected, so retain them explicitly.
# ============================================================

meta_use <- meta %>%
  dplyr::filter(
    dplyr::case_when(
      
      treatment %in% c(
        "WT",
        "WWv",
        "DHHKO"
      ) &
        selection == "Unselected" ~ TRUE,
      
      treatment == "Experimental_infertile" &
        selection == "Unselected" ~ TRUE,
      
      treatment == "Experimental_fertile" &
        selection == "GFP_neg" ~ TRUE,
      
      treatment == "Experimental_fertile" &
        orig.ident %in% c(
          "WWv_2_GFP1",
          "WWv_2_GFP2"
        ) ~ TRUE,
      
      TRUE ~ FALSE
    )
  )


# ============================================================
# 7. ASSIGN BIOLOGICAL MOUSE IDs
#
# Technical replicates from the same animal are given the
# same mouse_id.
#
# Anything not explicitly mapped below retains orig.ident.
# Therefore DHHKO15 and DHHKO18 remain separate mice.
# ============================================================

meta_use <- meta_use %>%
  dplyr::mutate(
    
    mouse_id = dplyr::case_when(
      
      orig.ident %in% c(
        "WWv_2_control",
        "WWv_2_GFP1",
        "WWv_2_GFP2"
      ) ~ "WWv2",
      
      orig.ident %in% c(
        "WWv3control",
        "WWv3GFPnegative"
      ) ~ "WWv3",
      
      orig.ident %in% c(
        "WWv5control",
        "WWv5GFPnegative"
      ) ~ "WWv5",
      
      orig.ident %in% c(
        "WWv6control",
        "WWv6GFPnegative"
      ) ~ "WWv6",
      
      orig.ident %in% c(
        "WWv7control",
        "WWv7GFPnegative"
      ) ~ "WWv7",
      
      orig.ident %in% c(
        "multiome1CTL",
        "multiome1GFPnegA",
        "multiome1GFPnegB"
      ) ~ "multiome1",
      
      orig.ident %in% c(
        "multiome2CTL",
        "multiome2GFPnegA",
        "multiome2GFPnegB"
      ) ~ "multiome2",
      
      TRUE ~ as.character(orig.ident)
    )
  )


# ============================================================
# 8. CHECK EXACTLY WHICH SAMPLES ARE INCLUDED
# ============================================================

sample_table <- meta_use %>%
  dplyr::distinct(
    treatment,
    mouse_id,
    orig.ident,
    selection
  ) %>%
  dplyr::arrange(
    factor(treatment, levels = treatment_order),
    mouse_id,
    orig.ident
  )

print(
  sample_table
  # n = Inf
)


# ============================================================
# 9. CALCULATE LEYDIG FRACTION PER BIOLOGICAL MOUSE
#
# Leydig fraction =
#
#       Leydig cells
# -------------------------
#    all somatic cells
#
#
# IMPORTANT:
# Technical replicates are combined at the CELL COUNT level
# before the percentage is calculated.
# ============================================================

leydig_mouse <- meta_use %>%
  dplyr::filter(
    treatment %in% treatment_order,
    final.cell.type %in% somatic_celltypes
  ) %>%
  
  dplyr::group_by(
    treatment,
    mouse_id
  ) %>%
  
  dplyr::summarise(
    
    total_somatic = dplyr::n(),
    
    Leydig_cells = sum(
      final.cell.type == "Leydig_cells"
    ),
    
    Leydig_fraction =
      Leydig_cells / total_somatic,
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    treatment = factor(
      treatment,
      levels = treatment_order
    )
  )


# ============================================================
# 10. INSPECT VALUES
# ============================================================

cat("\nLeydig abundance by biological mouse:\n")

print(
  leydig_mouse %>%
    dplyr::arrange(
      treatment,
      mouse_id
    ),
  n = Inf
)


cat("\nBiological replicates per treatment:\n")

print(
  leydig_mouse %>%
    dplyr::count(
      treatment,
      name = "n_mice"
    )
)


# ============================================================
# 11. IDENTIFY MATCHED EXPERIMENTAL PAIRS
# ============================================================

paired_mice <- leydig_mouse %>%
  dplyr::filter(
    treatment %in% c(
      "Experimental_fertile",
      "Experimental_infertile"
    )
  ) %>%
  
  dplyr::group_by(mouse_id) %>%
  
  dplyr::filter(
    dplyr::n_distinct(treatment) == 2
  ) %>%
  
  dplyr::ungroup()


cat("\nMatched experimental mice:\n")

print(
  paired_mice %>%
    dplyr::arrange(
      mouse_id,
      treatment
    ),
  n = Inf
)


# ============================================================
# 12. PLOT
#
# Experimental fertile / infertile samples from the same
# mouse are connected by a grey line.
# ============================================================

p_leydig <- ggplot(
  leydig_mouse,
  aes(
    x = treatment,
    y = Leydig_fraction
  )
) +
  
  # Boxplots
  geom_boxplot(
    aes(fill = treatment),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  
  # Paired experimental lines
  geom_line(
    data = paired_mice,
    aes(
      x = treatment,
      y = Leydig_fraction,
      group = mouse_id
    ),
    inherit.aes = FALSE,
    colour = "grey55",
    linewidth = 0.7,
    alpha = 0.8
  ) +
  
  # Individual biological samples
  geom_point(
    aes(fill = treatment),
    shape = 21,
    colour = "black",
    size = 3,
    stroke = 0.5,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  
  scale_fill_manual(
    values = treatment_cols,
    drop = FALSE
  ) +
  
  # Explicitly force desired order
  scale_x_discrete(
    limits = treatment_order,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    x = NULL,
    y = "Leydig cells (% of somatic cells)"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      colour = "black"
    )
  )


p_leydig


# ============================================================
# 13. PAIRED EXPERIMENTAL STATISTICS
#
# Experimental fertile vs experimental infertile
#
# Paired because the two testes are derived from the same mouse.
# ============================================================

paired_wide <- paired_mice %>%
  dplyr::select(
    mouse_id,
    treatment,
    Leydig_fraction
  ) %>%
  
  tidyr::pivot_wider(
    names_from = treatment,
    values_from = Leydig_fraction
  )


cat("\nPaired experimental values:\n")

print(
  paired_wide,
  n = Inf
)


cat("\nPaired Wilcoxon: Experimental fertile vs infertile\n")

print(
  wilcox.test(
    paired_wide$Experimental_fertile,
    paired_wide$Experimental_infertile,
    paired = TRUE,
    exact = FALSE
  )
)


# ============================================================
# 14. WT VS WWv
#
# Independent biological samples
# ============================================================

wt_values <- leydig_mouse %>%
  dplyr::filter(
    treatment == "WT"
  ) %>%
  dplyr::pull(Leydig_fraction)

wwv_values <- leydig_mouse %>%
  dplyr::filter(
    treatment == "WWv"
  ) %>%
  dplyr::pull(Leydig_fraction)


cat("\nWilcoxon: WT vs WWv\n")

print(
  wilcox.test(
    wt_values,
    wwv_values,
    exact = FALSE
  )
)


# ============================================================
# 15. WT VS DHHKO
#
# Independent biological samples
# ============================================================

dhhko_values <- leydig_mouse %>%
  dplyr::filter(
    treatment == "DHHKO"
  ) %>%
  dplyr::pull(Leydig_fraction)


cat("\nWilcoxon: WT vs DHHKO\n")

print(
  wilcox.test(
    wt_values,
    dhhko_values,
    exact = FALSE
  )
)


# ============================================================
# 16. OPTIONAL OVERALL TEST ACROSS ALL FIVE GROUPS
# ============================================================

cat("\nOverall Kruskal-Wallis test:\n")

print(
  kruskal.test(
    Leydig_fraction ~ treatment,
    data = leydig_mouse
  )
)



# ============================================================
# REVIEWER 2: LEYDIG ANDROGEN PRODUCTION / SERTOLI RESPONSE
# Paired experimental infertile vs GFP-transplanted fertile
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(msigdbr)
})

fertility_colors <- c(
  "Infertile" = "#c05127",
  "Fertile"   = "#48b1e0"
)

# ------------------------------------------------------------
# 1. Set up experimental object
# ------------------------------------------------------------

androgen_obj <- somatic.integrated.somatic6
DefaultAssay(androgen_obj) <- "RNA"

# Rename conditions for clarity
androgen_obj$fertility_status <- dplyr::recode(
  androgen_obj$treatment,
  "control" = "Infertile",
  "GFP transplanted" = "Fertile"
)

androgen_obj$fertility_status <- factor(
  androgen_obj$fertility_status,
  levels = c("Fertile", "Infertile")
)

table(
  androgen_obj$replicate,
  androgen_obj$cell.type3,
  androgen_obj$fertility_status
)


# ============================================================
# 2. LEYDIG CELLS: ANDROGEN/STEROIDOGENIC PROGRAM
# ============================================================

leydig <- subset(
  androgen_obj,
  subset = cell.type3 == "Leydig_cells"
)

androgen_synthesis_genes <- c(
  "Scarb1",
  "Star",
  "Cyp11a1",
  "Hsd3b1",
  "Hsd3b6",
  "Cyp17a1",
  "Hsd17b3",
  "Lhcgr"
)

# Keep only genes actually present in the dataset
androgen_synthesis_genes <- intersect(
  androgen_synthesis_genes,
  rownames(leydig)
)

cat(
  "\nLeydig steroidogenesis genes used:\n",
  paste(androgen_synthesis_genes, collapse = ", "),
  "\n"
)

# Cell-level module score
leydig <- AddModuleScore(
  leydig,
  features = list(androgen_synthesis_genes),
  name = "AndrogenSynthesis"
)

# Collapse to one value per mouse/testis
leydig_mouse <- leydig@meta.data %>%
  group_by(replicate, fertility_status) %>%
  summarise(
    androgen_synthesis_score = mean(
      AndrogenSynthesis1,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

print(leydig_mouse)


# ------------------------------------------------------------
# Paired plot: Leydig steroidogenic score
# ------------------------------------------------------------

p_leydig_score <- ggplot(
  leydig_mouse,
  aes(
    x = fertility_status,
    y = androgen_synthesis_score,
    group = replicate
  )
) +
  geom_line(alpha = 0.45) +
  geom_point(
    aes(fill = fertility_status),
    shape = 21,
    size = 3,
    stroke = 0.7
  ) +
  scale_fill_manual(values = fertility_colors)+
  theme_classic() +
  labs(
    x = NULL,
    y = "Leydig steroidogenesis score"
  ) +
  theme(
    legend.position = "none"
  )

print(p_leydig_score)


# ------------------------------------------------------------
# Paired test: Leydig steroidogenic score
# ------------------------------------------------------------

leydig_test <- leydig_mouse %>%
  select(
    replicate,
    fertility_status,
    androgen_synthesis_score
  ) %>%
  pivot_wider(
    names_from = fertility_status,
    values_from = androgen_synthesis_score
  ) %>%
  filter(
    !is.na(Fertile),
    !is.na(Infertile)
  )

cat("\nPaired Wilcoxon test: Leydig steroidogenesis score\n")

print(
  wilcox.test(
    leydig_test$Infertile,
    leydig_test$Fertile,
    paired = TRUE,
    exact = FALSE
  )
)


# ============================================================
# 3. INDIVIDUAL LEYDIG STEROIDOGENIC GENES
# ============================================================

leydig_expr <- FetchData(
  leydig,
  vars = c(
    "replicate",
    "fertility_status",
    androgen_synthesis_genes
  )
)

leydig_gene_mouse <- leydig_expr %>%
  group_by(replicate, fertility_status) %>%
  summarise(
    across(
      all_of(androgen_synthesis_genes),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(androgen_synthesis_genes),
    names_to = "gene",
    values_to = "expression"
  )

p_leydig_genes <- ggplot(
  leydig_gene_mouse,
  aes(
    x = fertility_status,
    y = expression,
    group = replicate
  )
) +
  geom_line(alpha = 0.35) +
  geom_point(
    aes(fill = fertility_status),
    shape = 21,
    size = 2
  ) +
  scale_fill_manual(values = fertility_colors)+
  facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 4
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean normalized expression"
  ) +
  theme(
    legend.position = "none"
  )

print(p_leydig_genes)


# ============================================================
# 4. SERTOLI CELLS: AR EXPRESSION
# ============================================================

sertoli <- subset(
  androgen_obj,
  subset = cell.type3 == "Sertoli"
)

sertoli_ar <- FetchData(
  sertoli,
  vars = c(
    "replicate",
    "fertility_status",
    "Ar"
  )
) %>%
  group_by(replicate, fertility_status) %>%
  summarise(
    Ar = mean(Ar, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

print(sertoli_ar)

p_sertoli_ar <- ggplot(
  sertoli_ar,
  aes(
    x = fertility_status,
    y = Ar,
    group = replicate
  )
) +
  geom_line(alpha = 0.45) +
  geom_point(
    aes(fill = fertility_status),
    shape = 21,
    size = 3,
    stroke = 0.7
  ) +
  scale_fill_manual(values = fertility_colors)+
  theme_classic() +
  labs(
    x = NULL,
    y = "Sertoli Ar expression"
  ) +
  theme(
    legend.position = "none"
  )

print(p_sertoli_ar)


# ============================================================
# 5. SERTOLI CELLS: ANDROGEN-RESPONSE PROGRAM
# ============================================================

# Get Hallmark Androgen Response gene set.
# tryCatch allows compatibility with old/new msigdbr syntax.

hallmark_mouse <- tryCatch(
  msigdbr(
    species = "Mus musculus",
    collection = "H"
  ),
  error = function(e) {
    msigdbr(
      species = "Mus musculus",
      category = "H"
    )
  }
)

androgen_response_genes <- hallmark_mouse %>%
  filter(
    gs_name == "HALLMARK_ANDROGEN_RESPONSE"
  ) %>%
  pull(gene_symbol) %>%
  unique()

androgen_response_genes <- intersect(
  androgen_response_genes,
  rownames(sertoli)
)

cat(
  "\nNumber of Hallmark Androgen Response genes used:",
  length(androgen_response_genes),
  "\n"
)

# Cell-level androgen response score
sertoli <- AddModuleScore(
  sertoli,
  features = list(androgen_response_genes),
  name = "AndrogenResponse"
)

# Mouse/testis-level values
sertoli_mouse <- sertoli@meta.data %>%
  group_by(replicate, fertility_status) %>%
  summarise(
    androgen_response_score = mean(
      AndrogenResponse1,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

print(sertoli_mouse)


# ------------------------------------------------------------
# Paired plot: Sertoli androgen-response score
# ------------------------------------------------------------

p_sertoli_response <- ggplot(
  sertoli_mouse,
  aes(
    x = fertility_status,
    y = androgen_response_score,
    group = replicate
  )
) +
  geom_line(alpha = 0.45) +
  geom_point(
    aes(fill = fertility_status),
    shape = 21,
    size = 3,
    stroke = 0.7
  ) +
  scale_fill_manual(values = fertility_colors)+
  theme_classic() +
  labs(
    x = NULL,
    y = "Sertoli androgen-response score"
  ) +
  theme(
    legend.position = "none"
  )

print(p_sertoli_response)


# ------------------------------------------------------------
# Paired test: Sertoli androgen-response score
# ------------------------------------------------------------

sertoli_test <- sertoli_mouse %>%
  select(
    replicate,
    fertility_status,
    androgen_response_score
  ) %>%
  pivot_wider(
    names_from = fertility_status,
    values_from = androgen_response_score
  ) %>%
  filter(
    !is.na(Fertile),
    !is.na(Infertile)
  )

cat("\nPaired Wilcoxon test: Sertoli androgen-response score\n")

print(
  wilcox.test(
    sertoli_test$Infertile,
    sertoli_test$Fertile,
    paired = TRUE,
    exact = FALSE
  )
)


# ============================================================
# 6. SENSITIVITY ANALYSIS
# Require >=20 Sertoli cells in BOTH testes from a mouse
# ============================================================

sertoli_mouse_20 <- sertoli_mouse %>%
  filter(n_cells >= 20) %>%
  group_by(replicate) %>%
  filter(n() == 2) %>%
  ungroup()

print(sertoli_mouse_20)

sertoli_test_20 <- sertoli_mouse_20 %>%
  select(
    replicate,
    fertility_status,
    androgen_response_score
  ) %>%
  pivot_wider(
    names_from = fertility_status,
    values_from = androgen_response_score
  ) %>%
  filter(
    !is.na(Fertile),
    !is.na(Infertile)
  )

cat(
  "\nPaired Wilcoxon test: Sertoli androgen response",
  "(>=20 cells per condition)\n"
)

print(
  wilcox.test(
    sertoli_test_20$Infertile,
    sertoli_test_20$Fertile,
    paired = TRUE,
    exact = FALSE
  )
)







#make a euler diagram

# ============================================================
# 7. OVERLAP OF SERTOLI DEGs WITH ANDROGEN-RESPONSE GENES
# ============================================================




# ------------------------------------------------------------
# 7A. Differential expression:
#     Experimental infertile vs fertile Sertoli cells
#
# Positive log2FC = higher in infertile
# Negative log2FC = higher in fertile
# ------------------------------------------------------------

sertoli_de <- sertoli

DefaultAssay(sertoli_de) <- "RNA"

# sertoli_de <- JoinLayers(
#   sertoli_de,
#   assay = "RNA"
# )

Idents(sertoli_de) <- "fertility_status"

sertoli_markers <- FindMarkers(
  sertoli_de,
  ident.1 = "Infertile",
  ident.2 = "Fertile",
  test.use = "LR",
  latent.vars = "nCount_RNA",
  logfc.threshold = 0,
  min.pct = 0.1
)

sertoli_markers$gene <- rownames(sertoli_markers)


# ------------------------------------------------------------
# 7B. Define DEGs
#
# Change log2fc_cutoff here if you want a stricter definition.
# 0.25 is useful for the broad DEG set.
# ------------------------------------------------------------

log2fc_cutoff <- 1
padj_cutoff <- 0.05

sertoli_deg_genes <- sertoli_markers %>%
  dplyr::filter(
    p_val_adj < padj_cutoff,
    abs(avg_log2FC) >= log2fc_cutoff
  ) %>%
  dplyr::pull(gene) %>%
  unique()


# ------------------------------------------------------------
# 7C. Make sure both sets are restricted to genes actually
#     present/testable in the Sertoli dataset
# ------------------------------------------------------------

androgen_genes_testable <- intersect(
  androgen_response_genes,
  rownames(sertoli_de)
)

sertoli_deg_genes <- intersect(
  sertoli_deg_genes,
  rownames(sertoli_de)
)


# ------------------------------------------------------------
# 7D. Calculate overlap
# ------------------------------------------------------------

overlap_genes <- intersect(
  sertoli_deg_genes,
  androgen_genes_testable
)

n_deg <- length(sertoli_deg_genes)
n_androgen <- length(androgen_genes_testable)
n_overlap <- length(overlap_genes)

percent_deg_androgen <- 100 * n_overlap / n_deg

cat("\n========================================\n")
cat("SERTOLI DEG / ANDROGEN RESPONSE OVERLAP\n")
cat("========================================\n")

cat("Sertoli DEGs:", n_deg, "\n")
cat("Androgen-response genes:", n_androgen, "\n")
cat("Overlap:", n_overlap, "\n")

cat(
  "Percent of Sertoli DEGs in androgen-response set:",
  round(percent_deg_androgen, 2),
  "%\n"
)

cat("\nOverlapping genes:\n")
print(sort(overlap_genes))


# ------------------------------------------------------------
# 7E. Euler plot
# ------------------------------------------------------------

fit <- eulerr::euler(
  list(
    "Sertoli DEGs" = sertoli_deg_genes,
    "Androgen response" = androgen_genes_testable
  )
)

plot(
  fit,
  fills = c(
    "#b12625",
    "#48b1e0"
  ),
  edges = "black",
  quantities = TRUE
)





########################################################################
########################################################################
########################################################################
########################################################################
### LIANA ###
########################################################################
########################################################################
########################################################################
########################################################################



setwd("/Users/ewhelan/Library/CloudStorage/Box-Box/Somatic_Multiome_Project")
# Load the data from the Excel file

df <- read_csv("LIANA_CellType2_leydig.csv")
df <- read_csv("LIANA_CellType2_myoid.csv")
df <- read_csv("LIANA_CellType2_sertoli.csv") #this one is buggered, fix it first



target_groups <- c("SSCs", "Progenitors", "DiffSpermatogonia", "PrelepSpermatocytes", "EarlySpermatocytes")

# Subset the rows where the target matches the specified values
subset_df <- df %>% filter(target %in% target_groups)

# Extract unique ligand.complex values
unique_ligands <- unique(subset_df$ligand.complex)

# View the unique ligand.complex results
print(unique_ligands)

# If you want to save the unique ligand.complex values back to a CSV file
write.csv(unique_ligands, "unique_ligands.csv", row.names = FALSE)








# Get all unique ligand/receptor combinations
unique_combinations <- df %>%
  select(ligand.complex, receptor.complex) %>%
  distinct()

# Get all unique cell types
cell_types <- unique(df$target)

# Initialize a new dataframe with unique ligand/receptor combinations
result <- drop_na(unique_combinations)

# Add a column for each cell type, initialized to FALSE
for (cell_type in cell_types) {
  result[[cell_type]] <- FALSE
}

# Fill in the TRUE values where the combinations are found in the cell types
for (i in 1:nrow(result)) {
  for (cell_type in cell_types) {
    if (any(df$ligand.complex == result$ligand.complex[i] & df$receptor.complex == result$receptor.complex[i] & df$target == cell_type)) {
      result[[cell_type]][i] <- TRUE
    }
  }
}

result_tibble <- result %>%
  select(ligand.complex, receptor.complex, SSCs, Progenitors, DiffSpermatogonia, PrelepSpermatocytes, EarlySpermatocytes, everything())

length(unique(result_tibble$receptor.complex))

# View the result
print(result_tibble)




write.csv(result_tibble, "myoid_output.csv")

setwd("/Users/ewhelan/Desktop/Somatic Project/somatic_samples")
load("somatic.integrated.germ4.Robj")
genes.of.interest <- result_tibble$receptor.complex
genes.of.interest <- c("Fgfr1", "Aplp1", "Aplp2", "Lrp10", "Ncstn", "Notch2", "Lrp6",
                       "Rpsa", "Lrp1", "Ar", "Sort1", "Itga1", "Itgb1", "Itga9", "Sdc1",
                       "Ptch1", "Cdon", "Ryr2", "Itgav", "Itgb5", "Tgfbr3", "Acvr1", "Itgav")
# genes.of.interest <- c("Sdc4", "Itga9", "Itgb1", "Bmpr1a", "Bmpr1b", "Bmpr2",
#                        "Robo1", "Scd1", "Cd93", "Cd47", "Itga1", "Itga9", "Itga6", "Dag1")
DotPlot(object = somatic.integrated.germ4, features = unique(genes.of.interest), cols = c("gray", "#890600")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16, face = "italic"),   # Increase x-axis label size and rotate
        axis.text.y = element_text(size = 14),                    # Increase y-axis label size and italicize
        axis.title = element_text(size = 12)) #+                                    # Increase axis title size
# scale_y_discrete(position = "right")


setwd("/Users/ewhelan/Documents/PROJECTS/Somatic cell analysis/")
load("somatic.integrated.somatic6.Robj")


DE.genes <- c("Esr1", "Rora", "Nfia", "Essrg",
              "Lgals1", "Adam12", "Bmp1", "Slit3", "Col4a5", "Lama4")

Idents(somatic.integrated.somatic6) <- "cell.type3"
DimPlot(somatic.integrated.somatic6)
DefaultAssay(somatic.integrated.somatic6) <- "RNA"

subset_Leydig <- subset(somatic.integrated.somatic6, idents = "Leydig_cells")
subset_T_immunomodulatory <- subset(somatic.integrated.somatic6, idents = "T_cells_immunomodulatory")
subset_T_effector <- subset(somatic.integrated.somatic6, idents = "T_cells_effector")
subset_MacrophageInter <- subset(somatic.integrated.somatic6, idents = "Macrophages_interstitial")
subset_MacrophagePeri <- subset(somatic.integrated.somatic6, idents = "Macrophages_peritubular")
subset_Dendritic <- subset(somatic.integrated.somatic6, idents = "Dendritic_cells")
subset_Epithelial <- subset(somatic.integrated.somatic6, idents = "Epithelial_cells")
subset_B <- subset(somatic.integrated.somatic6, idents = "B_cells")
subset_Mesenchymal <- subset(somatic.integrated.somatic6, idents = "Mesenchymal_progenitors")
subset_Myoid_peritubular <- subset(somatic.integrated.somatic6, idents = "Peritubular_myoid")
subset_Myoid_perivascular <- subset(somatic.integrated.somatic6, idents = "Perivascular_smooth_muscle")
subset_Sertoli <- subset(somatic.integrated.somatic6, idents = "Sertoli")
subset_Endothelial <- subset(somatic.integrated.somatic6, idents = "Endothelial")

Idents(subset_Leydig) <- "treatment"
Idents(subset_T_immunomodulatory) <- "treatment" 
Idents(subset_T_effector) <- "treatment" 
Idents(subset_MacrophageInter) <- "treatment" 
Idents(subset_MacrophagePeri) <- "treatment" 
Idents(subset_Dendritic) <- "treatment" 
Idents(subset_Epithelial) <- "treatment" 
Idents(subset_B) <- "treatment"
Idents(subset_Mesenchymal) <- "treatment" 
Idents(subset_Myoid_peritubular) <- "treatment" 
Idents(subset_Myoid_perivascular) <- "treatment" 
Idents(subset_Sertoli) <- "treatment"
Idents(subset_Endothelial) <- "treatment"


DE.genes <- unique(c(df$TF, df$ligand.complex))


markers.Leydig <- FindMarkers(subset_Leydig, features = DE.genes, logfc.threshold = 0, min.pct = 0,
                              ident.1 = "GFP transplanted", ident.2 = "control")
# markers.T_immunomodulatory <- FindMarkers(subset_T_immunomodulatory, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.T_effector <- FindMarkers(subset_T_effector, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.MacrophageInter <- FindMarkers(subset_MacrophageInter, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.MacrophagePeri <- FindMarkers(subset_MacrophagePeri, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.Dendritic <- FindMarkers(subset_Dendritic, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.Epithelial <- FindMarkers(subset_Epithelial, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.B <- FindMarkers(subset_B, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.Mesenchymal <- FindMarkers(subset_Mesenchymal, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.Myoid_peritubular <- FindMarkers(subset_Myoid_peritubular, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.Myoid_perivascular <- FindMarkers(subset_Myoid_perivascular, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.Sertoli <- FindMarkers(subset_Sertoli, ident.1 = "GFP transplanted", ident.2 = "control")
# markers.Endothelial <- FindMarkers(subset_Endothelial, ident.1 = "GFP transplanted", ident.2 = "control")



leydig.curated.GRN.results <- read.csv("/Users/ewhelan/Desktop/Somatic Project/somatic_SCENIC/Leydig_curated.csv")



# Read the CSV file
df <- leydig.curated.GRN.results



# Invert the sign of the values in Gene_log2FC and TF_log2FC
df$Gene_log2FC <- -df$Gene_log2FC
df$TF_log2FC <- -df$TF_log2FC

# Create the Gene heatmap
gene_heatmap_data <- df %>% 
  select(Gene, Gene_log2FC) %>%
  distinct()

gene_matrix <- as.matrix(gene_heatmap_data[,2])
rownames(gene_matrix) <- gene_heatmap_data$Gene

# Set color breaks for the Gene heatmap
gene_breaks <- seq(-0.9, 0.9, length.out = 100)

pheatmap(gene_matrix, cluster_rows = TRUE, cluster_cols = FALSE, 
         display_numbers = TRUE, main = "Gene Heatmap (Inverted log2FC)",
         breaks = gene_breaks, color = colorRampPalette(c("blue", "white", "red"))(100))

# Create the TF heatmap
tf_heatmap_data <- df %>% 
  select(TF, TF_log2FC) %>%
  distinct()

tf_matrix <- as.matrix(tf_heatmap_data[,2])
rownames(tf_matrix) <- tf_heatmap_data$TF

# Set color breaks for the TF heatmap
# tf_breaks <- seq(-0.9, 0.9, length.out = 100)

pheatmap(tf_matrix, cluster_rows = TRUE, cluster_cols = FALSE, 
         display_numbers = TRUE, main = "TF Heatmap (Inverted log2FC)",
         breaks = gene_breaks, color = colorRampPalette(c("blue", "white", "red"))(100))




# Assuming df is your data frame
# Aggregate the sca.LRscore by ligand.complex and receptor.complex
df_aggregated <- df %>%
  group_by(ligand.complex, receptor.complex) %>%
  summarise(sca.LRscore = mean(sca.LRscore, na.rm = TRUE))

# Create a unique identifier for each ligand-receptor pair
df_aggregated <- df_aggregated %>%
  mutate(pair = paste(ligand.complex, receptor.complex, sep = "_"))

# Convert to matrix format for heatmap
heatmap_matrix <- as.matrix(df_aggregated$sca.LRscore)
rownames(heatmap_matrix) <- df_aggregated$pair

# Replace NA values with zero
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Define the breaks for the color scale, ensuring 0 is mapped to a light color
max_val <- max(heatmap_matrix, na.rm = TRUE)
min_val <- min(heatmap_matrix, na.rm = TRUE)
breaks <- c(seq(min_val, 0, length.out = 50), seq(0, max_val, length.out = 51)[-1])

# Define the color palette
color_palette <- colorRampPalette(c("white", "lightyellow", "lightyellow", "yellow" ,"green3"))(length(breaks) - 1)

# Plot the heatmap with the adjusted color scale
pheatmap(heatmap_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # No need to cluster columns as there's only one column
         display_numbers = TRUE, 
         main = "Heatmap of LRscore by Ligand-Receptor Complexes", 
         color = color_palette,
         breaks = breaks)






###



# Ensure ligand and ligand.complex are character type
df$ligand <- as.character(df$ligand)
df$ligand.complex <- as.character(df$ligand.complex)

# Filter and extract unique TFs
unique_TFs <- df %>%
  filter(ligand %in% ligand.complex) %>%
  pull(TF) %>%
  unique()

print(unique_TFs)



FeaturePlot(somatic.integrated.somatic6, features = "Aldh1a1")
VlnPlot(somatic.integrated.somatic6, features = "Aldh1a1", split.by = "treatment")
VlnPlot(somatic.integrated.somatic6, features = "Aldh1a1")


tf_combinations <- read_csv("myoid_TF_combinations.csv")

all.ligands <- unique_ligands


subset_tf_combinations <- tf_combinations %>% filter(Gene %in% all.ligands)

# Get unique combinations
unique_combinations <- drop_na(subset_tf_combinations %>% distinct())

# View the unique combinations
print(unique_combinations)


########################################################################
########################################################################
########################################################################
########################################################################
### CHIPSEEKER ###
########################################################################
########################################################################
########################################################################
########################################################################




#ChIPseeker

library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(dplyr)
library(tidyverse)
library(GenomicRanges)
# all.markers.peaks.mouse



DefaultAssay(somatic.integrated.somatic7_ATAConly) <- "peaks"
Idents(somatic.integrated.somatic7_ATAConly) <- "cell.type3"
DimPlot(somatic.integrated.somatic7_ATAConly)

Idents(somatic.integrated.somatic7_ATAConly) <- "treatment"
table(Idents(somatic.integrated.somatic7_ATAConly))
somatic.integrated.somatic7_ATAConly <- RenameIdents(somatic.integrated.somatic7_ATAConly, 
                                                     'GFP transplanted' = "fertile",
                                                     'control' = "infertile")
somatic.integrated.somatic7_ATAConly$fertility <- Idents(somatic.integrated.somatic7_ATAConly)

somatic.integrated.somatic7_ATAConly$cell.type3 <- factor(somatic.integrated.somatic7_ATAConly$cell.type3, levels = c(
  "Sertoli", 
  "Leydig_cells",
  "Peritubular_myoid",
  "Perivascular_smooth_muscle",
  "Mesenchymal_progenitors",
  "Epithelial_cells",
  "Endothelial", 
  "Macrophages_peritubular",
  "Macrophages_interstitial",
  "Dendritic_cells",
  "T_cells_immunomodulatory",
  "T_cells_effector",
  "B_cells"
))

# list.of.cell.types <- levels(Idents(somatic.integrated.somatic7_ATAConly))
list.of.cell.types <- c("Sertoli", 
                        "Epithelial_cells",
                        "Leydig_cells",
                        "Peritubular_myoid",
                        "Perivascular_smooth_muscle",
                        "Mesenchymal_progenitors",
                        
                        "Endothelial", 
                        "Macrophages_peritubular",
                        "Macrophages_interstitial",
                        "Dendritic_cells",
                        "T_cells_immunomodulatory",
                        "T_cells_effector",
                        "B_cells")

all.markers.peaks.mouse <- FindAllMarkers(somatic.integrated.somatic7_ATAConly)



iteration <- 1
# Loop through each cell type in list.of.cell.types
for (cell_type in list.of.cell.types) {
  # cell_type <- "Leydig_cells"
  Idents(somatic.integrated.somatic7_ATAConly) <- "cell.type3"
  current.seurat <- subset(somatic.integrated.somatic7_ATAConly, idents = cell_type)
  Idents(current.seurat) <- "fertility"
  all.markers.fertility <- FindMarkers(current.seurat, ident.1 = "infertile")
  all.markers.fertility$region <- rownames(all.markers.fertility)
  
  var_name <- paste0("Granges_", gsub("[ /]", "", cell_type))
  
  # Filter all.markers.peaks.mouse for the current cell type
  top_peaks <- all.markers.peaks.mouse %>%
    dplyr::filter(cluster == cell_type) %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::arrange(desc(avg_log2FC)) 
  
  top_peaks_fertility <- all.markers.fertility %>%
    # dplyr::filter(cluster == cell_type) %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::arrange(desc(avg_log2FC)) 
  
  peak_data <- top_peaks %>%
    separate(gene, into = c("chr", "start", "end"), sep = "-")
  peak_data_fertility <- top_peaks_fertility %>%
    separate(region, into = c("chr", "start", "end"), sep = "-")
  
  granges_object <- GRanges(
    seqnames = peak_data$chr,
    ranges = IRanges(start = as.numeric(peak_data$start), end = as.numeric(peak_data$end)),
    p_val = peak_data$p_val,
    avg_log2FC = peak_data$avg_log2FC,
    cluster = peak_data$cluster
  )
  
  granges_object_fertility <- GRanges(
    seqnames = peak_data_fertility$chr,
    ranges = IRanges(start = as.numeric(peak_data_fertility$start), end = as.numeric(peak_data_fertility$end)),
    p_val = peak_data_fertility$p_val,
    avg_log2FC = peak_data_fertility$avg_log2FC,
    cluster = cell_type
  )
  
  library(org.Mm.eg.db)
  peakAnno.edb <- annotatePeak(peak = granges_object, tssRegion=c(-3000, 3000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno.edb.fertility <- annotatePeak(peak = granges_object_fertility, tssRegion=c(-3000, 3000),
                                         TxDb=txdb, annoDb="org.Mm.eg.db")
  # pie_plot <- plotAnnoPie(peakAnno.edb)
  # plotAnnoBar(peakAnno.edb)
  # plotDistToTSS(peakAnno.edb)
  var_name2 <- paste0("peakAnno_", gsub("[ /]", "", cell_type))
  
  annotation_summary <- as.data.frame(peakAnno.edb@annoStat)
  annotation_summary_fertility <- as.data.frame(peakAnno.edb.fertility@annoStat)
  
  # Assign the filtered peaks to the dynamically created variable
  assign(var_name, granges_object)
  assign(var_name2, peakAnno.edb)
  
  # rownames(annotation_summary) <- annotation_summary[,1]
  # annotation_summary2 <- as.data.frame(annotation_summary[,2])
  colnames(annotation_summary) <- c("feature", cell_type)
  rownames(annotation_summary) <- annotation_summary[,1]
  
  colnames(annotation_summary_fertility) <- c("feature", cell_type)
  rownames(annotation_summary_fertility) <- annotation_summary[,1]
  
  if(iteration == 1){
    annotation_summary_all <- annotation_summary
    annotation_summary_fertility_all <- annotation_summary_fertility
  }else{
    annotation_summary_all <- full_join(annotation_summary_all, annotation_summary, by = "feature")
    annotation_summary_fertility_all <- full_join(annotation_summary_fertility_all, annotation_summary_fertility, by = "feature")
  }
  
  
  iteration <- iteration + 1
  
}

annotation_summary_all
annotation_summary_fertility_all

# Load necessary libraries
library(ggplot2)
library(tidyr)

# Assuming your data is stored in `annotation_summary_all`
# Convert to long format
long_data <- pivot_longer(annotation_summary_all, 
                          cols = -feature, 
                          names_to = "Cell_Type", 
                          values_to = "Percentage")

long_data_fertility <- pivot_longer(annotation_summary_fertility_all, 
                                    cols = -feature, 
                                    names_to = "Cell_Type", 
                                    values_to = "Percentage")

# Ensure 'Cell_Type' is a factor with levels in the desired order
long_data$Cell_Type <- factor(long_data$Cell_Type, levels = list.of.cell.types)
long_data_fertility$Cell_Type <- factor(long_data_fertility$Cell_Type, levels = list.of.cell.types)

feature_order <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "3' UTR", "5' UTR", "1st Intron", "Other Intron", "1st Exon", "Other Exon", "Downstream (<=300)", "Distal Intergenic")


long_data <- long_data %>%
  mutate(feature = factor(feature, levels = feature_order))
long_data_fertility <- long_data_fertility %>%
  mutate(feature = factor(feature, levels = feature_order))

# Create the stacked bar plot with the specified order
stacked.bar.plot.mouse <- ggplot(long_data, aes(x = Cell_Type, y = Percentage, fill = feature)) +
  geom_bar(stat = "identity") +
  labs(title = "Mouse", 
       x = "Cell Type", 
       y = "Percentage of Features by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

stacked.bar.plot.mouse

stacked.bar.plot.mouse_fertility <- ggplot(long_data_fertility, aes(x = Cell_Type, y = Percentage, fill = feature)) +
  geom_bar(stat = "identity") +
  labs(title = "Infertile vs Fertile", 
       x = "Cell Type", 
       y = "Percentage of Features by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

stacked.bar.plot.mouse_fertility


# plot_grid(p1, p2, p3, p4, p5, p6, p7, 
#           p8, p9, p10, p11, p12, p13, p14)

peak.combined=GenomicRanges::GRangesList(SSCs=Granges_SSCs,
                                         undiff=Granges_Undiff.spermatogonia,
                                         earlydiff=Granges_Earlydiff.spermatogonia,
                                         latediff=Granges_Latediff.spermatogonia,
                                         prelep=Granges_Preleptotenespermatocytes,
                                         lepzyg=Granges_LepZygspermatocytes,
                                         pachy=Granges_Pachytenespermatocytes,
                                         dipsecondary=Granges_Dip2ndspermatocytes,
                                         earlyround=Granges_Earlyroundspermatids,
                                         midround=Granges_Midroundspermatids,
                                         lateround=Granges_Lateroundspermatids,
                                         earlyelongating=Granges_Earlyelongatingspermatids,
                                         midelongating=Granges_Midelongatingspermatids,
                                         lateelongating=Granges_Lateelongatingspermatids)



# View the GRanges object
# granges_object
table(Idents(mouse.germ.cells.multiome))
library(paletteer)
colors <- paletteer_d("ggthemes::manyeys", n = 14)
colors <- paletteer_c("grDevices::topo.colors", 14)
library(ggplot2)
col <- c(SSCs=colors[1], "Undiff. spermatogonia"=colors[2], "Early diff. spermatogonia"=colors[3],
         "Late diff. spermatogonia"=colors[4],"Preleptotene spermatocytes"=colors[5],"Lep/Zyg spermatocytes"=colors[6],
         "Pachytene spermatocytes"=colors[7],"Dip/2nd spermatocyte"=colors[8],"Early round spermatids"=colors[9],
         "Mid round spermatids"=colors[10],"Late round spermatids"=colors[11],"Early elongating spermatids"=colors[12],
         "Mid elongating spermatids"=colors[13],"Late elongating spermatids"=colors[14])

cov <- covplot(peak.combined, weightCol = "avg_log2FC",  title = "Mouse") 
# cov <- covplot(peak)
# cov + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
cov_mouse <- cov + 
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 15) + 
  theme(
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    plot.title = element_text(color = "white"),
    panel.grid = element_blank() 
  )






