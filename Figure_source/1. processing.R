# .h5data can be load from https://doi.org/10.5281/zenodo.18308841
# Packages must be installed before running this script (see readme for versions).
# If needed, run once: install.packages("Seurat")
library(Seurat)
library(hdf5r)

h5_so_SCT <- function(path,project){
data <- Read10X_h5(path)
so <- CreateSeuratObject(data,project = project)
so <- SCTransform(so)
return(so)
}

P3ctrl_n1 <- h5_so_SCT("P3_ctr_n1b1_filtered_feature_bc_matrix.h5",project = "P3_ctr_n1b1")
P3ctrl_n2L1 <- h5_so_SCT("P3_ctr_n1l1_filtered_feature_bc_matrix.h5",project = "P3_ctr_n1l1")
P3ctrl_n2L2 <- h5_so_SCT("P3_ctr_n1l2_filtered_feature_bc_matrix.h5",project = "P3_ctr_n1l2")

P3cKO_n1 <- h5_so_SCT("P3_cKO_n1b1_filtered_feature_bc_matrix.h5",project = "P3_cKO_n1b1")
P3cKO_n2L1 <- h5_so_SCT("P3_cKO_n1l1_filtered_feature_bc_matrix.h5",project = "P3_cKO_n1l1")
P3cKO_n2L2 <- h5_so_SCT("P3_cKO_n1l2_filtered_feature_bc_matrix.h5",project = "P3_cKO_n1l2")

P10ctrl_n1 <- h5_so_SCT("P10_ctr_n1b1_filtered_feature_bc_matrix.h5",project = "P10_ctr_n1b1")
P10ctrl_n2L1 <- h5_so_SCT("P10_ctr_n1l1_filtered_feature_bc_matrix.h5",project = "P10_ctr_n1l1")
P10ctrl_n2L2 <- h5_so_SCT("P10_ctr_n1l2_filtered_feature_bc_matrix.h5",project = "P10_ctr_n1l2")

P10cKO_n1 <- h5_so_SCT("P10_cKO_n1b1_filtered_feature_bc_matrix.h5",project = "P10_cKO_n1b1")
P10cKO_n2L1 <- h5_so_SCT("P10_cKO_n1l1_filtered_feature_bc_matrix.h5",project = "P10_cKO_n1l1")
P10cKO_n2L2<- h5_so_SCT("P10_cKO_n1l2_filtered_feature_bc_matrix.h5",project = "P10_cKO_n1l2")

P3ctrl_n1$stage <- "P3"
P3ctrl_n1$genotype <- "ctr"
P3ctrl_n1$bio_rep <- "P3_ctr_rep1"
P3ctrl_n1$library <- "P3_ctr_n1b1"

P3ctrl_n2L1$stage <- "P3"
P3ctrl_n2L1$genotype <- "ctr"
P3ctrl_n2L1$bio_rep <- "P3_ctr_rep2"
P3ctrl_n2L1$library <- "P3_ctr_n1l1"

P3ctrl_n2L2$stage <- "P3"
P3ctrl_n2L2$genotype <- "ctr"
P3ctrl_n2L2$bio_rep <- "P3_ctr_rep2"
P3ctrl_n2L2$library <- "P3_ctr_n1l2"

P3cKO_n1$stage <- "P3"
P3cKO_n1$genotype <- "cKO"
P3cKO_n1$bio_rep <- "P3_cKO_rep1"
P3cKO_n1$library <- "P3_cKO_n1b1"

P3cKO_n2L1$stage <- "P3"
P3cKO_n2L1$genotype <- "cKO"
P3cKO_n2L1$bio_rep <- "P3_cKO_rep2"
P3cKO_n2L1$library <- "P3_cKO_n1l1"

P3cKO_n2L2$stage <- "P3"
P3cKO_n2L2$genotype <- "cKO"
P3cKO_n2L2$bio_rep <- "P3_cKO_rep2"
P3cKO_n2L2$library <- "P3_cKO_n1l2"

P10ctrl_n1$stage <- "P10"
P10ctrl_n1$genotype <- "ctr"
P10ctrl_n1$bio_rep <- "P10_ctr_rep1"
P10ctrl_n1$library <- "P10_ctr_n1b1"

P10ctrl_n2L1$stage <- "P10"
P10ctrl_n2L1$genotype <- "ctr"
P10ctrl_n2L1$bio_rep <- "P10_ctr_rep2"
P10ctrl_n2L1$library <- "P10_ctr_n1l1"

P10ctrl_n2L2$stage <- "P10"
P10ctrl_n2L2$genotype <- "ctr"
P10ctrl_n2L2$bio_rep <- "P10_ctr_rep2"
P10ctrl_n2L2$library <- "P10_ctr_n1l2"

P10cKO_n1$stage <- "P10"
P10cKO_n1$genotype <- "cKO"
P10cKO_n1$bio_rep <- "P10_cKO_rep1"
P10cKO_n1$library <- "P10_cKO_n1b1"

P10cKO_n2L1$stage <- "P10"
P10cKO_n2L1$genotype <- "cKO"
P10cKO_n2L1$bio_rep <- "P10_cKO_rep2"
P10cKO_n2L1$library <- "P10_cKO_n1l1"

P10cKO_n2L2$stage <- "P10"
P10cKO_n2L2$genotype <- "cKO"
P10cKO_n2L2$bio_rep <- "P10_cKO_rep2"
P10cKO_n2L2$library <- "P10_cKO_n1l2"


so.list <- list(P3ctrl_n1,P3ctrl_n2L1,P3ctrl_n2L2,P3cKO_n1,P3cKO_n2L1,P3cKO_n2L2,
                P10ctrl_n1,P10ctrl_n2L1,P10ctrl_n2L2,P10cKO_n1,P10cKO_n2L1,P10cKO_n2L2)

features <- SelectIntegrationFeatures(object.list = so.list)
anchors <- FindIntegrationAnchors(so.list,anchor.features = features)
so.combined <- IntegrateData(anchorset = anchors)
