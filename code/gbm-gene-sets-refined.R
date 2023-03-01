
# Script description ------------------------------------------------------

# In this script we build on our previous work gathering gene lists to create a
# refined list, with gene sets relevant to the biology of glioma cell states



# Load R packages ---------------------------------------------------------

library(tidyverse)
library(fgsea)



# Load data ---------------------------------------------------------------

# Load GBM gene sets we already have
gbm <- readRDS(url("https://github.com/WalterMuskovic/gene_sets_and_gbm_markers/blob/master/data/MSigDB_and_gbm_gene_lists.rds?raw=true"))
names(gbm) <- c("all_gene_sets", "hallmark_gene_sets", "positional_gene_sets",
                "perturbations_gene_sets", "canonical_pathways_gene_sets", "biocarta_gene_sets",
                "kegg_gene_sets", "pid_gene_sets","reactome_gene_sets",
                "tf_gene_sets", "cancer_gene_gene_sets", "cancer_module_gene_sets",
                "go_bp_gene_sets", "oncogenic_gene_sets", "immunologic_gene_sets",
                "gbm_lists")
gbm <- gbm[["gbm_lists"]]

#Import .gmt files downloaded from MSigDB:

# Hallmark gene sets - Coherently expressed signatures derived by aggregating
# many MSigDB gene sets to represent well-defined biological states or processes.
hallmark_gene_sets <- gmtPathways("data/h.all.v2022.1.Hs.symbols.gmt")

# c2: curated gene sets - from online pathway databases, publications in
# PubMed, and knowledge of domain experts;
# chemical and genetic perturbations
perturbations_gene_sets <- gmtPathways("data/c2.cgp.v2022.1.Hs.symbols.gmt")
# all canonical pathways
canonical_pathways_gene_sets <- gmtPathways("data/c2.cp.v2022.1.Hs.symbols.gmt")
# KEGG
kegg_gene_sets <- gmtPathways("data/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
# reactome
reactome_gene_sets <- gmtPathways("data/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

# c4: computational gene sets - defined by mining large collections of
# cancer-oriented microarray data.
# cancer gene neighborhoods
cancer_gene_gene_sets <- gmtPathways("data/c4.cgn.v2022.1.Hs.symbols.gmt")

# c5: GO gene sets -  genes annotated by the same GO terms.
# biological processes
go_bp_gene_sets <- gmtPathways("data/c5.go.bp.v2022.1.Hs.symbols.gmt")
go_cc_gene_sets <- gmtPathways("data/c5.go.cc.v2022.1.Hs.symbols.gmt")
go_mf_gene_sets <- gmtPathways("data/c5.go.mf.v2022.1.Hs.symbols.gmt")


# c6: oncogenic gene sets - defined directly from microarray gene
# expression data from cancer gene perturbations.
oncogenic_gene_sets <- gmtPathways("data/c6.all.v2022.1.Hs.symbols.gmt")

# C8: cell type signature gene sets - Gene sets that contain curated cluster
# markers for cell types identified in single-cell sequencing studies of human
# tissue.
cell_type_gene_sets <- gmtPathways("data/c8.all.v2022.1.Hs.symbols.gmt")



# Save data ---------------------------------------------------------------

gbm <- list(
  gbm = gbm,
  hallmark = hallmark_gene_sets,
  perturbations = perturbations_gene_sets,
  canonical_pathways = canonical_pathways_gene_sets,
  kegg = kegg_gene_sets,
  reactome = reactome_gene_sets,
  cancer_gene = cancer_gene_gene_sets,
  go_biological_processes = go_bp_gene_sets,
  go_cellular_components = go_cc_gene_sets,
  go_molecular_functions = go_mf_gene_sets,
  oncogenic = oncogenic_gene_sets,
  cell_type = cell_type_gene_sets
)
saveRDS(gbm, "data/MSigDB_and_gbm_gene_lists_refined.rds")

