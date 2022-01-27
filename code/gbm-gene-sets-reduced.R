
# Load R packages ---------------------------------------------------------

library(tidyverse)
library(fgsea)



# Get data ----------------------------------------------------------------

# Load GBM gene sets we already have
gbm <- readRDS(url("https://github.com/WalterMuskovic/gene_sets_and_gbm_markers/blob/master/data/MSigDB_and_gbm_gene_lists.rds?raw=true"))
names(gbm) <- c("all_gene_sets", "hallmark_gene_sets", "positional_gene_sets",
                     "perturbations_gene_sets", "canonical_pathways_gene_sets", "biocarta_gene_sets",
                     "kegg_gene_sets", "pid_gene_sets","reactome_gene_sets",
                     "tf_gene_sets", "cancer_gene_gene_sets", "cancer_module_gene_sets",
                     "go_bp_gene_sets", "oncogenic_gene_sets", "immunologic_gene_sets",
                     "gbm_lists")
gbm <- gbm[["gbm_lists"]]

# Download gene sets from MSigDB - downloaded 27/01/2022
hallmark <- gmtPathways("~/Downloads/h.all.v7.5.1.symbols.gmt")
cell.types <- gmtPathways("~/Downloads/c8.all.v7.5.1.symbols.gmt")
C5.GOBP <- gmtPathways("~/Downloads/c5.go.bp.v7.5.1.symbols.gmt")
C5.GOCC <- gmtPathways("~/Downloads/c5.go.cc.v7.5.1.symbols.gmt")
C5.GOMF <- gmtPathways("~/Downloads/c5.go.mf.v7.5.1.symbols.gmt")
pid <- gmtPathways("~/Downloads/c2.cp.pid.v7.5.1.symbols.gmt")


# Save data ---------------------------------------------------------------

gbm <- list(
  gbm = gbm,
  hallmark = hallmark,
  cell_types = cell.types,
  go_biological_processes = C5.GOBP,
  go_cellular_components = C5.GOCC,
  go_molecular_functions = C5.GOMF,
  pathway_interaction_database = pid
  )
saveRDS(gbm, "data/MSigDB_and_gbm_gene_lists_subset.rds")

