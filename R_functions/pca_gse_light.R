# Custom Principal component gene set enrichment analysis.
#  Pipeline for identifying gene-sets that explain major PCs and may reflect technical variance


# Inputs:
#  - expr: matrix or data.frame, rows = genes (symbols or Ensembl), cols = samples (libraries)
#  - rownames(expr) should be gene symbols OR Ensembl IDs; see mapping step below
#
# Outputs:
#  - camera_res_pr: list object with PCGSE results


# Notes: 
# Currently pulls the Hallmark, Gene Ontology Biological Process, and HOUNKPE_HOUSEKEEPING_GENES genesets. This could be
#  made more flexible in the future, or internals can be manually tweaked if desired. 


#  There are a few existing approaches to this, including PCGSE (Frost et al.) and 
#  GO-PCA , though here we run similar analysis denovo using limma's competetive gene-set enrichment "camera".

#  We include user-defined technical gene-sets (mitochondrial, ribosomal, histones, hemoglobin, housekeeping) in order
# to locate PCs that likely reflect technical covariance. 

#  In practice most of these are under-represented in the data and so only a few are useful.

#  - Adjust gene set sources (MSigDB categories or GO) as needed.

# # ---------------------------
# # Load packages
# # ---------------------------
# 
# require(PCGSE) 
# require(msigdbr)  
# require(biomaRt)  
# require(tidyverse)
# require(pheatmap)
# require(matrixStats)
# require(rrvgo)
# require(GOSemSim)
# require(org.Hs.eg.db)
# require(GO.db)
# require(GSEABase)

pca_gse_light<-function(
  expr=NULL,        # input normalized expression matrix with genes/features as rows and samples/libraries as columns. rownames must be ensembl gene ids or HGNC symbols
  gene_key =NULL,   # input dataframe defining genes in expr. must include a column "symbol" denoting HGNC symbols
  norm_data=NULL,     # optionally an EList object which includes normalized expression matrix "E" and gene key "genes" which conform to specifications above.
  min_geneset_size = 10,   # min gene-set size to include
  max_geneset_size = 2000, # max gene-set size to include
  species = "Homo sapiens", # Species for msigdbr. Default is human.
  nPC=NULL, # The number of PCs to include in the competetive gene-set test. Default is any PCs explaining 1% or greater of variance.
  topN=250 # The number of top GO gene-sets to retain for each PC prior to reduction based on semantic similarity.
  ){
  
  # ---------------------------
  # Load packages
  # ---------------------------
  
  require(PCGSE) 
  require(msigdbr)  
  require(biomaRt)  
  require(tidyverse)
  require(pheatmap)
  require(matrixStats)
  require(rrvgo)
  require(GOSemSim)
  require(org.Hs.eg.db)
  require(GO.db)
  require(GSEABase)

# ---------------------------
# User inputs: replace / set
# ---------------------------
# expr <- <your matrix>  # rows = genes; columns = samples
# Example: expr <- readRDS("normalized_expression_matrix.rds")
if(!is.null(norm_data) & "E" %in% names(norm_data) & "genes" %in% names(norm_data)){
  print("Pulling expression and gene data from normalized data input.")
  expr=norm_data[["E"]]
  gene_key=norm_data[["genes"]]
} else if(!is.null(norm_data) & (!("E" %in% names(norm_data)) | !("genes" %in% names(norm_data)))){
  print("It appears you have entered a normalized data object, but that internal components are misformatted. See documentation.")
} else if (is.null(norm_data) & (!is.null(expr) & !is.null(gene_key))){
  print("Expression matrix and gene info located.")
} else if(is.null(norm_data) & (is.null(expr) | is.null(gene_key))){
  print("Expression matrix or gene data missing. See documentation.")
}

# ---------------------------
# Basic checks & transformations
# ---------------------------

# Make sure rownames(expr) are gene symbols OR Ensembl IDs (no duplicates per rowname)

# If rows appear to be Ensembl IDs (start with ENSG), map to gene symbols
is_ensembl <- all(grepl("^ENS", rownames(expr)[1:min(30, nrow(expr))]) ) # check formatting

# if(is_ensembl) { # map ensembl ids to hgnc if necessary
#   message("Detected rownames as Ensemb IDs ; attempting mapping to gene symbols via biomaRt.")
#   mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")  # for human; adjust if needed
#   map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                filters = "ensembl_gene_id",
#                values = rownames(expr),
#                mart = mart)
#   
#   
#   # keep only mapped and non-empty symbols; collapse duplicates by choosing first
#   map <- map[map$hgnc_symbol != "", ]
#   symbol_order <- map$hgnc_symbol[match(rownames(expr), map$ensembl_gene_id)]
#   keep_idx <- which(!is.na(symbol_order))
#   expr <- expr[keep_idx, , drop=FALSE]
#   rownames(expr) <- symbol_order[keep_idx]
#   # remove duplicate gene symbols by keeping the row with highest variance
#   dup_sym <- rownames(expr)[duplicated(rownames(expr))]
#   if(length(dup_sym)) {
#     message("Collapsing duplicate gene symbols by max row-variance.")
#     keep_rows <- sapply(unique(dup_sym), function(g) {
#       rows <- which(rownames(expr) == g)
#       rows[which.max(rowVars(expr[rows, , drop=FALSE]))]
#     })
#     nondup_idx <- which(!rownames(expr) %in% dup_sym)
#     keep_idx <- c(nondup_idx, keep_rows)
#     expr <- expr[keep_idx, , drop=FALSE]
#   }
# }
# 
# # Ensure gene names are non-empty and unique
# stopifnot(!any(rownames(expr) == ""))

if(any(duplicated(rownames(expr)))) stop("Duplicate rownames remain after mapping.")

##############################
#    Build gene-set list
##############################

# this currently defines/pulls gene lists internally. This is inefficient if running multiple times, and if/when this is
# added to any package development these can be called from a separate function.

print("Compiling reference genesets")

#--------------------------------------------------------------------------------------------
# 1) Pull broad curated gene sets from MSigDB (example: GO BP or Hallmark)
#    msigdbr provides many collections; by default we pull GO:BP and Hallmarks as they are typically most useful for our biological enrichments.
#--------------------------------------------------------------------------------------------

# Get all GO annotations for all genes
go_bp_sets <- select(org.Hs.eg.db,
                  keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
                  columns = c("GO","ONTOLOGY", "ENSEMBL", "SYMBOL"),
                  keytype = "ENSEMBL")%>%
  filter(ONTOLOGY=="BP")

go_terms<-
  select(GO.db,
         keys = keys(GO.db, keytype = "GOID"),
         columns = c("GOID","ONTOLOGY", "TERM"),
         keytype = "GOID")%>%
  filter(ONTOLOGY=="BP")

go_bp_sets_gene_map<-
  go_terms%>%
  left_join(go_bp_sets, by=c("GOID"="GO", "ONTOLOGY"))%>%
  mutate(gs_name=paste0("GOBP_",str_replace_all(str_to_upper(TERM), "\\ ", "_")))%>%
  group_by(GOID, TERM, gs_name) %>%
  dplyr::summarize(genes = list(unique(SYMBOL)), 
                   ensembl_gene_ids = list(unique(ENSEMBL)), .groups = "drop")

# Load GO Slim definition
# slim_obo <- GSEABase::getOBOCollection("https://current.geneontology.org/ontology/subsets/goslim_generic.obo")
# bp_collection <- GSEABase::GOCollection(go_bp_sets_mapped$GOID)
# 
# gobp_slimmed<- # reduce to GOslim terms
#   GSEABase::goSlim(bp_collection, 
#                    slimCollection = slim_obo, 
#                    ontology =  "BP")
# 
# go_bp_sets_mapped<-
#   filter(go_bp_sets_gene_map, GOID %in% rownames(gobp_slimmed))
go_bp_sets_mapped<-go_bp_sets_gene_map


msig_hallmark <- msigdbr(species = species, collection = "H") %>%
  dplyr::select(gs_name, gene_symbol, ensembl_gene) %>%
  group_by(gs_name) %>%
  dplyr::summarize(genes = list(unique(gene_symbol)), 
                   ensembl_gene_ids = list(unique(ensembl_gene)), .groups = "drop")


# Hounkpe housekeeping (MSigDB name available as HOUNKPE_HOUSEKEEPING_GENES)
# msigdbr can fetch it too; attempt to pull by gs_name filter
hk_msig <- msigdbr(species = "Homo sapiens") %>%
  filter(grepl("HOUNKPE_HOUSEKEEPING_GENES", gs_name, ignore.case=TRUE)) %>%
  pull(ensembl_gene) %>% unique()


# Combine into named list
gs_list <- c(setNames(go_bp_sets_mapped$ensembl_gene_ids, go_bp_sets_mapped$gs_name),
             setNames(msig_hallmark$ensembl_gene_ids, msig_hallmark$gs_name))
gs_list[["HOUNKPE_HOUSEKEEPING_GENES"]] = hk_msig

# 2.Custom defined "technical" gene sets ------------------------------------
#----------------------------------------------------------------------------
#   Add curated "technical" gene sets likely to capture technical variance:
#    - mitochondrial (MT- prefix for symbols)
#    - ribosomal proteins (RPS*, RPL*)
#    - histones (HIST*)
#    - hemoglobins (HBA1, HBB, etc.)
#    - common housekeeping/spike-ins (e.g. ACTB, GAPDH, ERCC if present)
#-----------------------------------------------------------------------------

# Mitochondrial genes by symbol pattern (human genes usually begin MT-)
mito<-gene_key[grepl("^MT-", gene_key$symbol, ignore.case = FALSE), ] 

# Ribosomal proteins

# cytosolic ribosomal proteins.
ribosomal_cytosolic <- gene_key[grepl("^RPS|^RPL", gene_key$symbol),]

#  mitochondrial ribosomal proteins.
ribosomal_mito <- gene_key[grepl("^MRPL|^MRPS", gene_key$symbol),]

# selected GO mitochondrial gene-sets
# here specifies a small number of mechanistic axes that are known to dominate PCs when mitochondrial signal is technical.
go_mito <-go_bp_sets_mapped%>%
  dplyr::filter(gs_name %in% c(
    "GOBP_MITOCHONDRIAL_GENE_EXPRESSION", #mtRNA abundance, RNA integrity, rRNA contamination
    "GOBP_MITOCHONDRIAL_TRANSLATION", # mitoribosome activity;  proxy for mitochondrial RNA load
    "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN", # * note, may be very sensitive to normalization and depth
    "GOBP_OXIDATIVE_PHOSPHORYLATION", # tightly co-expressed nuclear OXPHOS genes; often drivers of main PCs
    "GOBP_MITOCHONDRION_ORGANIZATION" # mitochondrial mass / turnover
  )) %>%
  dplyr::pull(ensembl_gene_ids) %>%
  unlist()%>%
  unique()
  
# Hallmark oxphos
hallmark_oxphos <- msig_hallmark %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  dplyr::pull(ensembl_gene_ids) %>%
  unique()

# Histones (HIST1H, HIST2H, HIST3H etc.)
histones <- gene_key[grepl("^HIST", gene_key$symbol), ]

# Hemoglobins (common symbols)
globins <- 
  filter(gene_key, symbol %in% c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBZ"))

# Housekeeping example
hk <- filter(gene_key, symbol %in% c("ACTB","GAPDH","TUBB","EEF1A1","RPLP0"))

# ERCC genes (sometimes included as spike-in)
ercc <- gene_key[grepl("^ERCC-|^ERCC", gene_key$symbol, ignore.case = TRUE), ] 

# Build the custom list (only include sets with >= min_geneset_size)
custom_sets <- list(
  TECHSIG_MTchr = unique(mito$ensembl_gene_id),
  TECHSIG_GOBP_MITOCHONDRIAL = unique(go_mito),
  
  TECHSIG_RIBOSOMAL_CYTO = unique(ribosomal_cytosolic$ensembl_gene_id),
  TECHSIG_RIBOSOMAL_MITO = unique(ribosomal_mito$ensembl_gene_id),
  
  TECHSIG_HALLMARK_OXPHOS=unique(unlist(hallmark_oxphos)),
  
  TECHSIG_HISTONES = unique(histones$ensembl_gene_id),
  TECHSIG_GLOBINS = unique(globins$ensembl_gene_id),
  TECHSIG_HOUSEKEEPING = unique(hk$ensembl_gene_id),
  TECHSIG_ERCC = unique(ercc$ensembl_gene_id)
)

custom_sets_names<-c(names(custom_sets), "HOUNKPE_HOUSEKEEPING_GENES")# save names for downstream use

# Filter empty sets
custom_sets <- custom_sets[lengths(custom_sets) > 0]

# Combine with msig sets into a named list and filter by membership vs expression genes
all_gs <- c(gs_list, custom_sets)

# Intersect each set to genes present in expr and filter by size
all_gs <- lapply(all_gs, function(g) intersect(g, gene_key$ensembl_gene_id))

all_gs <- all_gs[sapply(all_gs, length) >= min_geneset_size & sapply(all_gs, length) <= max_geneset_size]

message(sprintf("Using %d gene sets after filtering by size and overlap.", length(all_gs)))

##############################
#          Run PCA
##############################

# We'll use svd on scaled data (genes centered) — typical is to center genes (rows).
# We'll run PCA on the transposed matrix so that PCs are sample-level PCs and loadings are per-gene.

# Option: row-center genes (remove gene mean) and optionally scale by sd
mat_centered <- t(scale(t(expr), center = TRUE, scale = TRUE)) # genes standardized across samples

# If many zeros/near-zero sd, scale=TRUE may produce NAs — check:
if(any(is.na(mat_centered))) {
  mat_centered <- t(scale(t(expr), center = TRUE, scale = TRUE))
  warning("Some genes had zero variance; scaled without dividing by sd.")
}
# fit PCA. Note that this is equivalent to having prcomp center and scale.
pca <- prcomp(t(mat_centered), center = FALSE, scale. = FALSE)%>%  # samples in rows
  summary() # add variance estimates

# save loadings
#loadings <- pca$rotation    

# save variance summary
pca_vars<-
  pca$importance%>%
  t()%>%as.data.frame()%>%  # transpose and frame
  janitor::clean_names()%>% # improve column names
  dplyr::rename_with(~ paste0("PC_", .x), everything())%>%
  rownames_to_column("PC")%>%  # preserve PC names
  mutate(PC=factor(PC, levels=gtools::mixedsort(.$PC)))%>%
  arrange(PC)

# If loadings rownames are missing, set them
if(is.null(rownames(pca$rotation))){ rownames(pca$rotation ) <- rownames(mat_centered)}

# Inspect PCA Variances to guide selection of PCs
# Here we find the first PC that explains below 1% of variance, and we display 2X that number.
max_pc<- dplyr::first(which(pca_vars$PC_proportion_of_variance<0.01))

pca_scree<-
  head(pca_vars, 2*max_pc)%>%
  mutate(plot_col=ifelse(PC_proportion_of_variance>=0.01, "tomato", "grey60"))%>%
  {ggplot(data=.,aes(x=PC, y=PC_proportion_of_variance))+
      geom_path(aes(group="1"))+
      geom_point(aes(color=plot_col), size=3)+
      geom_hline(yintercept = 0.01, color="red3", linetype="dashed")+
      scale_y_continuous(limits=c(0, max(.$PC_proportion_of_variance)), 
                                  breaks=seq(0, max(.$PC_proportion_of_variance), 0.01))+
      scale_x_discrete(labels=parse_number(as.character(.$PC)))+
      labs(y="Proportion of Variance Explained.")+
      scale_color_identity()+
      theme_bw()
  }

# Prepare gene-level test statistic-----------------------------------------------------
# --------------------------------------------------------------------------------------
# Prepare gene-level test statistic for PCA correlation
# Compute correlation between each gene (row of mat_centered) and PC scores (pca$x[,k])
# Then Fisher-transform correlations for stability
# ---------------------------------------------------------------------------------------

# number of PCs to test
if(is.null(nPC)){nPC<-max_pc-1}

gene_stats <- sapply(1:nPC, function(k) {
  pc_scores <- pca$x[, k]
  cor_vec <-  # compute Pearson correlation between each gene and the PC
    apply(mat_centered, MARGIN = 1, FUN = function(f){ cor(f, pc_scores, method = "spearman")})
  atanh(cor_vec)}) # Z-Score transformation

colnames(gene_stats) <- paste0("PC", 1:nPC) # preserve appropriate column naming

# 3. Gene Set Analysis ----
##########################################################################################
#          Run Competitive Gene Set Analysis on PC scaled loading correlations
##########################################################################################


# --------------------------- --------------------------- --------------------------- ------
# Run limma::camera, a competetive gene-set enrichment analysis, to associate PCs to GS
# --------------------------- --------------------------- --------------------------- ------

# estimate assumed gene-gene correlation
estimate_intergene_cor <- function(mat){
  cor_sample <- cor(t(mat[ sample(nrow(mat), 2000), ])) # subsample to avoid memory problems
  mean(cor_sample[upper.tri(cor_sample)], na.rm=TRUE)
}

# calculate pairwise gene cors after centering
"Estimating centered pairwise gene correlation."
set.seed(32)
intergene_cor<-cor(t(mat_centered[ sample(nrow(mat_centered), 2500), ]))
inter_cor_est <- mean(intergene_cor[upper.tri(intergene_cor)], na.rm=TRUE)

#convert gene sets to index list
gs_indexed<-
  limma::ids2indices(all_gs, rownames(expr), remove.empty=TRUE)

# Run pre-ranked camera
camera_res_pr<-data.frame()

for(i in 1:nPC){
  # run competetive enrichment test using limma camera
  camera_res_pr<-
    bind_rows(camera_res_pr, 
              limma::cameraPR(statistic = abs(gene_stats[,i]),  
                              index=gs_indexed, 
                              directional=F, # run adirectionally by abs correlation to PC
                              use.ranks=TRUE,
                              inter.gene.cor=inter_cor_est)%>%
                mutate(PC=paste0("PC", i))%>%
                rownames_to_column("gene_set"))
}


# Ensure typing is correct
camera_res_pr <- camera_res_pr %>%
  mutate(
    PC = factor(PC, levels = paste0("PC", seq_len(nPC))),
    gene_set = as.character(gene_set)
  )%>%
  rowwise()%>%
  mutate(gene_set_db=case_when( grepl("TECHSIG", gene_set) ~"TECHSIG", 
                                grepl("GOBP", gene_set) ~"GOBP", 
                               grepl("HALLMARK", gene_set) ~"HALLMARK"))%>%
  left_join(
    dplyr::select(go_bp_sets_mapped, gs_name, GOID), by=c("gene_set"="gs_name")
  )

##########################################
# 4. Summarize and Visualize Results ----
##########################################


# Summarize top enrichments per PC

tech_and_top_per_PC <- 
  filter(camera_res_pr, gene_set %in% custom_sets_names) %>%
  bind_rows(
    filter(camera_res_pr, !(gene_set %in% custom_sets_names))%>%
      group_by(PC)%>%
      slice_min(PValue, n = 5))%>%
  ungroup()%>%
  mutate(gene_set_factor=
           factor(gene_set, 
                  levels = c(custom_sets_names, 
                             unique(.$gene_set[which( !(.$gene_set%in% custom_sets_names))]))))%>%
  arrange(PC, gene_set_factor)


## VISUAL SUMMARY
PC_gene_set_summary_plot<-
  tech_and_top_per_PC%>%
  mutate(plot_col=ifelse(gene_set_factor %in% custom_sets_names, "dodgerblue", "black"))%>%
  {ggplot(data=., aes(x=-log10(FDR), y=gene_set_factor))+
      geom_vline(xintercept = 0)+
      geom_vline(xintercept = -log10(0.05), linetype="dashed", color="red3")+
      geom_point(aes(col=plot_col))+
      scale_color_identity()+
      labs(y="")+
      facet_grid(rows=vars(PC),  space="free_y", scales="free_y")+
      theme_bw()+theme(axis.text.y = element_text(size=5.5), strip.background = element_rect(fill="white"))
  }


# Extract technical enrichment per PC
tech_scores <- camera_res_pr %>%
  filter(gene_set %in% custom_sets_names) %>%
  arrange(PC, FDR)%>%
  group_by(PC) %>%
  mutate(tech_score_group=gene_set[1])%>%
  summarise(tech_score = max(-log10(FDR)), 
            tech_group=ifelse(tech_score> -log10(0.05), 
                              unique(tech_score_group)%>%str_remove("TECHSIG_"), 
                              "No significant technical enrichments."),
            .groups="drop")%>%
  left_join(dplyr::select(pca_vars, PC, PC_proportion_of_variance))

# Panel 2: Technical score
tech_scores_plot <- ggplot(tech_scores, aes(x = PC, y = tech_score)) +
  geom_point(size = 3, aes(color = tech_group)) +
  geom_hline(yintercept = -log10(0.05), 
             color="red3", linetype="dashed") +
  scale_color_manual(values=c("No significant technical enrichments."="grey40",
                              "GOBP_MITOCHONDRIAL"="#E69F00",
                              "MTchr"="red3",
                              "RIBOSOMAL_CYTO"="#56B4E9",
                              "RIBOSOMAL_MITO" = "darkblue",
                              "HALLMARK_OXPHOS" ="#F0E442" ,
                              "HOUNKPE_HOUSEKEEPING_GENES" = "#009E73",
                              "HOUSEKEEPING" = "green3",
                              "HISTONES"= "violet", 
                              "GLOBINS"= "purple", 
                              "ERCC"="tomato"))+

  labs(title = "Technical Signal per PC", color="Technical Gene-Set\nw/ Max Enrichment",
       y = "Maximum enrichment among custom technical sets.\n[ –log10(FDR) ]",
       x = "") +
  theme_bw()


tech_scores_plot_2 <- 
  ggplot(tech_scores, aes(x = PC_proportion_of_variance, y = tech_score)) +
  geom_hline(yintercept = -log10(0.05), 
             color="red3", linetype="dashed") +
  geom_point(size = 6, aes(color = tech_group))+
  geom_text(aes(label=str_remove(PC, "PC")), color="white", size=3)+
  scale_color_manual(values=c("No significant technical enrichments."="grey40",
                              "GOBP_MITOCHONDRIAL"="#E69F00",
                              "MTchr"="red3",
                              "RIBOSOMAL_CYTO"="#56B4E9",
                              "RIBOSOMAL_MITO" = "darkblue",
                              "HALLMARK_OXPHOS" ="#F0E442" ,
                              "HOUNKPE_HOUSEKEEPING_GENES" = "#009E73",
                              "HOUSEKEEPING" = "green3",
                              "HISTONES"= "violet", 
                              "GLOBINS"= "purple", 
                              "ERCC"="tomato"))+
  
  labs(title = "Technical Signal per PC", color="Technical Gene-Set\nw/ Max Enrichment",
       y = "Maximum enrichment among custom technical sets.\n[ –log10(FDR) ]",
       x = "Prop. Variance per PC") +
  theme_bw()+theme(aspect.ratio = 1)

# ---------------------------
# Save/export results
# ---------------------------

return(
  list("PC_enrichments_top"=tech_and_top_per_PC,
       "PC_enrichments_full" = camera_res_pr, 
       "plot_PC_enrichments"= PC_gene_set_summary_plot,
       "plot_PC_techscores"=tech_scores_plot,
       "plot_PC_techscores_scatter"=tech_scores_plot_2,
       "plot_scree" = pca_scree, 
       "PCA_vals"=pca$x,
       "PCA_loadings"=pca$rotation,
       "PCA_var"= pca_vars,
       "nPC"=nPC
       )
)


}





