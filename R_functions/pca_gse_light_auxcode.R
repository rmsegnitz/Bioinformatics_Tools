# auxillary code removed from pca_gse_light.R

##########################################################
#     Reduce gene-sets to cluster representatives.
##########################################################

# camera_res_pr is your results data.frame from cameraPR
# It must have columns: gene_set (name/id), PValue or FDR

# prefer FDR, fallback to PValue

camera_res_pr<-camera_res_pr%>%
  mutate(score_val=ifelse(!is.na(FDR),
                          -log10(FDR + 1e-300),
                          -log10(PValue + 1e-300)))

# Reduce to top scoring genesets per PC, this should improve efficiency downstream.
top_GO <- 
  camera_res_pr %>% 
  filter(!is.na(GOID)) %>%
  #summarize(score = max(-log10(FDR + 1e-300), -log10(PValue + 1e-300))) %>% 
  arrange(PC, PValue) %>% 
  group_by(PC) %>% 
  slice_head(n = topN) %>% 
  pull(GOID)%>%
  unique()

# reduce terms
camera_res_pr <- camera_res_pr %>%
  filter(GOID %in% top_GO | gene_set_db %in% c("HALLMARK", "TECHSIG"))


# Build scores named vector (names = GO IDs)
scores_vec <- setNames(camera_res_pr$score_val, camera_res_pr$GOID)
# optionally restrict to GO IDs with non-zero score:
scores_vec <- scores_vec[!is.na(names(scores_vec)) & !is.na(scores_vec)]


# semData for human BP
#semData <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = TRUE)

# calculate similarity matrix for the GO IDs of interest (names of scores_vec)
simMatrix <- rrvgo::calculateSimMatrix(names(scores_vec),
                                       orgdb = "org.Hs.eg.db",
                                       ont = "BP",
                                       method = "Rel")  # Rel or Wang; Rel works well with scores

simMatrix_reduced <- rrvgo::reduceSimMatrix(simMatrix,
                                            scores = scores_vec,
                                            threshold = 0.7,
                                            orgdb = "org.Hs.eg.db")%>%
  group_by(parent)%>%
  add_tally(name = "GOID_collapsed")

simMatrix_reduced_alt <- rrvgo::reduceSimMatrix(simMatrix,
                                                scores = "uniqueness",
                                                threshold = 0.7,
                                                orgdb = "org.Hs.eg.db")%>%
  group_by(parent)%>%
  add_tally(name = "GOID_collapsed")

# representatives
rep_GOID <- unique(simMatrix_reduced$parent)  # these are GO IDs chosen as representatives


# filter results to reduced gene-sets and recalculate FDR
print("Reducing to non-redundant gene-sets.")
camera_res_pr<-
  filter(camera_res_pr, gene_set_db !="GOBP" | GOID  %in% rep_GOID)#%>%
# group_by(PC)%>%
# mutate(FDR=p.adjust(PValue, method="BH"))
