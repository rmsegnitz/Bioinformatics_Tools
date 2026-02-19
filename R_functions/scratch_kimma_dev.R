 

load("~/Documents/WORK/UW/ALTMAN/GITHUB/GRADUATE_batch2/data_clean/P662_2_voom_meta.RData")
temp_data<-dat_meta
#rownames(temp_data$targets)<-temp_data$targets$libID

tempE<-
  subset_voom(temp_data, gene_keep = rownames(temp_data$E)[1:100])

tempE_wex<-
  tempE$targets%>%
  dplyr::left_join(
    tempE$E[1:100, ]%>%t()%>%as.data.frame()%>%tibble::rownames_to_column("libID"))


lm_eg<-lm(ENSG00000000457 ~ visit, data=tempE_wex)

lmer_erg<-lme4::lmer(ENSG00000000457 ~ visit + (1|ptID), data=tempE_wex)

lmer_erg_ixn<-lme4::lmer(ENSG00000000457 ~ treatment*visit + (1|ptID), data=tempE_wex)


# fit equivalent kmFits
km_lm_eg<-kmFit(dat=tempE, 
                model =  "~ visit", weights = F,
                run_lm = T, patientID = "ptID", libraryID = "libID")

km_lm_eg_contrast<-kmFit(dat=tempE, 
                model =  "~ visit", weights = F, run_contrast = T, contrast_var = "visit",
                run_lm = T, patientID = "ptID", libraryID = "libID")

lm_eg_contrast<-pairs(emmeans::emmeans(lm_eg, specs = "visit"), adjust="none")


lm_eg_contrast_formtted<-
emmeans::emmeans(lm_eg, adjust="none",
                 stats::as.formula(paste("pairwise~", "visit", sep="")))$contrasts%>%
  broom::tidy()%>%
  dplyr::mutate(term = gsub(":","*", "visit")) %>%
  tidyr::separate(contrast, into=c("contrast_ref","contrast_lvl"),
                  sep=" - ")



View()
# test FDR correction
test_FDR<-km_lm_eg_contrast$lm.contrast%>%
  group_by(variable, contrast_ref, contrast_lvl)%>%
  mutate(FDR_by_contrast=p.adjust(pval, method = "BH"))%>%
  group_by(variable)%>%
  mutate(FDR_across_contrast=p.adjust(pval, method = "BH"))

# Work out how to specify contrast matrix. 
# default is pairwise
contrast_var<-"visit"
contrast_lvls<-levels(dplyr::pull(tempE_wex, get(contrast_var)))

dummy_pairs_list<-
  gtools::combinations(n=length(contrast_lvls), r = 2, 
                       levels(dplyr::pull(tempE_wex, get(contrast_var))))%>%
  as.data.frame()%>%
  arrange(match(V1, contrast_lvls), match(V2, contrast_lvls))%>%
  dplyr::rename(contrast_ref=V1, contrast_lvl=V2)%>%
  mutate(contrast_val=1, 
         contrast_ref_code=ifelse(contrast_val==1, -1, 0), 
         contrast_lvl_code=ifelse(contrast_val==1, 1, 0))
                       
# Coerce dataframe into list of contrasts.
contrasts_ap<-
  dummy_pairs_list[1:5,]

default_vec<-rep(0, length(contrast_lvls))
names(default_vec)<-contrast_lvls
  
contrasts_list<-list()

for(i in 1:nrow(contrasts_ap)){
  # specify contrast vector
  contrast_vec<-default_vec
  contrast_vec[which(names(contrast_vec)==contrasts_ap$contrast_ref[i])]<- -1
  contrast_vec[which(names(contrast_vec)==contrasts_ap$contrast_lvl[i])]<- 1
  # save to list
  contrasts_list[[paste(contrasts_ap[i, "contrast_lvl"], contrasts_ap[i, "contrast_ref"], sep="_vs_")]]<-
    contrast_vec
}

  
# Apply contrast
start<-Sys.time()
lm_eg_contrast_ap<-
  emmeans::contrast(emmeans::emmeans(lm_eg, specs = "visit"), contrasts_list, adjust="none")%>%
  as.data.frame()%>%
  tidyr::separate(contrast, c("contrast_lvl", "contrast_ref"), sep = "_vs_")
Sys.time()-start

# compare to pairs
start<-Sys.time()
test_pairs<- 
  pairs(emmeans::emmeans(lm_eg, specs = "visit"), adjust="none")%>%
  as.data.frame()%>%
  tidyr::separate(contrast, c("contrast_lvl", "contrast_ref"), sep = " - ")
Sys.time()-start

View(as.data.frame(lm_eg_contrast_ap))
View(as.data.frame(lm_eg_contrast))


# Work out package code to handle contrast specification

test_pre<-
  kmFit(dat=tempE, 
        model =  "~ treatment*visit + (1|ptID)", 
        use_weights = F, run_contrast = T, 
        contrast_var =  c("visit"),
        run_lm = T, patientID = "ptID", libraryID = "libID")

test_contrasts<-dummy_pairs_list[1:10, ]%>%
  dplyr::select(contrast_ref, contrast_lvl)

test_dev<-
  kmFit_dev(dat=tempE, 
        model =  "~ treatment*visit + (1|ptID)", 
        use_weights = F, run_contrast = T, 
        contrast_var =  c("visit"),
        contrast_spec=test_contrasts,
        run_lm = T, patientID = "ptID", libraryID = "libID")


ixn_lvls<-tempE$targets[,c("treatment", "visit")]%>%interaction(sep = " ")%>%levels()

test_contrasts_ixn<-
  gtools::combinations(n=length(ixn_lvls), r = 2, ixn_lvls)%>%
  as.data.frame()%>%
  dplyr::rename(contrast_ref=V1, contrast_lvl=V2)%>%
  separate(contrast_ref, c("treatment_ref", "visit_ref"), sep=" ", remove = F)%>%
  separate(contrast_lvl, c("treatment_lvl", "visit_lvl"), sep=" ", remove = F)%>%
  mutate(treatment_ref=factor(treatment_ref, levels(tempE$targets$treatment)), 
         visit_ref =factor(visit_ref, levels=levels(tempE$targets$visit)), 
         treatment_lvl=factor(treatment_lvl, levels=levels(tempE$targets$treatment)), 
         visit_lvl=factor(visit_lvl, levels=levels(tempE$targets$visit)))%>%
  mutate(contrast_type=case_when(treatment_ref==treatment_lvl ~ "within-treatment", 
                                 visit_ref==visit_lvl ~ "cross-sectional"))%>%
  filter(!is.na(contrast_type))%>%
  mutate(swap_trt=ifelse(contrast_type=="cross-sectional" & treatment_lvl=="Placebo", T,F))%>%
  mutate(swap_visit=ifelse(contrast_type=="within-treatment" & as.numeric(visit_lvl)<=as.numeric(visit_ref), T,F))%>%
  gggenomes::swap_if(swap_trt | swap_visit, contrast_ref, contrast_lvl)%>%
  separate(contrast_ref, c("treatment_ref", "visit_ref"), sep=" ", remove = F)%>%
  separate(contrast_lvl, c("treatment_lvl", "visit_lvl"), sep=" ", remove = F)%>%
  mutate(NAC_contrast= case_when(contrast_type=="cross-sectional" & visit_ref %in% c("-1", "S4", "S9", "S14") ~ T,
                                 contrast_type=="within-treatment" & 
                                 visit_ref %in% c("-1", "S4", "S9", "S14") & visit_lvl %in% c("-1", "S4", "S9", "S14") ~ T, 
                                 .default = F))

# Ok so we have reduced contrast combinations from 
(length(ixn_lvls)*(length(ixn_lvls)-1))/2
# to 
nrow(test_contrasts_ixn)
# Sensible contrasts.
# and

# NAC contrasts that we are actually interested in.
test_contrasts_ixn_NAC<-
  filter(test_contrasts_ixn, NAC_contrast)
nrow(test_contrasts_ixn_NAC)


NAC_contrasts_for_km<-
  test_contrasts_ixn_NAC%>%
  dplyr::select(contrast_ref, contrast_lvl)

# Test contrast fitting with NAC contrasts

test_dev_NAC<-
  kmFit_dev(
  dat=temp_data, 
  model =  "~ treatment*visit + (1|ptID)" ,
  use_weights = F,
  run_contrast = T ,
  contrast_var =  c("treatment:visit"),
  contrast_spec = NAC_contrasts_for_km,
  run_lme = T,
  patientID = "ptID",
  libraryID = "libID")

test_og_NAC<-
  kmFit(
    dat=temp_data, 
    model =  "~ treatment*visit + (1|ptID)" ,
    use_weights = F,
    run_contrast = T ,
    contrast_var =  c("treatment:visit"),
    #contrast_spec = NAC_contrasts_for_km,
    run_lme = T,
    patientID = "ptID",
    libraryID = "libID")

  