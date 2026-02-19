# Format specified contrasts lists of different lengths
contrast_lists_variable<-list()
for( i in c(10, 20, 30, 40, 50)){
  
  contrasts_ap<-
    dummy_pairs_list[1:i,]
  
  contrasts_list<-list()
  
  for(j in 1:i){
    # specify contrast vector
    contrast_vec<-default_vec
    contrast_vec[which(names(contrast_vec)==contrasts_ap$contrast_ref[j])]<- -1
    contrast_vec[which(names(contrast_vec)==contrasts_ap$contrast_lvl[j])]<- 1
    # save to list
    contrasts_list[[paste(contrasts_ap[j, "contrast_lvl"], contrasts_ap[j, "contrast_ref"], sep="_vs_")]]<-
      contrast_vec
    
  }
  contrast_lists_variable[[paste0(i, "_contrasts")]]<-contrasts_list
}




# run some timing of methods to compare ###

run_times<-
  data.frame("run"=1:3000, 
             "method"=c(rep("pairwise_66_contrast", 500), 
                        rep("specified_10_contrasts", 500), 
                        rep("specified_20_contrasts", 500), 
                        rep("specified_30_contrasts", 500), 
                        rep("specified_40_contrasts", 500), 
                        rep("specified_50_contrasts", 500)), runtime=NA)

for(i in run_times$run){
  if(run_times$method[i]=="pairwise_66_contrast"){
    start<-Sys.time()
    test_pairs<- 
      pairs(emmeans::emmeans(lm_eg, specs = "visit"), adjust="none")%>%
      as.data.frame()%>%
      tidyr::separate(contrast, c("contrast_lvl", "contrast_ref"), sep = " - ")
    run_times$runtime[i]<-Sys.time()-start
  } else if(grepl("specified", run_times$method[i])){
    start<-Sys.time()
    lm_eg_contrast_ap<-
      emmeans::contrast(emmeans::emmeans(lm_eg, specs = "visit"), 
                        contrast_lists_variable[[str_remove(run_times$method[i], "specified_")]], adjust="none")%>%
      as.data.frame()%>%
      tidyr::separate(contrast, c("contrast_lvl", "contrast_ref"), sep = "_vs_")
    run_times$runtime[i]<-Sys.time()-start
    
  }
}

runtime_plot<-
run_times%>%
mutate(contrast_method=ifelse(grepl("pairwise", method), "All pairwise.", "Specified Contrasts."))%>%
  mutate(contrast_number=parse_number(method))%>%
ggplot( aes(x=contrast_number, y=runtime))+
  ggbeeswarm::geom_quasirandom(aes(color=contrast_method), alpha=0.5)+
  scale_x_continuous(breaks=c(10,20,30,40,50,66))+
  theme_bw()+
  labs(y= "Runtime for single gene.", x="Number of contrasts computed.", color="Contrast Method")


runtime_plot_scaled_up<-
  run_times%>%
  mutate(contrast_method=ifelse(grepl("pairwise", method), "All pairwise.", "Specified Contrasts."))%>%
  mutate(contrast_number=parse_number(method))%>%
  ggplot( aes(x=contrast_number, y=runtime*15000))+
  ggbeeswarm::geom_quasirandom(aes(color=contrast_method), alpha=0.5)+
  scale_x_continuous(breaks=c(10,20,30,40,50,66))+
  scale_y_continuous(breaks=seq(40, 400, 20))+
  theme_bw()+
  labs(y= "Approx runtime for 15,000 genes.", x="Number of contrasts computed.", color="Contrast Method")

