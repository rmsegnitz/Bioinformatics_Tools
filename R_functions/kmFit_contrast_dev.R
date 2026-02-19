#' Run pairwise model comparisons with emmeans
#'
#' @param fit model fit from lm( ) or lmer( )
#' @param contrast_var Character vector of variable in model to run contrasts of
#' #' @param contrast_spec 2-Column dataframe specifying contrasts of variable in model to run contrasts of.
#' @param to_model_gene Formatted data from kimma_cleaning( ), subset to gene of interest
#' @param genotype_name Character string. Used internally for kmFit_eQTL
#'
#' @return data frame with contrast model results
#' @keywords internal

kmFit_contrast_dev <- function(fit, contrast_var, contrast_spec, to_model_gene, genotype_name){
  contrast.i <- term <- p.value <- contrast_ref <- contrast_lvl <- contrast <- null.value <- estimate <- NULL
  contrast.result <- data.frame()

  #Fit variables in contrast_var
  for(contrast.i in contrast_var){
    contrast.result.temp <- NULL
    #Class
    i.split <- strsplit(contrast.i, split=":")[[1]]
    contrast.is.numeric <- unlist(lapply(to_model_gene[,i.split], is.numeric))

    #If any numeric
    if(any(contrast.is.numeric)){
      contrast.result.temp <- tryCatch({
        emmeans::emtrends(fit, adjust="none",var=i.split[contrast.is.numeric],
                          stats::as.formula(paste("pairwise~", i.split[!contrast.is.numeric],
                                                  sep="")))$contrasts %>%
          broom::tidy() %>%
          dplyr::mutate(term = gsub(":","*", contrast.i)) %>%
          tidyr::separate(contrast, into=c("contrast_ref","contrast_lvl"),
                          sep=" - ")
      }, error=function(e){ return(NULL) })

      #if model ran, add to results
      if(is.data.frame(contrast.result.temp)){

        contrast.result <- contrast.result.temp %>%
          dplyr::bind_rows(contrast.result)
      }
    } else if(is.null(contrast_spec)){
      contrast.result.temp <- tryCatch({
        emmeans::emmeans(fit, adjust="none",
                         stats::as.formula(paste("pairwise~", contrast.i, sep="")))$contrasts %>%
          broom::tidy() %>%
          dplyr::mutate(term = gsub(":","*", contrast.i)) %>%
          tidyr::separate(contrast, into=c("contrast_ref","contrast_lvl"),
                          sep=" - ")


      }, error=function(e){ return(NULL) })} else if(!is.null(contrast_spec)){

        if(ncol(contrast_spec)!=2){stop("Improper dimensions; Check formatting of specified contrasts.")}
        if(!("contrast_ref"%in% colnames(contrast_spec)) | !("contrast_lvl"%in% colnames(contrast_spec))){
          stop("Improper contrast formatting; Check contrast_spec dataframe.")}

        if(grepl(":", contrast.i)){

          contrast.i.split<-stringr::str_split(contrast.i, ":")
           contrast.i.a<-unlist(contrast.i.split)[1]
            contrast.i.b<-unlist(contrast.i.split)[2]

          contrast_lvls<-model.frame(fit)[,c(contrast.i.a, contrast.i.b)]%>%interaction(sep = " ")%>%levels()
            default_vec<-rep(0, length(contrast_lvls))
              names(default_vec)<-contrast_lvls
        } else {
          contrast_lvls<-model.frame(fit)%>%pull(contrast.i)%>%levels()
          default_vec<-rep(0, length(contrast_lvls))
          names(default_vec)<-contrast_lvls
        }

        if(all(contrast_spec$contrast_ref %in% contrast_lvls) &
           all(contrast_spec$contrast_lvl %in% contrast_lvls)){

          contrasts_list<-list()

              for(i in 1:nrow(contrast_spec)){
                # specify contrast vector
                contrast_vec<-default_vec
                contrast_vec[which(names(contrast_vec)==contrast_spec$contrast_ref[i])]<- -1
                contrast_vec[which(names(contrast_vec)==contrast_spec$contrast_lvl[i])]<- 1
                # save to list
                contrasts_list[[paste(contrast_spec[i, "contrast_lvl"], contrast_spec[i, "contrast_ref"], sep=" - ")]]<-
                  contrast_vec
              }


          } else {stop("Check contrasts formatting. Contrast levels do not correspond to those of variable specified.")}




        contrast.result.temp <- tryCatch({
          emmeans::emmeans(fit, adjust="none", spec=as.formula(paste("~", contrast.i.a, "*", contrast.i.b)))%>%
            emmeans::contrast(contrasts_list, adjust="none")%>%
            broom::tidy() %>%
            dplyr::mutate(term = gsub(":","*", contrast.i)) %>%
            tidyr::separate(contrast, into=c("contrast_lvl","contrast_ref"),
                            sep=" - ")
      }, error=function(e){ return(NULL) })}

      #if model ran, add to results
      if(is.data.frame(contrast.result.temp)){
        #fix genotype names
        if(!is.null(genotype_name)){
          if(grepl(genotype_name, contrast.i)){
          contrast.result.temp <- contrast.result.temp %>%
            dplyr::mutate(dplyr::across(c(contrast_ref, contrast_lvl),
                                 ~gsub(genotype_name, paste0(genotype_name,"_"), .)))
        }}

        contrast.result <- contrast.result.temp %>%
          dplyr::bind_rows(contrast.result)
      }
}


  if(is.null(contrast_spec)){
  contrast.result.format <- contrast.result %>%
    dplyr::rename(variable=term, pval=p.value) %>%
    dplyr::select(-null.value) %>%
    #Switch estimate sign to match lvl minus ref calculation
    #THE REF AND LVL VALUES ARE INCORRECT UNTIL YOU DO THIS
    dplyr::mutate(estimate = -estimate)%>%
    # add flag for interpretation
    dplyr::mutate(higher_in=ifelse(estimate>0, contrast_lvl, contrast_ref))%>%
    dplyr::relocate(higher_in, .before = estimate)

  } else if(!is.null(contrast_spec)){
    contrast.result.format <- contrast.result %>%
      dplyr::rename(variable=term, pval=p.value) %>%
      dplyr::select(-null.value)%>%
      dplyr::mutate(higher_in=ifelse(estimate>0, contrast_lvl, contrast_ref))%>%
      dplyr::relocate(higher_in, .before = estimate)

  }


  return(contrast.result.format)
}

