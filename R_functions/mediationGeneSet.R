
# Mediation Analysis for Gene Sets 
# Author: Max Segnitz, msegnitz@uw.edu
# Started April 2020
#
# © Richard M Segnitz 2020
# License: This software is licensed under GNU General Public License, and may
# be modified and/or redistributed under GNU GPL version 3 or later. License details
# can be found in the accompanying this script, or at  (http://www.gnu.org/licenses/).
#
# DESCRIPTION:
# Contains function to run mediation analysis over a set of genes or gene modules
# (or any set of mediators, for that matter). 
#
# mediationGeneSet() runs automated model fitting and mediation analysis over a set 
#         of genes with a given outcome, model structure, and specified contrasts.
#
# The core function takes a number of inputs required to specify the desired input 
# models for mediation analysis, run mediation analysis over the desired sets of 
# contrasts, and automate these processes over a set of specified genes or modules.

#'##############################
### 1) mediationGeneSet()   ####
#'##############################

# REQUIRED
#  model.data = (data.frame) Data frame with design metadata and gene expression
#  gene.list = (character) Vector of genes (or modules) to run mediation analysis on.
#  outcome = (character string) Outcome that we are looking for mediation effect on. 
#  treatment = (character string) Treatment variable used in the models (which has direct effect on outcome)
#  t.c.contrasts = (list) A list of vectors specifying the desired "control" vs "treatment" contrasts. Takes the form: list(c("treatment", "control"), c("treatment", "control"))


# OPTIONAL
#  covariates = (character) Vector specifying any additional covariates to be included in the models.
#  random = (logical) Whether or not to fit a random effecs model. If TRUE, function used "lme4" and model formula must be specified as such. 
#  random.effect = (character string) Factor used as RE in the models.  Currently only supports single, random intercept RE.
#  interaction = (logical) Whether or not interaction between mediator and treatment is to be specified. 
#  out.dir = (character string) Filepath to directory in which to save outputs.
#  save.output = (logical) Whether to write analysis output to file.
#  plots = (logical) Indicates whether or not to produce coefficient plots. 
#  save.plots = (logical) Whether to write plots to file. Must also have plots=TRUE to function.
#  plot.dir = (character string) Filepath to directory where figures are to be saved (will create subdirectory witin this directory) 
#  plot.height = (numeric) Height in inches of saved plot (deafaults to 4)
#  plot.width = (numeric) Height in inches of saved plot (deafaults to 6)
#  color.groups = (character) A named vector of plotting colors. Should include all the levels of the contrasts requested. Defaults to NULL, in which case
#                  ("#000000","#E69F00","#85C0F9","#601A4A") is the default pallette. If you wish to override the color specification for "Average" or "Total"
#                  they can be included in the vector of specified colors. 
#  boot = (logical) Whether to generate CIs by non-parametric bootrapping. default is FALSE, which will produce CIs via quasi Bayesian approximation.
#  boot.ci.type = (character) Which type of CIs to generate if boot=TRUE. Defaults to "perc" (percentile-based), but can be "bca"(bias-coorected & accellerated)
#  sims = (numeric) Number of simulations to use when generating confidence intervals. Defaults to 1000. 
#  robustSE = (logical) Whether or not to use robust, heteroscedasticity consistent errors when usuing quasi-Bayesian estimation. Default = FALSE.
#  seed = (numeric) Set.seed prior to running mediate() and random bootstraps in order to reproduce results. Default = NULL

# EXAMPLE USAGE
# mediation.analysis<-
#   mediationGeneSet(model.data = testDesignSubwEx,
#                    gene.list = gene.subset,
#                    outcome = "TNSSMAX",
#                    treatment = "group",
#                    covariates = c("visit"),
#                    interaction = TRUE,
#                    t.c.contrasts = list(c("AMG 157/SCIT" , "Placebo/Placebo"),
#                                         c("AMG 157/SCIT", "Placebo/SCIT" ),
#                                         c("Placebo/SCIT", "Placebo/Placebo")),
#                    save.output = TRUE,
#                    out.dir = res.dir,
#                    random = TRUE ,
#                    random.effect = "PID",
#                    plots=TRUE,
#                    color.groups = c( "Placebo/Placebo" = "#2B2725","Placebo/SCIT" = "#EE9A00", 
#                                       "AMG 157/Placebo" ="#3A5FCD", "AMG 157/SCIT" ="#FF0460"),
#                    plot.dir = fig.dir)



########### DEFINE INPUTS ###############
mediationGeneSet<- function(model.data,
                            gene.list,
                            outcome,
                            treatment,
                            covariates=NULL,
                            interaction = FALSE,
                            t.c.contrasts,
                            save.output=FALSE,
                            out.dir,
                            plots=FALSE,
                            color.groups = NULL,
                            plot.height = 4,
                            plot.width = 6,
                            save.plots = FALSE,
                            plot.dir,
                            random=FALSE,
                            random.effect,
                            boot=FALSE,
                            sims=1000,
                            seed=NULL,
                            boot.ci.type="perc",
                            robustSE=FALSE
                            ){

  ########## LOAD PACKAGES ############# 
  set.seed(2828)
  
  # Required custom dependencies
  # Extract mediation summary
  source("https://raw.githubusercontent.com/rmsegnitz/Bioinformatics_Tools/master/R_functions/extractMediationSummary.R")
  

  # Data manipulation and figures
  require(tidyverse, quietly = TRUE)
  require(stringi, quietly = TRUE)
  require(stringr, quietly = TRUE)
  # mediation analysis
  require(mediation, quietly = TRUE)
  # Progress tracking
  require(svMisc, quietly = TRUE)

  # Plotting
  if(plots){require(cowplot, quietly = TRUE)
            require(DiagrammeR, quietly = TRUE)
            require(DiagrammeRsvg, quietly = TRUE)
            require(rsvg, quietly = TRUE)}
  
  
  ###### SETUP OUTPUT STORAGE ######
  mediation.anovas<-list()
  mediation.models<-list()
  mediation.model.sums<-list()
  mediation.model.singularity<-list()
  mediation.output<-list()
  mediation.estimate.simulations<-list()
  mediation.summary<-list()
  mediation.summary.mat<-list()
  mediation.plots<-list()

 ##### SETUP PLOTTING AES #####

if(!is.null(color.groups) & !all(unlist(t.c.contrasts) %in% names(color.groups))){print("WARNING: not all contrast levels present in plot colors provided.")}

if(is.null(color.groups)){plot.colors = c("#000000","#E69F00","#85C0F9","#601A4A")}
  
if(!is.null(color.groups)){if("Average" %in% names(color.groups) & "Total" %in% names(color.groups)){plot.colors<-color.groups} else
  if("Average" %in% names(color.groups)){plot.colors<- c(color.groups, c("Total" = "#601A4A"))} else if("Total" %in% names(color.groups)){
    plot.colors<- c(color.groups, c("Average" = "#85C0F9"))} else {plot.colors<- c(color.groups, c("Average" = "#85C0F9" , "Total" = "#601A4A"))}}


  
###########  RUN MEDIATION IF RANDOM = TRUE ####################
  if(random){
    # load lme4 for random effects model fitting
    require(lme4, quietly = TRUE)
    
    # for each gene in list, construct formulas for component models.
    for(g in 1:length(gene.list)) { 
      progress(g, length(gene.list))
      i<-gene.list[g]
      # Specify formula for mediator model
      if(is.null(covariates)){ med.formula<-as.formula(
        paste(i, "~",treatment, "+", paste("(1|",random.effect , ")", sep=""), sep=" "))
        } else {med.formula<-as.formula(
        paste(i, "~",treatment,"+",paste(covariates, collapse=" + "), "+", paste("(1|",random.effect , ")", sep=""), sep=" "))}
      
      # specify interaction if indicated.
      if(interaction){

      # specify formula for outcome model
      if(is.null(covariates)){      out.formula<-as.formula(
        paste(paste(outcome, "~", sep=" "),
              paste(i, treatment, sep="*"),"+" , 
              paste("(1|",random.effect , ")", sep=""), sep=" "))}else {
      out.formula<-as.formula(
        paste(paste(outcome, "~", sep=" "),
              paste(i, treatment, sep="*"),"+", 
              paste(covariates, collapse=" + "),"+" , 
              paste("(1|",random.effect , ")", sep=""), sep=" "))}} else{ # else no interaction if not specified.
                if(is.null(covariates)){out.formula<-as.formula(
                  paste(paste(outcome, "~", sep=" "),
                        i,"+", treatment,"+" , 
                        paste("(1|",random.effect , ")", sep=""), sep=" "))}else{  
                  out.formula<-as.formula(
                  paste(paste(outcome, "~", sep=" "),
                        i,"+", treatment, "+",
                        paste(covariates, collapse=" + "),"+" , 
                        paste("(1|",random.effect , ")", sep=""), sep=" "))}}
      
      # Fit mediator model
      med.fit<-lmer( med.formula, model.data)
      
      # Fit outcome (dependent variable) model
      out.fit<-lmer(out.formula, model.data)
      
      # Save input model output
      mediation.model.sums[[i]][[paste(i,"mediator",sep="_")]]<-summary(med.fit)
      mediation.model.sums[[i]][[paste(i,"outcome",sep="_")]]<-summary(out.fit)
      
      mediation.models[[i]][[paste(i,"mediator",sep="_")]]<-broom.mixed::tidy(med.fit)
      mediation.models[[i]][[paste(i,"outcome",sep="_")]]<-broom.mixed::tidy(out.fit)
      
      mediation.anovas[[i]][[paste(i,"mediator_Anova",sep="_")]]<-broom::tidy(car::Anova(med.fit))
      mediation.anovas[[i]][[paste(i,"outcome_Anova",sep="_")]]<-broom::tidy(car::Anova(out.fit))
      
      mediation.model.singularity[[i]][[paste(i,"mediator",sep="_")]]<-isSingular(med.fit)
      mediation.model.singularity[[i]][[paste(i,"outcome",sep="_")]]<-isSingular(out.fit)
      
      # Run moderation analysis for each set of contrasts
      for (j in 1:length(t.c.contrasts)){
        
        print(paste("Running Mediation Analysis |", i,"|", "Treatment:",t.c.contrasts[[j]][1],", Control:",t.c.contrasts[[j]][2], sep=" "))
        
        # run mediation analysis
        set.seed(seed)
        results = mediate(med.fit, out.fit, treat=treatment, mediator=i,
                          control.value = t.c.contrasts[[j]][2], 
                          treat.value = t.c.contrasts[[j]][1], boot = boot, sims = sims, boot.ci.type = boot.ci.type, robustSE = robustSE)
        
        #######  PRODUCE MEDIATION PLOTS ###############
        
        if(plots) { 
          # temporary internal functions for massaging the output
          temp.fun<- function(x){ifelse(stringi::stri_detect_fixed(x, "(control)"), paste(t.c.contrasts[[j]][2]), 
                                        ifelse(stringi::stri_detect_fixed(x, "(treated)"), paste(t.c.contrasts[[j]][1]), 
                                               ifelse(stringi::stri_detect_fixed(x, "(average)"), "Average", "Total")))
          }
          
          temp.fun2<- function(x){ifelse(stringi::stri_detect_fixed(x, "Prop."), "Proportion", "Effect")}
          
          ## FOREST PLOTS ##
          
          # create plotting data frame from mediation output 
          plot_df<-as.data.frame(extractMediationSummary(results))%>%
            rownames_to_column("Effect")%>%
            mutate(effect.type = temp.fun2(Effect))%>%
            mutate(group = temp.fun(Effect))%>%
            separate(Effect, sep="\\(", into=c("Effect", "temp.label"))%>%
            mutate(group=fct_relevel(group, paste(t.c.contrasts[[j]][1]), paste(t.c.contrasts[[j]][2]), "Average", "Total"))%>%
            mutate(Effect=fct_relevel(Effect, "ADE", "ACME", "Total Effect"))%>%
            suppressWarnings() # suppresses unimportant warning from separate command
          
          # Define x axis range (use extreme values for all coef/proportion columns so that axes match)
          x.axis.lims<-plot_df%>%
            select_if(is.numeric)%>%
            dplyr::select(-`p-value`)%>%
            gather()%>%
            dplyr::select(value)%>%
            range(na.rm = T)
          
          # Plot the mediation coefficient comparisons
          plot1<- plot_df%>%
            filter(effect.type=="Effect")%>%
            ggplot(aes(color=group, fill=group))+
            geom_vline(xintercept = 0, color=gray(1/2), lty=2)+
            geom_linerange(aes(y=Effect, xmin=`95% CI Lower`, xmax = `95% CI Upper`), lwd=1, position= position_dodge(width=1/2))+
            geom_pointrange(aes(y = Effect, x= Estimate, xmin=`95% CI Lower`, xmax = `95% CI Upper`), 
                            lwd = 1/2, position=position_dodge(width=1/2), shape=21)+
            scale_color_manual(values = plot.colors )+
            scale_fill_manual(values = plot.colors )+
            ylab("")+
            xlim(x.axis.lims)+
            ggtitle(paste(outcome, "Mediation by", i, sep= " "))+
            theme_bw()+
            theme(plot.title = element_text(hjust = 0.5),
                  legend.title = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
          
          # plot proportion of effect explained by mediation.
          plot2<- plot_df%>%
            filter(effect.type=="Proportion")%>%
            ggplot(aes(color=group, fill=group))+
            geom_vline(xintercept = 0, color=gray(1/2), lty=2)+
            geom_linerange(aes(y=Effect, xmin=`95% CI Lower`, xmax = `95% CI Upper`), lwd=1, position= position_dodge(width=1/2))+
            geom_pointrange(aes(y = Effect, x= Estimate, xmin=`95% CI Lower`, xmax = `95% CI Upper`), 
                            lwd = 1/2, position=position_dodge(width=1/2), shape=21)+
            scale_color_manual(values = plot.colors )+
            scale_fill_manual(values = plot.colors )+
            xlab("Proportion of Effect Mediated")+
            ylab("")+
            xlim(x.axis.lims)+
            theme_bw()+
            theme(axis.text.y=element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.title = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))  
          
          # combine plots into multipanel
          plot3<- plot_grid(plot1, 
                            plot2 + theme(legend.position = "none"),
                            ncol=1, align = "v", axis="lr", rel_heights = c(2,1))
          
          mediation_plot<-plot3
          
          
          ### DENSITY PLOTS ###
          mediation.estimate.sims<-data.frame(
            ACME.control = results$d0.sims,
            ACME.treatment = results$d1.sims,
            ACME.avg = results$d.avg.sims,
            ADE.control = results$z0.sims,
            ADE.treatment = results$z1.sims,
            ADE.avg = results$z.avg.sims,
            prop_mediated.control = results$n0.sims,
            prop_mediated.treatment = results$n1.sims,
            prop_mediated.avg = results$n.avg.sims,
            total_effect.total = results$tau.sims)%>%
            gather(effect, estimate)%>%
            separate(effect, c("effect", "group","sim"), sep="[.]")%>%
            mutate(group = plyr::revalue(group, c("control" =  paste(t.c.contrasts[[j]][2]),  "treatment" =  paste(t.c.contrasts[[j]][1]))))
          
          # plot densities
         
          plot.ACME.dens<-
            mediation.estimate.sims%>%
            filter(effect == "ACME" & group !="avg")%>%
            ggplot(aes(x = estimate, color = group, fill = group))+
            geom_density(alpha = 0.4)+
            geom_vline(xintercept = 0, linetype = "dashed")+
            xlab("Avg. Causal Mediation Effect")+
            theme_bw()+
            theme(axis.text.y=element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.title = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
            
            
          plot.ACME.ADE.dens<-
            mediation.estimate.sims%>%
            filter(effect %in% c("ACME", "ADE") & group !="avg")%>%
            ggplot(aes(x = estimate, color = group, linetype = effect))+
            geom_density(size=1)+
            geom_vline(xintercept = 0, linetype = "dashed")+
            xlab("Estimated Effect")+
            theme_bw()+
            theme(axis.text.y=element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.title = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
          
          # save plots
          mediation.plots[[
            paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
            ]][["plot.forest"]]<-mediation_plot
          mediation.plots[[
            paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]][["plot.dens.ACME"]]<-plot.ACME.dens
          mediation.plots[[
            paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]][["plot.dens.ACME.ADE"]]<-plot.ACME.ADE.dens
          
          
        }
        
        
        # Save outputs
        mediation.output[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]]<-results
        
        mediation.summary[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]]<-summary(results)
        
        mediation.summary.mat[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]]<-extractMediationSummary(results)
        
        mediation.estimate.simulations[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
        ]]<-mediation.estimate.sims
        
      
        print("Analysis complete.")
        }}} else {
  
  ###########  RUN MEDIATION IF RANDOM = FALSE ####################

    # for each gene, construct model formulas
    for(g in 1:length(gene.list)) { 
      progress(g, length(gene.list))
      i<-gene.list[g]
      
      # Specify formula for mediator model
      
      if(is.null(covariates)){med.formula<-as.formula(
        paste(i, "~", treatment, sep=" "))} else{med.formula<-as.formula(
        paste(i, "~", treatment,"+",paste(covariates, collapse=" + "), sep=" "))}
      
      # specify interaction if indicated.
      if(interaction){
        # specify formula for outcome model
        if(is.null(covariates)){out.formula<-as.formula(
          paste(paste(outcome, "~", sep=" "),
                paste(i, treatment, sep="*"), sep=" "))} else {out.formula<-as.formula(
          paste(paste(outcome, "~", sep=" "),
                paste(i, treatment, sep="*"),"+", 
                paste(covariates, collapse=" + "), sep=" "))}} else {
      
      # Specify behavior if interaction = FALSE

        if(is.null(covariates)){out.formula<-as.formula(
          paste(paste(outcome, "~", sep=" "),
                i,"+", treatment, sep=" "))} else {out.formula<-as.formula(
          paste(paste(outcome, "~", sep=" "),
                i,"+", treatment, "+",
                paste(covariates, collapse=" + "), sep=" "))
      }}
      
      # Fit mediator model
      med.fit<-lm(med.formula, model.data)
      
      # Fit outcome (dependent variable) model
      out.fit<-lm(out.formula, model.data)
      
      # Save input model output
      mediation.model.sums[[i]][[paste(i,"mediator",sep="_")]]<-summary(med.fit)
      mediation.model.sums[[i]][[paste(i,"outcome",sep="_")]]<-summary(out.fit)
      
      mediation.models[[i]][[paste(i,"mediator",sep="_")]]<-broom::tidy(med.fit)
      mediation.models[[i]][[paste(i,"outcome",sep="_")]]<-broom::tidy(out.fit)
      
      mediation.model.singularity[[i]][[paste(i,"mediator",sep="_")]]<-NA # isSingular check not relevent for class "lm"
      mediation.model.singularity[[i]][[paste(i,"outcome",sep="_")]]<-NA
      
      mediation.anovas[[i]][[paste(i,"mediator_Anova",sep="_")]]<-broom::tidy(car::Anova(med.fit))
      mediation.anovas[[i]][[paste(i,"outcome_Anova",sep="_")]]<-broom::tidy(car::Anova(out.fit))
      
      
      
      # Run moderation analysis for each set of contrasts
      for (j in 1:length(t.c.contrasts)){
        print(paste("Running Mediation Analysis |", i,"|", "Treatment:",t.c.contrasts[[j]][1],", Control:",t.c.contrasts[[j]][2], sep=" "))
        
        # Run mediation analysis.
        set.seed(seed)
        results = mediate(med.fit, out.fit, treat=treatment, mediator=i,
                          control.value = t.c.contrasts[[j]][2], 
                          treat.value = t.c.contrasts[[j]][1], boot = boot, sims = sims, boot.ci.type = boot.ci.type)
        
        #######  PRODUCE MEDIATION PLOTS ###############
        
        if(plots) { 
          # temporary internal functions for massaging the output
          temp.fun<- function(x){ifelse(stringi::stri_detect_fixed(x, "(control)"), paste(t.c.contrasts[[j]][2]), 
                                        ifelse(stringi::stri_detect_fixed(x, "(treated)"), paste(t.c.contrasts[[j]][1]), 
                                               ifelse(stringi::stri_detect_fixed(x, "(average)"), "Average", "Total")))
          }
          
          temp.fun2<- function(x){ifelse(stringi::stri_detect_fixed(x, "Prop."), "Proportion", "Effect")}
          
          
           
            # Construct plot data frame from mediation output.
            plot_df<-as.data.frame(extractMediationSummary(results))%>%
              rownames_to_column("Effect")%>%
              mutate(effect.type = temp.fun2(Effect))%>%
              mutate(group = temp.fun(Effect))%>%
              separate(Effect, sep="\\(", into=c("Effect", "temp.label"))%>%
              mutate(group=fct_relevel(group, paste(t.c.contrasts[[j]][1]), paste(t.c.contrasts[[j]][2]), "Average", "Total"))%>%
              mutate(Effect=fct_relevel(Effect, "ADE", "ACME", "Total Effect"))%>%
              suppressWarnings() # suppresses unimportant warning from separate command
            
            # Define x axis limits using extreme values from all estimate/proportion columns so that axes match 
            x.axis.lims<-plot_df%>%
              select_if(is.numeric)%>%
              dplyr::select(-`p-value`)%>%
              gather()%>%
              dplyr::select(value)%>%
              range(na.rm = T)
            
            # mediation coefficient plot
            plot1<- plot_df%>%
              filter(effect.type=="Effect")%>%
              ggplot(aes(color=group, fill=group))+
              geom_vline(xintercept = 0, color=gray(1/2), lty=2)+
              geom_linerange(aes(y=Effect, xmin=`95% CI Lower`, xmax = `95% CI Upper`), lwd=1, position= position_dodge(width=1/2))+
              geom_pointrange(aes(y = Effect, x= Estimate, xmin=`95% CI Lower`, xmax = `95% CI Upper`), 
                              lwd = 1/2, position=position_dodge(width=1/2), shape=21)+
              scale_color_manual(values = plot.colors)+
              scale_fill_manual(values = plot.colors)+
              ylab("")+
              xlim(x.axis.lims)+
              ggtitle(paste(outcome, "Mediation by", i, sep= " "))+
              theme_bw()+
              theme(plot.title = element_text(hjust = 0.5),
                    legend.title = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1))
            
            # proportion of effect explained by mediation.
            plot2<- plot_df%>%
              filter(effect.type=="Proportion")%>%
              ggplot(aes(color=group, fill=group))+
              geom_vline(xintercept = 0, color=gray(1/2), lty=2)+
              geom_linerange(aes(y=Effect, xmin=`95% CI Lower`, xmax = `95% CI Upper`), lwd=1, position= position_dodge(width=1/2))+
              geom_pointrange(aes(y = Effect, x= Estimate, xmin=`95% CI Lower`, xmax = `95% CI Upper`), 
                              lwd = 1/2, position=position_dodge(width=1/2), shape=21)+
              scale_color_manual(values = plot.colors )+
              scale_fill_manual(values = plot.colors )+
              xlab("Proportion of Effect Mediated")+
              ylab("")+
              xlim(x.axis.lims)+
              theme_bw()+
              theme(axis.text.y=element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.title = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1))  
            
            # combine plots in to multipanel
            plot3<- plot_grid(plot1, 
                              plot2 + theme(legend.position = "none"),
                              ncol=1, align = "v", axis="lr", rel_heights = c(2,1))
            
            mediation_plot<-plot3
            
            
            
            
            ### DENSITY PLOTS ###
            mediation.estimate.sims<-data.frame(
              ACME.control = results$d0.sims,
              ACME.treatment = results$d1.sims,
              ACME.avg = results$d.avg.sims,
              ADE.control = results$z0.sims,
              ADE.treatment = results$z1.sims,
              ADE.avg = results$z.avg.sims,
              prop_mediated.control = results$n0.sims,
              prop_mediated.treatment = results$n1.sims,
              prop_mediated.avg = results$n.avg.sims,
              total_effect.total = results$tau.sims)%>%
              gather(effect, estimate)%>%
              separate(effect, c("effect", "group","sim"), sep="[.]")%>%
              mutate(group = plyr::revalue(group, c("control" =  paste(t.c.contrasts[[j]][2]),  "treatment" =  paste(t.c.contrasts[[j]][1]))))
            
            # plot densities
            
            plot.ACME.dens<-
              mediation.estimate.sims%>%
              filter(effect == "ACME" & group !="avg")%>%
              ggplot(aes(x = estimate, color = group, fill = group))+
              geom_density(alpha = 0.4)+
              geom_vline(xintercept = 0, linetype = "dashed")+
              xlab("Avg. Causal Mediation Effect")+
              theme_bw()+
              theme(axis.text.y=element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.title = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1))
            
            
            plot.ACME.ADE.dens<-
              mediation.estimate.sims%>%
              filter(effect %in% c("ACME", "ADE") & group !="avg")%>%
              ggplot(aes(x = estimate, color = group, linetype = effect))+
              geom_density(size=1)+
              geom_vline(xintercept = 0, linetype = "dashed")+
              xlab("Avg. Causal Mediation Effect")+
              theme_bw()+
              theme(axis.text.y=element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.title = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1))
            
            # save plots
            mediation.plots[[
              paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
            ]][["plot.forest"]]<-mediation_plot
            mediation.plots[[
              paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
            ]][["plot.dens.ACME"]]<-plot.ACME.dens
            mediation.plots[[
              paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
            ]][["plot.dens.ACME.ADE"]]<-plot.ACME.ADE.dens
            
            
        
        }
        
        
        
        # Save outputs
        mediation.output[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]]<-results
        
        mediation.summary[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]]<-summary(results)
        
        mediation.summary.mat[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
          ]]<-extractMediationSummary(results)
        
        mediation.estimate.simulations[[
          paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
        ]]<-mediation.estimate.sims
        
 
        print("Analysis complete.")
      }
    }}
  

  
  
  ####### FORMAT OUTPUT AND SAVE TO DISK #########
  
  if(save.output){
    
    # Save mediation analysis summaries
    dir.create(paste(out.dir, "mediation_output", sep="/"), showWarnings = FALSE)
    
    for(s in 1:length(mediation.summary.mat)){
      write.csv(mediation.summary.mat[[s]],
                     paste(paste(out.dir, "mediation_output", sep="/"),
                           gsub(" ","_",gsub("/","_", paste(names(mediation.summary[s]),
                                 ".csv", sep=""))), sep="/"))
      write.csv(mediation.estimate.simulations[[s]],
                paste(paste(out.dir, "mediation_output", sep="/"),
                      gsub(" ","_",gsub("/","_", paste(names(mediation.summary[s]),
                                                       ".sims.csv", sep=""))), sep="/"))
      
    }
    
 
  ## Save input model summaries and anovas
    dir.create(paste(out.dir, "input_models", sep="/"), showWarnings = FALSE)
    
    for(m in gene.list){
      
      mod.filename= paste(paste(out.dir, "input_models", sep="/"), paste(m, "_models.csv", sep=""), sep="/")
      anova.filename= paste(paste(out.dir, "input_models", sep="/"), paste(m, "_anovas.csv", sep=""), sep="/")
      summary.filename= paste(paste(out.dir, "input_models", sep="/"), paste(m, "_model_summaries.txt", sep=""), sep="/")
      
      
      mods.df<-
        mediation.models[[m]][[1]]%>%
        mutate(model = "mediator")%>%
        mutate(DV = m)%>%
        mutate(isModelSingular = mediation.model.singularity[[m]][1])%>%
        bind_rows(
          mediation.models[[m]][[2]]%>%
            mutate(model = "outcome")%>%
            mutate(DV = outcome)%>%
            mutate(isModelSingular = mediation.model.singularity[[m]][2]))%>%
        apply(2,as.character) # flatten lists
      
      mod.anovas.df<-
        mediation.anovas[[m]][[1]]%>%
        mutate(model = "mediator")%>%
        mutate(DV = m)%>%
        mutate(isModelSingular = mediation.model.singularity[[m]][1])%>%
        bind_rows(
          mediation.anovas[[m]][[2]]%>%
            mutate(model = "outcome")%>%
            mutate(DV = outcome)%>%
            mutate(isModelSingular = mediation.model.singularity[[m]][2]))%>%
        apply(2,as.character) # flatten lists
      
      # Save model info
      write.csv(mods.df, mod.filename)
      write.csv(mod.anovas.df, anova.filename)
      
      sink(summary.filename)
      print("MEDIATOR MODEL")
      print(mediation.model.sums[[m]][1])
      print("OUTCOME MODEL")
      print(mediation.model.sums[[m]][2])
      sink()
      
    }

  
  }
  
  if(plots & save.plots){
    # Save mediation analysis plots
    dir.create(paste(plot.dir, "mediation_output_plots", sep="/"), showWarnings = FALSE)
    
    for(p in 1:length(mediation.plots)){
      
      ggsave(plot=mediation.plots[[p]][["plot.forest"]],
             filename= paste(paste(plot.dir, "mediation_output_plots", sep="/"), gsub(" ","_", gsub("/","_", paste(names(mediation.plots[p]), 
                                                                                                                   "_forest.png", sep=""))), sep="/"),
             dpi=300, height = plot.height, width = plot.width)
      
      ggsave(plot=mediation.plots[[p]][["plot.dens.ACME.ADE"]],
             filename= paste(paste(plot.dir, "mediation_output_plots", sep="/"), gsub(" ","_", gsub("/","_", paste(names(mediation.plots[p]), 
                                                                                                                   "_density.png", sep=""))), sep="/"),
             dpi=300, height = plot.height, width = plot.width)
      
    }
  }

  
  
  #################  RETURN OUTPUT ##############

  if(plots){list(mediation.output=mediation.output,
                 mediation.estimate.sims = mediation.estimate.simulations,
                 summary=mediation.summary.mat,
                 input.anovas=mediation.anovas,
                 plots=mediation.plots, 
                 contrasts = t.c.contrasts)}else{
                          list(mediation.output=mediation.output,
                               mediation.estimate.sims = mediation.estimate.sims,
                               summary=mediation.summary.mat,
                               input.anovas=mediation.anovas, 
                               contrasts = t.c.contrasts)}
}




