
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
# Contains two functions to run mediation analysis over a set of genes or gene modules
# (or any set of mediatiors, for that matter). 
#
# mediationGeneSet() runs automated model fitting and mediation analysis over a set 
#         of genes with a given outcome, model structure, and specified contrasts.
#
# extractMediationSummary() is required by mediationGeneSet() to extract mediation 
#         summary in useful matrix format.
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
#                    plot.dir = fig.dir)




########### DEFINE INPUTS ###############
mediationGeneSet<- function(model.data,
                            gene.list,
                            outcome,
                            treatment,
                            covariates,
                            interaction = FALSE,
                            t.c.contrasts,
                            save.output=FALSE,
                            out.dir,
                            plots=FALSE,
                            plot.height = 4,
                            plot.width = 6,
                            save.plots = FALSE,
                            plot.dir,
                            random=FALSE,
                            random.effect){
  
  ########## LOAD PACKAGES ############# 
  set.seed(2828)
  
  # Extract mediation summary
  source("https://github.com/rmsegnitz/Bioinformatics_Tools/tree/master/R_functions/extractMediationSummary.R")
  
  # Data manipulation and figures
  library(tidyverse)
  library(stringi)
  library(stringr)
  # mediation analysis
  library(mediation)
  # Progress tracking
  library(svMisc)

  # Plotting
  if(plots){library(cowplot)}
  
  
  ###### SETUP OUTPUT STORAGE ######
  mediation.anovas<-list()
  mediation.models<-list()
  mediation.output<-list()
  mediation.summary<-list()
  mediation.summary.mat<-list()
  mediation.plots<-list()

# Load data

  
  
###########  RUN MEDIATION IF RANDOM = TRUE ####################
  if(random){
    # load lme4 for random effects model fitting
    library(lme4)
    
    # for each gene in list, construct formulas for component models.
    for(g in 1:length(gene.list)) { 
      progress(g, length(gene.list))
      i<-gene.list[g]
      # Specify formula for mediator model
      med.formula<-as.formula(
        paste(i, "~",treatment,"+",paste(covariates, collapse=" + "), "+", paste("(1|",random.effect , ")", sep=""), sep=" "))
      
      # specify interaction if indicated.
      if(interaction){

      # specify formula for outcome model
      out.formula<-as.formula(
        paste(paste(outcome, "~", sep=" "),
              paste(i, treatment, sep="*"),"+", 
              paste(covariates, collapse=" + "),"+" , 
              paste("(1|",random.effect , ")", sep=""), sep=" "))} else{ # else no interaction if not specified.
        out.formula<-as.formula(
                  paste(paste(outcome, "~", sep=" "),
                        i,"+", treatment, "+",
                        paste(covariates, collapse=" + "),"+" , 
                        paste("(1|",random.effect , ")", sep=""), sep=" "))}
      
      # Fit mediator model
      med.fit<-lmer( med.formula, model.data)
      
      # Fit outcome (dependent variable) model
      out.fit<-lmer(out.formula, model.data)
      
      mediation.models[[paste(i,"mediator",sep="_")]]<-summary(med.fit)
      mediation.models[[paste(i,"outcome",sep="_")]]<-summary(out.fit)
      
      mediation.anovas[[paste(i,"mediator_Anova",sep="_")]]<-car::Anova(med.fit)
      mediation.anovas[[paste(i,"outcome_Anova",sep="_")]]<-car::Anova(out.fit)
      
      # Run moderation analysis for each set of contrasts
      for (j in 1:length(t.c.contrasts)){
        
        print(paste("Running Mediation Analysis |", i,"|", "Treatment:",t.c.contrasts[[j]][1],", Control:",t.c.contrasts[[j]][2], sep=" "))
        
        # run mediation analysis
        results = mediate(med.fit, out.fit, treat=treatment, mediator=i,
                          control.value = t.c.contrasts[[j]][2], 
                          treat.value = t.c.contrasts[[j]][1])
        
        #######  PRODUCE MEDIATION PLOTS ###############
        
        if(plots) { 
          # temporary internal functions for massaging the output
          temp.fun<- function(x){ifelse(stringi::stri_detect_fixed(x, "(control)"), paste(t.c.contrasts[[j]][2]), 
                                        ifelse(stringi::stri_detect_fixed(x, "(treated)"), paste(t.c.contrasts[[j]][1]), 
                                               ifelse(stringi::stri_detect_fixed(x, "(average)"), "Average", "Total")))
          }
          
          temp.fun2<- function(x){ifelse(stringi::stri_detect_fixed(x, "Prop."), "Proportion", "Effect")}
          
          
          
          # create plotting data frame from mediation output 
          plot_df<-as.data.frame(extractMediationSummary(results))%>%
            rownames_to_column("Effect")%>%
            mutate(effect.type = temp.fun2(Effect))%>%
            mutate(group = temp.fun(Effect))%>%
            separate(Effect, sep="\\(", into=c("Effect", "temp.label"))%>%
            mutate(group=fct_relevel(group, paste(t.c.contrasts[[j]][1]), paste(t.c.contrasts[[j]][2]), "Average", "Total"))%>%
            mutate(Effect=fct_relevel(Effect, "ADE", "ACME", "Total Effect"))
          
          # Define x axis range (use extreme values for all coef/proportion columns so that axes match)
          x.axis.lims<-plot_df%>%
            select_if(is.numeric)%>%
            dplyr::select(-`p-value`)%>%
            gather()%>%
            dplyr::select(value)%>%
            range(na.rm = T)
          
          # Plot the mdeiation coefficient comparisons
          plot1<- plot_df%>%
            filter(effect.type=="Effect")%>%
            ggplot(aes(color=group, fill=group))+
            geom_vline(xintercept = 0, color=gray(1/2), lty=2)+
            geom_linerange(aes(y=Effect, xmin=`95% CI Lower`, xmax = `95% CI Upper`), lwd=1, position= position_dodge(width=1/2))+
            geom_pointrange(aes(y = Effect, x= Estimate, xmin=`95% CI Lower`, xmax = `95% CI Upper`), 
                            lwd = 1/2, position=position_dodge(width=1/2), shape=21)+
            scale_color_colorblind()+
            scale_fill_colorblind()+
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
            scale_color_colorblind()+
            scale_fill_colorblind()+
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
          
          # save plot
          mediation.plots[[
            paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
            ]]<-mediation_plot
          
          
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
        
      
        print("Analysis complete.")
        }}} else {
  
  ###########  RUN MEDIATION IF RANDOM = FALSE ####################

    # for each gene, construct model formulas
    for(g in 1:length(gene.list)) { 
      progress(g, length(gene.list))
      i<-gene.list[g]
      
      # Specify formula for mediator model
      med.formula<-as.formula(
        paste(i, "~", treatment,"+",paste(covariates, collapse=" + "), sep=" "))
      
      # specify interaction if indicated.
      if(interaction){
        # specify formula for outcome model
        out.formula<-as.formula(
          paste(paste(outcome, "~", sep=" "),
                paste(i, treatment, sep="*"),"+", 
                paste(covariates, collapse=" + "), sep=" "))} else {
      
      # Specify behavior if interaction = FALSE

        out.formula<-as.formula(
          paste(paste(outcome, "~", sep=" "),
                i,"+", treatment, "+",
                paste(covariates, collapse=" + "), sep=" "))
      }
      
      # Fit mediator model
      med.fit<-lm(med.formula, model.data)
      
      # Fit outcome (dependent variable) model
      out.fit<-lm(out.formula, model.data)
      
      mediation.anovas[[paste(i,"mediator_Anova",sep="_")]]<-car::Anova(med.fit)
      mediation.anovas[[paste(i,"outcome_Anova",sep="_")]]<-car::Anova(out.fit)
      
      # Run moderation analysis for each set of contrasts
      for (j in 1:length(t.c.contrasts)){
        print(paste("Running Mediation Analysis |", i,"|", "Treatment:",t.c.contrasts[[j]][1],", Control:",t.c.contrasts[[j]][2], sep=" "))
        
        # Run mediation analysis.
        results = mediate(med.fit, out.fit, treat=treatment, mediator=i,
                          control.value = t.c.contrasts[[j]][2], 
                          treat.value = t.c.contrasts[[j]][1])
        
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
              mutate(Effect=fct_relevel(Effect, "ADE", "ACME", "Total Effect"))
            
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
              scale_color_colorblind()+
              scale_fill_colorblind()+
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
              scale_color_colorblind()+
              scale_fill_colorblind()+
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
            
            # save plot
            mediation.plots[[
              paste(i,paste("T",t.c.contrasts[[j]][1], sep="-"),paste("C",t.c.contrasts[[j]][2], sep="-"), sep="_")
              ]]<-mediation_plot
            
        
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
        
 
        print("Analysis complete.")
      }
    }}
  

  
  
  ####### FORMAT OUTPUT AND SAVE TO DISK #########
  
  if(save.output){
    
    # Save mediation analysis summaries
    dir.create(paste(out.dir, "mediation_output", sep="/"), showWarnings = FALSE)
    
    for(s in 1:length(mediation.summary.mat)){
      write.csv(mediation.summary.mat[s],
                     paste(paste(out.dir, "mediation_output", sep="/"), 
                           gsub(" ","_",gsub("/","_", paste(names(mediation.summary[s]), 
                                 ".csv", sep=""))), sep="/"))
    }
    
 
  ## Save input model summaries and anovas
    dir.create(paste(out.dir, "input_models", sep="/"), showWarnings = FALSE)
    
    for(m in 1:length(mediation.models)){
      filename= paste(paste(out.dir, "input_models", sep="/"), paste(names(mediation.models)[[m]], ".txt", sep=""), sep="/")
      # Setup capture file
      cat(paste("Input Model Summary & ANOVA:", names(mediation.models)[[m]], sep=" "), file=filename)
      # add 2 newlines
      cat("\n\n", file = filename, append = TRUE)
      # export anova test output
      cat("SUMMARY\n", file = filename, append = TRUE)
      capture.output(summary(mediation.models[[m]]), file = filename, append = TRUE)
      # add 2 newlines
      cat("\n\n", file = filename, append = TRUE)
      # export anova test output
      cat("ANOVA\n", file = filename, append = TRUE)
      capture.output(print(mediation.anovas[[m]]), file = filename, append = TRUE)

    }

  
  }
  
  if(plots & save.plots){
    # Save mediation analysis plots
    dir.create(paste(plot.dir, "mediation_output_plots", sep="/"), showWarnings = FALSE)
    
    for(p in 1:length(mediation.plots)){
      
      ggsave(plot=mediation.plots[[p]],
        filename= paste(paste(plot.dir, "mediation_output_plots", sep="/"), gsub(" ","_", gsub("/","_", paste(names(mediation.plots[p]), 
                                                                                                ".png", sep=""))), sep="/"),
        dpi=300, height = plot.height, width = plot.width)
    }
  }

  
  
  #################  RETURN OUTPUT ##############

  if(plots){list(output=mediation.output,
                 summary=mediation.summary.mat,
                 input.anovas=mediation.anovas,
                 plots=mediation.plots)}else{
                          list(output=mediation.output,
                                summary=mediation.summary.mat,
                                  input.anovas=mediation.anovas)}
}




#'########################################
####  2)  extractMediationSummary     ####
#'########################################

## This function comes from [https://stackoverflow.com/questions/41582486/how-to-convert-r-mediation-summary-to-data-frame ] 
## and is not my own. It has been modified in the name of the function only.
## Original Author: Stack Overflow user "hrbrmstr"
## This material is covered in Creative Commons, 

### DEFINE FUNCTION ####

extractMediationSummary <- function (x) { 
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                   (inherits(x$model.y, "glm") && x$model.y$family$family == 
                      "gaussian" && x$model.y$family$link == "identity") || 
                   (inherits(x$model.y, "survreg") && x$model.y$dist == 
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                        "ADE (control)", "ADE (treated)", "Total Effect", 
                        "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
  
}



