

load("~/Documents/WORK/UW/ALTMAN/GITHUB/GRADUATE_batch2/data_clean/P662_2_voom_meta.RData")

targets<-dat_meta$targets

geneLME_contrast_spec(targets, contrast_vars = c("treatment", "TNSS_AUC"))


library(emmeans)
library(lme4)
target_NAC<-dplyr::filter(targets, visit_nac)

test_mod<-lmer(TNSS_AUC ~ treatment*visit + (1|ptID), data=target_NAC)

contrast(emmeans(test_mod, pairwise ~ treatment | visit), "pairwise")

mvcontrast(emmeans(test_mod, ç), mult.name = c("visit"))


test_mod<-lmer(TNSS_AUC ~ age*treatment + (1|ptID), data=target_NAC)
emtrends(test_mod, "treatment", var = "age" )

pairs(contrast(emmeans(test_mod, ~ treatment | age), at = list(age =c (35, 40, 45))))

emmeans(test_mod, pairwise ~ treatment | age, at = list(age = c(35, 40, 45)))


# case where discrete*continuous and no levels given. 
# behavior: calculates the pairwise discrete contrasts at the standardard quartiles of the continuous
test_mod<-lmer(TNSS_AUC ~ age*treatment + (1|ptID), data=target_NAC)

test_out_disc_cont_emm<-
emmeans(test_mod, pairwise ~ treatment | age, at = list(age = quantile(test_mod@frame$age)), adjust="none")

test_out_disc_cont_contrasts<-
  as.data.frame(test_out_disc_cont_emm$contrasts)


# case where continuouse by discrete and no levels given. 
# behavior: calculates the pairwise discrete contrasts at the standardart quartiles of the continuous
test_mod<-lmer(TNSS_AUC ~ age*treatment + (1|ptID), data=target_NAC)

test_out_cont_disc_emtrend<-
  emtrends(test_mod,   "treatment", var= "age")

# second order by default returns pairs
test_out_cont_disc_emtrend_secondOrder<-
  pairs(test_out_cont_disc_emtrend, adjust="none")






# toy plots
library(tidyverse)
test_out_disc_cont_contrasts%>%
mutate(plot_sig=ifelse(p.value<0.05,"p<0.05", "p>=0.05"))%>%
ggplot(data=., aes(x=age, y=estimate))+
  geom_errorbar(aes(ymin=estimate-SE, ymax=estimate+SE))+
  geom_point(color="white")+
  geom_point(aes(shape=plot_sig))+
  scale_shape_manual(values=c("p<0.05"=16, "p>=0.05"=1))+
  theme_bw()+
  facet_wrap(~contrast)

