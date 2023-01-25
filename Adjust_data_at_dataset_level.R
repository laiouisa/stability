#Analyses for resubmission

#packages-------
library(tidyverse)#managing data
library(readxl) #reading excel files
library(lme4)#Mixed Effects models
library(MuMIn) #r.squared + confidence intervals with visreg (replaces predict.merMod?)
library(car) #for vif
#for correlations
library(psych)
library(corrplot)
#plots
library(ggplot2) 
library(coefplot) 
library(visreg)
library(sjPlot)
library(ggpubr)
library(RColorBrewer)
# library(sjmisc)
library(jtools)
library(huxtable)
library(lmerTest)
library(MASS)
library(cAIC4)

source("/Users/laiouisa/OneDrive - ČZU v Praze/R/R scripts/facetzoom2.R") #for facet zoom on the right 


#DATA------
load("/Users/laiouisa/OneDrive - ČZU v Praze/STABILITY Project/Stability_and_traits/final matrix.RData")
stability_selected<-as_tibble(stability_selected)
names(stability_selected)

#datasets names
ids_datasets<-read_xlsx("/Users/laiouisa/OneDrive - ČZU v Praze/STABILITY Project/New Phytologist/sup mat/Table S1 Datasets information.xlsx", sheet = "ids")
ids_datasets<-read_xlsx("C:/Users/conti/OneDrive - ČZU v Praze/STABILITY Project/New Phytologist/sup mat/Table S1 Datasets information.xlsx", sheet = "ids")
class(ids_datasets$Original_dataset)

# load("C:/Users/conti/OneDrive - ČZU v Praze/STABILITY Project/Stability_and_traits/boosted regression trees.RData")

# Detrended CV vs CV-----
scatter.smooth(stability_selected$sp_cv,stability_selected$sp_cv_t3)
hist(stability_selected$sp_cv)
hist(stability_selected$sp_cv_t3)
boxplot(stability_selected$sp_cv, stability_selected$sp_cv_t3, names = c("CV", "CV detrended"),
        main="All species")
plot(stability_selected$sp_cv~abs(stability_selected$cor_time_abu))
plot(stability_selected$sp_cv_t3~abs(stability_selected$cor_time_abu))

#how the different indexes relate to the abundance in time
names(stability_selected)
mod_cv<-lmer(sp_cv~abs(cor_time_abu) + (1|species) + (1|original_dataset), data=stability_selected)
mod_cvt3<- lmer(sp_cv_t3~abs(cor_time_abu) + (1|species) + (1|original_dataset), data=stability_selected)
summary(mod_cv) #POSITIVE correlation with abundance in time
summary(mod_cvt3) #NEGATIVE correlation with abundance in time
c<-coefplot(mod_cv, intercept = F, title = 'Normal CV', 
         ylab = '', xlab = '', newNames = c('abs(cor_time_abu)'='Abundance in time') )

coefplot(mod_cvt3, intercept = F, title = 'Detrended CV', 
         ylab = '', xlab = '', newNames = c('abs(cor_time_abu)'='Abundance in time'))
r.squaredGLMM(mod_cv) #R2m is 0.00005
r.squaredGLMM(mod_cvt3) #R2m is 0.06 - detrended is negatively related to abundance in time!?

names(stability_selected)

#get better plots?
visreg(mod_cv)
visreg(mod_cvt3)

#this not working with abs()
# mod_dat <- get_model_data(mod_cv, type = "pred", terms = "abs(cor_time_abu) [all]", pretty = FALSE)
# c<-ggplot(aes(y=log_sp_mean_abu_rel, x=mean_LDMC), data = stability_selected[woodiness_new=="non-woody"&
#                                                                                !Climate_Change==1&   
#                                                                                !Sugar_Addition==1&
#                                                                                !Burial==1&
#                                                                                !Shelter==1,])
# c<-c+ geom_point( colour="gray50", show.legend = F) + #geom_abline(slope = 1.3697, intercept = -4.9615, col="red", lwd=1.5) 
#   geom_ribbon(data = mod_dat, aes(x = x,  ymin = conf.low, ymax=conf.high, fill = "red"), 
#               inherit.aes = F, show.legend = F, alpha = .25)+
#   geom_line(data = mod_dat, aes(x = x, y = predicted, group = group, colour = "red"), 
#             size = 1.025, show.legend = F)
# c<-c+ theme_classic(base_size=22) + labs(y="LOG Species relative abundance", x="LDMC [g/g]")
# c+ annotate("text", x= 0.55, y = -13 ,  label = "paste(italic(R) ^ 2, \"m= 0.004\")", parse=T, size=6)+ 
#   annotate("text", x= 0.55, y = -14 ,  label = "paste(italic(R) ^ 2, \"c=0.602\")", parse=T, size=6)


#violin plot (instead of boxplot) - a bit better
ggplot(data = stability_selected%>%filter(!Climate_Change==1,
                                          !Sugar_Addition==1,
                                          !Sugar_Addition==1, 
                                          !Burial==1,
                                          !Shelter==1)%>%select(sp_cv,sp_cv_t3)%>%
         pivot_longer(sp_cv:sp_cv_t3,names_to = "index", values_to = "value"), 
       aes(y=value, x=index))+
  geom_violin(na.rm = T, trim = F) + 
  # geom_point(size = 0.5, colour="gray", show.legend = F, position = 'jitter')+
  theme_classic()+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")
                                                                                                                                                 

#Average CV values of species in each dataset - USING NORMAL CV-------
names(stability_selected)
summary(stability_selected$Burial==1)

# stability_selected_mean<-stability_selected%>%
#   group_by(species, dataset) %>%
#   select(sp_cv, sp_cv_t3) %>%
#   summarise(
#     mean_sp_cv = mean(sp_cv, na.rm = TRUE),
#     mean_sp_cv_t3 = mean(sp_cv_t3, na.rm = TRUE), .groups='keep'
#   )

#Already done for regression trees (z.scores at species-dataset level)!!! - but woody spp are missing
mean_stability<-as_tibble(mean_stability)
names(mean_stability)
mean_stability[,'residency_prop']

#add the woody species
names(mean_stability)
summary(mean_stability$woodyness)
names(mean_stability_wood)
mean_stability_wood <- as_tibble(mean_stability_wood)
names(mean_stability_wood) #categorical traits missing - adjust 
mean_stability_all<-full_join(mean_stability, as_tibble(mean_stability_wood))
names(mean_stability_all)

#woodyness
summary(mean_stability_all$woodyness_new)
mean_stability_all$woodyness_new[is.na(mean_stability_all$woodyness_new)]<-'woody'
mean_stability_all$woodyness[is.na(mean_stability_all$woodyness)]<-'woody'

#lifespan
summary(mean_stability_all$lifespan_new)
mean_stability_all$lifespan_new[is.na(mean_stability_all$lifespan_new)]<-'not-annual'

#lifeform
summary(mean_stability_all$lifeform)
mean_stability_all$lifeform[is.na(mean_stability_all$lifeform)]<-as.factor(LifeForm_average$LifeForm_AVE[match(mean_stability_all$Group.1[is.na(mean_stability_all$lifeform)], LifeForm_average$Species)])

#growthform
summary(mean_stability_all$growthform)
mean_stability_all$growthform[is.na(mean_stability_all$growthform)]<-GrowthForm_average$growthform_AVE[match(mean_stability_all$Group.1[is.na(mean_stability_all$growthform)], GrowthForm_average$Species)]


#<---------------->#########
# GLOBAL MIXED EFFECT MODEL ---------
#explaining species stability 
#at dataset level (averaging CV values of species in each dataset)
#with all traits considered as fixed effect and species nested in 
#datasets as random effect, only non-woody species(?)
names(mean_stability)

mean_stability%>%filter(woodyness_new=='non-woody')
summary(duplicated(mean_stability$Group.1))

attach(mean_stability)
hist(scale(mean_height))
hist(scale(log(mean_height)))
mean_stability$z.log.mean_height<-scale(log(mean_height))
hist(scale(mean_LeafN))
mean_stability$z.mean_LeafN<-scale(mean_LeafN)
hist(scale(mean_LeafP))
mean_stability$z.mean_LeafP<-scale(mean_LeafP)
hist(scale(log(mean_SeedMass)))
mean_stability$z.log.mean_SeedMass<-scale(log(mean_SeedMass))
hist(scale(log(mean_SLA)))
mean_stability$z.log.mean_SLA<-scale(log(mean_SLA))
hist(scale(mean_LDMC))
mean_stability$z.mean_LDMC<-scale(mean_LDMC)
hist(scale(mean_SSD))
mean_stability$z.mean_SSD<-scale(mean_SSD)
unique(growthform)
unique(woodyness)

#++Additive model quantitative traits----------------
#only non-woody species

names(mean_stability)
mean_stability[,"Group.1"] #species
mean_stability[,"Group.2"] #dataset


mod_mean_z.cv.halfFULL<- lmer(z.sp_cv ~ z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                              z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD +
                          (1|Group.1) + (1|Group.2), 
                        data= dataframe.non_woody, REML = F)

dataframe.non_woody<-mean_stability%>%filter(woodyness_new=='non-woody')

summary(mod_mean_z.cv.halfFULL)
r.squaredGLMM(mod_mean_z.cv.halfFULL) #0.0488 0009 0.1278361
vif(mod_mean_z.cv.halfFULL) 
# 1 = not correlated.
# Between 1 and 5 = moderately correlated.
# Greater than 5 = highly correlated.
# z.log.mean_height        z.mean_LeafN        z.mean_LeafP z.log.mean_SeedMass      z.log.mean_SLA 
# 1.072697            2.142954            1.511217            1.082380            1.559278 
# z.mean_LDMC          z.mean_SSD 
# 1.542247            1.189606 


c<-coefplot(mod_mean_z.cv.halfFULL, intercept = F, title = "z.scores")
c+ geom_text(x=-0.2, y=7, label="R2m=0.05\nR2c=0.13", size=2)


summary(mod_mean_z.cv.halfFULL)
c<-coefplot(mod_mean_z.cv.halfFULL, intercept = F, title= "normal")

c+ geom_text(x=-0.1, y=7, label="R2m=0.04\nR2c=0.39", size=2)

r.squaredGLMM(mod_mean_z.cv.halfFULL) #0.04009069 0.3937302
AIC(mod_mean_z.cv.halfFULL) #1939.179

#and then simplifying the model - 
#This would replace the boosted regression trees,
#and would allow us to define the main traits explaining the temporal variability of species, regardless of the treatment.

#method to test variables? 
drop1(mod_mean_z.cv.halfFULL , trace = T , na.action=na.omit) #not working 
# Error in drop1.merMod(mod_mean_cv.halfFULL, trace = T, na.action = na.omit) : 
#   number of rows in use has changed: remove missing values?

#compare manually?
coefplot(mod_mean_z.cv.halfFULL)
r.squaredGLMM(mod_mean_z.cv.halfFULL) 
AIC(mod_mean_z.cv.halfFULL) #1939.179 #cannot compare it because of NAs!
#scale(mean_LeafP) scale(mean_SSD) 

mod_mean_z.cv.halfFULL@call$formula

# R2m       R2c
# [1,] 0.04880009 0.1278361

r.squaredGLMM(update(mod_mean_z.cv.halfFULL, .~. - z.log.mean_height)) # higher -  out
# R2m      R2c
# [1,] 0.04883321 0.128101 #3

r.squaredGLMM(update(mod_mean_z.cv.halfFULL, .~. -   z.mean_LeafN )) #higher
# R2m       R2c
# [1,] 0.0482846 0.1288344 #5

r.squaredGLMM(update(mod_mean_z.cv.halfFULL, .~. - z.mean_LeafP)) #higher - out
# R2m       R2c
# [1,] 0.04958985 0.1258183 #2

r.squaredGLMM(update(mod_mean_z.cv.halfFULL, .~. - z.mean_SSD)) #much higher - out
# R2m       R2c
# [1,] 0.05445529 0.1541864  #1

r.squaredGLMM(update(mod_mean_z.cv.halfFULL, .~. - z.log.mean_SLA)) #higher
# R2m       R2c
# [1,] 0.04870645 0.1277245 #4


#final model------
mod_mean_z.cv<-update(mod_mean_z.cv.halfFULL, .~. -   z.mean_SSD) #1
r.squaredGLMM(mod_mean_z.cv)
r.squaredGLMM(update(mod_mean_z.cv, .~. -  z.mean_LeafP)) #2
mod_mean_z.cv<-update(mod_mean_z.cv, .~. -  z.mean_LeafP)
r.squaredGLMM(mod_mean_z.cv)
mod_mean_z.cv<-update(mod_mean_z.cv, .~. -  z.log.mean_height) #3
r.squaredGLMM(mod_mean_z.cv)
r.squaredGLMM(update(mod_mean_z.cv, .~. -  z.log.mean_SLA)) #4 - LOWER
r.squaredGLMM(update(mod_mean_z.cv, .~. -  z.mean_LeafN)) #5 - LOWER

# final model
summary(mod_mean_z.cv)
r.squaredGLMM(mod_mean_z.cv)
# R2m       R2c
# 0.06957713 0.1824155

#variables are not entirely independent?
vif(mod_mean_z.cv) #moderately correlated 
# z.mean_LeafN z.log.mean_SeedMass      z.log.mean_SLA         z.mean_LDMC 
# 1.402810            1.061450            1.306819            1.232314 

#cor test?? - to improve 
cor.test(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LeafN)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.log.mean_SLA)))
plot(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LeafN)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.log.mean_SLA)))

cor.test(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LDMC)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.log.mean_SLA)))
plot(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LDMC)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.log.mean_SLA)))

cor.test(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LDMC)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.mean_LeafN)))
plot(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LDMC)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.mean_LeafN)))

cor.test(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LDMC)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.log.mean_SeedMass)))
plot(unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.mean_LDMC)), unlist(mean_stability%>%filter(woodyness_new=='non-woody')%>%select( z.log.mean_SeedMass)))


#compare AIC for full and reduced model --------------
#fit reduced with dataset from full
mod_mean_z.cv.halfFULL #full model
mod_mean_z.cv #reduced model

frame.full<-mod_mean_z.cv.halfFULL@frame

mod_reduced.frame.full<- lmer(z.sp_cv ~ z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA +     
                                z.mean_LDMC + (1 | Group.1) + (1 | Group.2),
                              data= frame.full, REML = F)
summary(mod_reduced.frame.full) #AIC  1934.6
summary(mod_mean_z.cv.halfFULL) #AIC 1939.2

#stepAIC by hand - but there are still NAs??
AIC(mod_mean_z.cv.halfFULL) #1939.179
AIC(update(mod_mean_z.cv.halfFULL, .~. -   z.mean_SSD)) #3980.805 higher
AIC(update(mod_mean_z.cv.halfFULL, .~. -  z.mean_LeafP)) #1969.341 higher
AIC(update(mod_mean_z.cv.halfFULL, .~. -  z.log.mean_height)) #1937.189 #lowest
AIC(update(mod_mean_z.cv.halfFULL, .~. -  z.log.mean_SLA)) #1937.233 - lower
AIC(update(mod_mean_z.cv.halfFULL, .~. -  z.mean_LeafN)) #1937.324 - lower
AIC(update(mod_mean_z.cv.halfFULL, .~. -  z.mean_LDMC)) #1960.271 - higher
AIC(update(mod_mean_z.cv.halfFULL, .~. -  z.log.mean_SeedMass)) #1951.56 - higher

mod_step.aic<-update(mod_mean_z.cv.halfFULL, .~. -  z.log.mean_height)
AIC(mod_step.aic) #1937.189
AIC(update(mod_step.aic, .~. -  z.mean_SSD)) #3988.097
AIC(update(mod_step.aic, .~. -  z.mean_LeafP)) #1967.398
AIC(update(mod_step.aic, .~. -  z.log.mean_SLA)) #1935.24 #lowest
AIC(update(mod_step.aic, .~. -  z.mean_LeafN)) #1935.331 
AIC(update(mod_step.aic, .~. -  z.mean_LDMC)) #1958.47
AIC(update(mod_step.aic, .~. -  z.log.mean_SeedMass)) #1949.6

mod_step.aic<-update(mod_step.aic, .~. -  z.log.mean_SLA)
AIC(mod_step.aic) #1935.24
AIC(update(mod_step.aic, .~. -  z.mean_SSD)) #3987.352
AIC(update(mod_step.aic, .~. -  z.mean_LeafP)) #1965.502
AIC(update(mod_step.aic, .~. -  z.mean_LeafN)) #1933.487 #lowest
AIC(update(mod_step.aic, .~. -  z.mean_LDMC)) #1956.708
AIC(update(mod_step.aic, .~. -  z.log.mean_SeedMass)) #1947.727

mod_step.aic<-update(mod_step.aic, .~. -  z.mean_LeafN)
AIC(mod_step.aic) #1933.487
AIC(update(mod_step.aic, .~. -  z.mean_SSD)) #3990.761
AIC(update(mod_step.aic, .~. -  z.mean_LeafP)) #1980.44
AIC(update(mod_step.aic, .~. -  z.mean_LDMC)) #1957.914
AIC(update(mod_step.aic, .~. -  z.log.mean_SeedMass)) #1946.203
summary(mod_step.aic)

#tables-----
mod_mean_z.cv

tab<-export_summs(mod_mean_z.cv.halfFULL, mod_mean_z.cv, 
                  note=NULL) # results = 'asis'
#p values calculated using Satterthwaite d.f.
quick_xlsx(tab, file='table.xlsx')
# tab$names[1]<-"GEN_MPD_Multi_trait"
number_format(tab)<-gsub("%.2f", "%.4f", number_format(tab))
# tab_final<-tab
# 
# tab_final<- add_rows(tab_final, tab)

#AIC comparison table
tab_aic<-export_summs(mod_mean_z.cv.halfFULL, mod_reduced.frame.full, 
                  note=NULL) # results = 'asis'
#p values calculated using Satterthwaite d.f.
quick_xlsx(tab_aic, file='table_aic.xlsx')
# tab$names[1]<-"GEN_MPD_Multi_trait"
number_format(tab)<-gsub("%.2f", "%.4f", number_format(tab))

summary(mod_reduced.frame.full) #AIC  1934.6
summary(mod_mean_z.cv.halfFULL) #AIC 1939.2


#coefficent plot------


#How wide the confidence interval should be, normally inner=1 standard deviation, outer=2 sd
#on a normally distributed variable this corresponds to 68.3% and 95.5% of data, respectively

summary(mod_mean_z.cv)
c.mean_z.cv<-coefplot(mod_mean_z.cv, intercept = F, title='', color="black",
         pointSize = 5, lwdInner = 3, lwdOuter = 1.5, sort = 'natural')
c.mean_z.cv<-c.mean_z.cv+  labs(y="", x="Coefficent") + theme_classic() +
  scale_y_discrete(labels=c("z.mean_LeafN"="Leaf N content", 
                            "z.log.mean_SeedMass"="Seed Mass",
                            "z.log.mean_SLA"="SLA",
                            "z.mean_LDMC"="LDMC"))
# c_sm<-c_sm + xlim(-0.1,0.2)+ theme_classic() 




#effect of abundance
mod_mean_mix<- lmer(z.sp_cv ~  z.mean_LeafN +    
                       z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + scale(sp_mean_abu_rel) +
                       (1|Group.1) + (1|Group.2), 
                     data= mean_stability%>%filter(woodyness_new=='non-woody'), REML = F)

coefplot(mod_mean_mix)





#effect plot-----
names(mod_mean_z.cv@frame)
visreg(mod_mean_z.cv)

#Do it by hand

names(mean_stability)
data_full<-mean_stability%>%filter(woodyness_new=='non-woody')%>%select(Group.1, Group.2,  z.sp_cv,z.mean_LeafN,  z.log.mean_SeedMass, z.log.mean_SLA, z.mean_LDMC )
names(data_full)
str(mod_mean_z.cv)
data_z.cv<- mod_mean_z.cv@frame
# data_fens$dataset<-'fens'

mod_dat_z.cv <- get_model_data(mod_mean_z.cv, type = "pred", pretty = FALSE)

mod_dat_z.cv$z.mean_LDMC

names(data_full)
str(data_z.cv)

plot(data_full$z.sp_cv~data_full$z.mean_LDMC)
plot(data_z.cv$z.sp_cv~data_z.cv$z.mean_LDMC)

#create mean cv for each species (395 species)
data_z.cv<-mutate(data_z.cv%>%group_by(Group.1), mean_z.cv=mean(z.sp_cv))

#LDMC
str(mod_dat_z.cv$z.mean_LDMC)

p1<-ggplot(data =data_z.cv, aes(y = z.sp_cv, x = z.mean_LDMC))+
  xlab(label = "LDMC")+
  ylab(label = "CV")+
  theme_classic(base_size = 15)+
  geom_point(col=gray(0.25, 0.5))+
  # geom_point( aes(y = mean_z.cv, x = z.mean_LDMC), col='black') +
  #geom_smooth(method = "glm")
  # geom_ribbon(data = mod_dat_z.cv$z.mean_LDMC, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              # inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_z.cv$z.mean_LDMC, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, colour='red')
  # geom_text(x=0.35, y=100, label="R2m=0.09\nR2c=0.45", size=10)
p1

#SLA
str(mod_dat_cv$z.log.mean_SLA)

p2<-ggplot(data =data_z.cv, aes(y = z.sp_cv, x = z.log.mean_SLA))+
  xlab(label = "SLA")+
  ylab(label = "CV")+
  theme_classic(base_size = 15)+
  #fens
  geom_point(col=gray(0.25, 0.5))+
  # geom_point( aes(y = mean_z.cv, x = z.log.mean_SLA), col='black') +
  #geom_smooth(method = "glm")
  # geom_ribbon(data = mod_dat_z.cv$z.log.mean_SLA, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              # inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_z.cv$z.log.mean_SLA, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, colour='red')
# geom_text(x=0.35, y=100, label="R2m=0.09\nR2c=0.45", size=10)
p2

#SM
str(mod_dat_cv$z.log.mean_SeedMass)

p3<-ggplot(data =data_z.cv, aes(y = z.sp_cv, x = z.log.mean_SeedMass))+
  xlab(label = "Seed Mass")+
  ylab(label = "CV")+
  theme_classic(base_size = 15)+
  #fens
  geom_point(col=gray(0.25, 0.5))+
  # geom_point( aes(y = mean_z.cv, x = z.log.mean_SeedMass), col='black') +
  #geom_smooth(method = "glm")
  # geom_ribbon(data = mod_dat_z.cv$z.log.mean_SeedMass, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              # inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_z.cv$z.log.mean_SeedMass, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, colour='red')
# geom_text(x=0.35, y=100, label="R2m=0.09\nR2c=0.45", size=10)
p3


#Leaf N
str(mod_dat_cv$z.mean_LeafN)

p4<-ggplot(data =data_cv, aes(y = z.sp_cv, x = z.mean_LeafN))+
  xlab(label = "Leaf N content")+
  ylab(label = "CV")+
  theme_classic(base_size = 15)+
  #fens
  geom_point(col=gray(0.25, 0.5))+
  # geom_point( aes(y = mean_cv, x = z.mean_LeafN), col='black') +
  #geom_smooth(method = "glm")
  # geom_ribbon(data = mod_dat_cv$z.mean_LeafN, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              # inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_cv$z.mean_LeafN, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, colour='red')
# geom_text(x=0.35, y=100, label="R2m=0.09\nR2c=0.45", size=10)
p4


#put together using ggarrange
quartz(width = 7, height = 4) 
ggarrange(p1,p2,p3, p4, labels = 'auto', ncol=2, nrow=2)



#Caterpillar plots---------
#only random intercept - does it makes sense??
mod_mean_z.cv
plot_model(mod_mean_z.cv, type='re') #ROSE_ECN alwayes the weird one


#do it by hand
ranef(mod_mean_z.cv)$Group.1 #species (395)
ranef(mod_mean_z.cv)$Group.2 #Datasets (77) 

#in the model there is one dataset less than in the original data
randoms_mod_mean<-ranef(mod_mean_z.cv, condVar = TRUE)$Group.2
#condVar - an optional logical argument indicating if the conditional variance-covariance matrices 
#of the random effects should be added as an attribute.
#check Dataset names
qq_mod_mean <- attr(ranef(mod_mean_z.cv, condVar = TRUE)[[2]], "postVar")
qq_mod_mean[1,1,1:77] #for random intercept
# qq[2,2,1:77] #for random slope
# rand.interc<-randoms
df_mod_mean<-data.frame(Intercepts=randoms_mod_mean[,1],
                        sd.interc= 2*sqrt(qq_mod_mean[1,1,1:77]),
                        lev.names=rownames(randoms_mod_mean))
# df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)]) #order intercepts when plotting
# df$lev.names<-factor(df$lev.names, levels = rev(levels(df$lev.names)))


#CHECK NAMES
row.names(ranef(mod_mean_z.cv)$Group.2)
names_datasets<-read_xlsx("/Users/laiouisa/OneDrive - ČZU v Praze/STABILITY Project/Dataset info.xlsx" ,sheet = 3)

# setdiff(row.names(ranef(mod_mean_cv)$Group.2), names_datasets$Original_dataset)
# row.names(ranef(mod_mean_cv)$Group.2)<-gsub("Penuelas_Garraf_Peńuelas","Penuelas_Garraf_Peñuelas",  row.names(ranef(mod_mean_cv)$Group.2))

# library(readxl)
# LOTVS_traits_coverage_Aug2019 <- read_excel("/Users/laiouisa/OneDrive - ČZU v Praze/STABILITY Project/LOTVS_traits_coverage-Aug2019.xlsx", 
#                                             sheet = "Species_prop")
# LOTVS_traits_coverage_Aug2019$Original_dataset
# LOTVS_traits_coverage_Aug2019$Standardized_dataset_name[60]<-"Schutz_SNP"

# df_mod_mean$lev.names<-as.character(df_mod_mean$lev.names)
df_mod_mean$lev.names[51]<-"Penuelas_Garraf_Peñuelas" 
# df_mod_mean$lev.names[61]<-"Schütz_SNP" 
# df_mod_mean$lev.names[68]<-"ShmidtGöttingerForest" 

df_mod_mean$new.names<- names_datasets$Dataset_Name[match(df_mod_mean$lev.names, names_datasets$Original_dataset)]
# missing: 43, 74 
names_datasets[c(43,74),]


# df_mod_mean$new.names<-gsub("_", " ", df_mod_mean$new.names)

df_mod_mean$new.names<-factor(df_mod_mean$new.names, levels = rev(levels(as.factor(df_mod_mean$new.names))))

# #colors for positive or negative Intercepts
# summary(mod_mean_cv)
# df_mod_mean$new.slope<-c(-0.02807+df_mod_mean$Intercepts)
# df_mod_mean$min<-c(df_mod_mean$Slopes-df_mod_mean$sd.slopes)
# df_mod_mean$max<-c(df_mod_mean$Slopes+df_mod_mean$sd.slopes)
# df_mod_mean[df_mod_mean$Slopes>0&df_mod_mean$min>0|df_mod_mean$Slopes<0&df_mod_mean$max<0,"new.names"] #34 significative different slopes
# df_LDMC$direction[c(df_LDMC$Slopes>0&df_LDMC$min>0|df_LDMC$Slopes<0&df_LDMC$max<0 ) & df_LDMC$new.slope<0]<-"negative"
# df_LDMC$direction[c(df_LDMC$Slopes>0&df_LDMC$min>0|df_LDMC$Slopes<0&df_LDMC$max<0 ) & df_LDMC$new.slope>0]<-"positive"
# df_LDMC$direction[is.na(df_LDMC$direction)]<-"none"
# 
# LDMC_positive_slopes<-df_LDMC[c(df_LDMC$Slopes>0&df_LDMC$min>0|df_LDMC$Slopes<0&df_LDMC$max<0 ) & df_LDMC$new.slope>0,]
# write.table(LDMC_positive_slopes, "positive_slopes.txt")

names(df_mod_mean)
df_mod_mean$new.names
summary(mod_mean_z.cv) # Intercept -0.02807  std error 0.03609
#random Intercepts
p_intercept <- ggplot(df_mod_mean,aes(new.names,Intercepts))
#Added horizontal line at y=0, error bars to points and points with size two
p_intercept <- p_intercept  +geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0) + geom_point(size=2)
p_intercept <- p_intercept + geom_hline(yintercept=-0.02807, col='red') +  geom_hline(yintercept=c(-0.02807-0.03609), col='red', alpha=0.5) +  geom_hline(yintercept=c(-0.02807+0.03609), col='red', alpha=0.5)
p_intercept <- p_intercept + coord_flip()
#Removed legends 
p_intercept <- p_intercept + theme_bw() + theme(legend.position="none") #axis.text.x = element_text(hjust=1, vjust=1, angle = 45)
#Changed x and y axis labels
p_intercept <- p_intercept + labs( x="Dataset",  y="Random Intercept Effect")
#change colors
# p_intercept <- p_intercept + scale_color_manual(values=rep("black", 77)) 


# ggarrange(p_intercept, p_slopes, ncol=2, align = "h",widths =  c(2, 1.6))

#34  datasets with significatively different slope (out of 77) 44% 
#16 resulted in a positive slope 21%


#try with function from https://stackoverflow.com/questions/34344599/a-caterpillar-plot-of-just-the-significant-random-effects-from-a-mixed-effects
# ggCaterpillar(ranef(CV_LDMC, postVar = TRUE)[[3]], QQ=FALSE, likeDotplot=TRUE, reorder=FALSE) #not working

# pdf(file = "C:/Users/Luisi/Dropbox/Tesi_Luisa/PhD/Tübingen/Paper/Journal of Ecology/Re-submission/figures/caterpillar_reproduction_nozero.pdf", 
#     useDingbats = FALSE, width = 8.3, height =  6 )
# print(p)
# dev.off()


#<---------------->#########
#DETRENDED CV T3--------------------------------------
#Final model but using CV t3 - sp_cvt3
summary(mod_mean_cv)

mod_mean_cv@call$formula
names(mean_stability)

#zscores
mean_stability$z.sp_cv_t3<-
  unlist(tapply(mean_stability$'sp_cv_t3', mean_stability$'Group.2', function(x){scale.default(na.omit(x))}))
names(mean_stability)

mod_mean_cvt3<- lmer(z.sp_cv_t3 ~  z.mean_LeafN +    
                              z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC +
                              (1|Group.1) + (1|Group.2), 
                            data= mean_stability%>%filter(woodyness_new=='non-woody'), REML = F)

summary(mod_mean_cvt3)
c.mean_z.cvt3<-coefplot(mod_mean_cvt3, intercept = F, title='', color="black",
                      pointSize = 5, lwdInner = 3, lwdOuter = 1.5, sort = 'natural')
c.mean_z.cvt3<-c.mean_z.cvt3+  labs(y="", x="Coefficent") + theme_classic() +
  scale_y_discrete(labels=c("z.mean_LeafN"="Leaf N content", 
                            "z.log.mean_SeedMass"="Seed Mass",
                            "z.log.mean_SLA"="SLA",
                            "z.mean_LDMC"="LDMC"))
r.squaredGLMM(mod_mean_cvt3)


#<---------------->#########
#SINGLE TRAIT MODELS-----------------------------

# z.sp_cv ~ z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
#   z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD +
#   (1|Group.1) + (1|Group.2), 
# data= mean_stability%>%filter(woodyness_new=='non-woody')
#Intercorrelations between variables in the additive model (even though vif values are ok)

#z.scores for CV values------
boxplot(mean_stability$z.sp_cv~mean_stability$Group.2)
boxplot(mean_stability$sp_cv~mean_stability$Group.2)
boxplot(mean_stability$sp_cv_t3~mean_stability$Group.2)

boxplot(mean_stability$z.sp_cv_t3~mean_stability$Group.2)
boxplot(mean_stability$z.sp_cvt3~mean_stability$Group.2)

#what is the difference between doing this and what I did in line 435??? - this seems to better align mean values between datasets
mean_stability<-mean_stability%>%group_by(Group.2) %>%
            mutate(z.sp_cvt3 = scale.default(na.omit(sp_cv_t3)))

boxplot(mean_stability$z.sp_cvt3~mean_stability$Group.2)

hist(mean_stability$z.sp_cvt3)

length(mean_stability$z.sp_cvt3)
summary(mean_stability$z.sp_cvt3)
sum(is.na(mean_stability$z.sp_cvt3))
mean_stability[which(is.na(mean_stability$z.sp_cvt3)),"Group.2"]
mean_stability[mean_stability$Group.2=="LepsOhrazeni", c("Group.1","Group.2", "sp_cv_t3", "z.sp_cvt3", "z.sp_cv_t3"), ] 


prova<-mean_stability[mean_stability$Group.2=="LepsOhrazeni", ]

prova$prova<-c(prova$sp_cv_t3-mean(prova$sp_cv_t3))/sd(prova$sp_cv_t3)

prova$z.sp_cv_t3
prova$z.sp_cvt3

prova$prova%in%prova$z.sp_cv_t3
prova$prova%in%prova$z.sp_cvt3

leps<-tapply(mean_stability$'sp_cv_t3', mean_stability$'Group.2', function(x){scale.default(na.omit(x))})["LepsOhrazeni"]

prova$z.sp_cv_t3%in%as.vector(unlist(leps))



#CV and CVt3
names(mean_stability)
lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(sp_cv), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.sp_cv), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(sp_cv_t3), hist)



#Check correlations--------------

#using spearman correlations and original variables
names(mean_stability)
corpat <- as.data.frame(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.log.mean_height, z.mean_LeafN, z.mean_LeafP, z.log.mean_SeedMass, z.log.mean_SLA, z.mean_LDMC, z.mean_SSD))
dim(corpat)
## Spearman rs Correlations
print(corr.p(cor(corpat, method = c("spearman"), use='pairwise.complete.obs'),n=3397,adjust="holm"),short=FALSE)
M1 <-cor(corpat, method="spearman", use='pairwise.complete.obs')

# matrix of the p-value of the correlation
cor1 <- corr.p(cor(corpat, method = c("spearman"), use='pairwise.complete.obs'),3397,adjust="holm") 
p.mat1 <- cor1$p # Probability values for plotting
head(p.mat1[, 1:6])


# Plot

corrplot(M1, method = "number", type="upper",cl.pos="r", tl.pos="d", tl.col="black", 
         tl.srt=60, cl.cex= 0.75, number.cex= 0.75, tl.cex = 0.75, 
         col=brewer.pal(n=8, name="RdYlBu"), p.mat = p.mat1, sig.level = 0.05, 
         insig = "blank" )

# method = c("circle", "square", "ellipse", "number", "shade", "color", "pie"),
# type = c("full", "lower", "upper")

names(mean_stability)

#do it by hand
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(corpat, lower.panel = panel.smooth, upper.panel=panel.cor)


scatter.smooth()

#Height-----


#Outliers in Height?? - 10 points (7 spp) with height > 2 m - REMOVE?
lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(mean_height), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_height>2)%>%select(mean_height), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_height<2)%>%select(mean_height), hist)

mean_stability%>%filter(woodyness_new=='non-woody', mean_height>2)%>%select(Group.1, Group.2, mean_height)
mean_stability%>%filter(woodyness_new=='non-woody', mean_height<2)%>%select(Group.1, Group.2, mean_height)
mean_stability%>%filter(woodyness_new=='non-woody', mean_height<2, mean_height>1)%>%select(Group.1, Group.2, mean_height)

# transform with log
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_height<2)%>%select(mean_height), function(x){hist(log(x))})
mean_stability<-mutate(mean_stability, log_mean_height=log(mean_height) )
lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(log_mean_height), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(z.log.mean_height), hist)
#centered and scaled is better?

mean_stability<-mean_stability%>%ungroup()

# Check original values on original script
# head(Traits)
# hist(Traits$mean_Height[Traits$mean_Height<10&Traits$mean_Height>2])
# Traits[Traits$mean_Height>2,c('Species', 'mean_Height')]


mod_mean_cv.height<- lmer(z.sp_cv ~ 
                            z.log.mean_height + 
                                (1|Group.1) + 
                                (1|Group.2) , 
                              data= mean_stability%>%filter(woodyness_new=='non-woody', mean_height<2), REML = F)
# boundary (singular) fit: see ?isSingular - when using random slope
# Number of obs: 2809, groups:  Group.1, 1111; Group.2, 78
summary(mod_mean_cv.height) #AIC 7844.553 (random slope) #7841.7 (just random intercept)
coefplot(mod_mean_cv.height) #nope
r.squaredGLMM(mod_mean_cv.height) #2.044468e-05 0.1810405 (just random intercept)
visreg(mod_mean_cv.height)
hist(ranef(mod_mean_cv.height)$Group.2$`(Intercept)`)

mod_mean_cvt3.height<- lmer(z.sp_cvt3 ~ 
                              z.log.mean_height + 
                            (1|Group.1) , 
                          data= mean_stability%>%filter(woodyness_new=='non-woody', mean_height<2), REML = F)
#singular fit when using not-scaled trait
plot(mod_mean_cvt3.height)
hist(resid(mod_mean_cvt3.height))
summary(mod_mean_cvt3.height) #AIC 7874.4 (random slope) #  7870.7 (just random intercept)
coefplot(mod_mean_cvt3.height, intercept = F) #nope
r.squaredGLMM(mod_mean_cvt3.height)# 0.0009103496 0.174432 (random slope) # 0.001065647 0.1728588  (just random intercept)
visreg(mod_mean_cvt3.height) #, xtrans = log

#Random effects
hist(ranef(mod_mean_cvt3.height)$Group.2$`(Intercept)`)
hist(ranef(mod_mean_cvt3.height)$Group.1$`(Intercept)`)

row.names(ranef(mod_mean_cvt3.height)$Group.2)

#check random slopes
hist(ranef(mod_mean_cvt3.height)$Group.2$z.log.mean_height )
max(ranef(mod_mean_cvt3.height)$Group.2$z.log.mean_height )


# plot random-slope-intercept
plot_model(mod_mean_cvt3.height,  type="pred",
           terms=c("z.log.mean_height","Group.2"),pred.type="re",show.legend=F )
# vars = "c12hour", sample.n = 15
#problem with color palette!!



#Leaf N-----------
names(mean_stability)

lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(mean_LeafN), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafN<60)%>%select(mean_LeafN), hist)


#Outliers Leaf N >60 1 spp REMOVED
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafN>60)%>%select(mean_LeafN), hist)
mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafN>60)%>%select(Group.1, mean_LeafN)

mod_mean_cv.leafN<- lmer(z.sp_cv ~ 
                            mean_LeafN + 
                            (1|Group.1) + 
                            (1|Group.2) , 
                          data= mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafN<60), REML = F)
#using random slope it gives a SINGULAR FIT
#Number of obs: 1938, groups:  Group.1, 550; Group.2, 77


setdiff(row.names(ranef(mod_mean_cv.height)$Group.2), row.names(ranef(mod_mean_cv.leafN)$Group.2))



summary(mod_mean_cv.leafN)
coefplot(mod_mean_cv.leafN, intercept = F) #significant!
r.squaredGLMM(mod_mean_cv.leafN) #0.01990624 0.1726765
visreg(mod_mean_cv.leafN) #positive
#check random slopes
hist(ranef(mod_mean_cv.leafN)$Group.2$mean_LeafN ) #one dataset is high
max(ranef(mod_mean_cv.leafN)$Group.2$mean_LeafN )


# plot random-slope-intercept
plot_model(mod_mean_cv.leafN,  type="pred",
           terms=c("mean_LeafN","Group.2"),pred.type="re",show.legend=F )
# vars = "c12hour", sample.n = 15
#problem with color palette??


#using t3
mod_mean_cvt3.leafN<- lmer(z.sp_cvt3 ~ 
                             z.mean_LeafN + 
                           (1|Group.1)  , 
                         data= mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafN<60), REML = F)
# Number of obs: 1938, groups:  Group.1, 550; Group.2, 77
summary(mod_mean_cvt3.leafN) #AIC 5377.1 (random slope) #5373.6 (only intercept)
coefplot(mod_mean_cvt3.leafN, intercept = F) #significant!
r.squaredGLMM(mod_mean_cvt3.leafN) #0.01630545 0.155356 (random slope) #0.01703351 0.1516496 (only intercept)
visreg(mod_mean_cvt3.leafN) #positive
#check random slopes
hist(ranef(mod_mean_cvt3.leafN)$Group.2$z.mean_LeafN ) 
max(ranef(mod_mean_cvt3.leafN)$Group.2$z.mean_LeafN )


# plot random-slope-intercept
plot_model(mod_mean_cvt3.leafN,  type="pred",
           terms=c("z.mean_LeafN","Group.2"),pred.type="re",show.legend=F )
# vars = "c12hour", sample.n = 15
#problem with color palette??


#Leaf P-----------
names(mean_stability)

lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(mean_LeafP), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafP<6)%>%select(mean_LeafN), hist)


#Outliers Leaf P > 6 10 points (3 spp) REMOVED 
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafP>6)%>%select(mean_LeafP), hist)
mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafP>6)%>%select(Group.1, mean_LeafP)

mod_mean_cv.leafP<- lmer(z.sp_cv ~ 
                          mean_LeafP + 
                          (1|Group.1) + 
                          (mean_LeafP|Group.2) , 
                        data= mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafP<6), REML = F)
#using random slope it gives a SINGULAR FIT


summary(mod_mean_cv.leafP)
coefplot(mod_mean_cv.leafP, intercept = F) #significant!
r.squaredGLMM(mod_mean_cv.leafP) # 0.01006328 0.1424406
visreg(mod_mean_cv.leafP) #positive
#check random slopes
hist(ranef(mod_mean_cv.leafP)$Group.2$mean_LeafP ) #one dataset is high
max(ranef(mod_mean_cv.leafP)$Group.2$mean_LeafP ) #Rose_ECN

# plot random-slope-intercept - adjust to show just the extremes
plot_model(mod_mean_cv.leafP,  type="pred",
           terms=c("mean_LeafP","Group.2"),pred.type="re",show.legend=F , show.data = T) #
# vars = "c12hour", sample.n = 15
#problem with color palette??

#trying with t3
mod_mean_cvt3.leafP<- lmer(z.sp_cvt3 ~ 
                             z.mean_LeafP + 
                           (1|Group.1)  , 
                         data= mean_stability%>%filter(woodyness_new=='non-woody', mean_LeafP<6), REML = F)
# Number of obs: 1640, groups:  Group.1, 415; Group.2, 76
summary(mod_mean_cvt3.leafP) #AIC 4562.9 #
coefplot(mod_mean_cvt3.leafP, intercept = F) #nope!
r.squaredGLMM(mod_mean_cvt3.leafP) #0.002539381 0.1529264 (random slope) #
visreg(mod_mean_cvt3.leafP) #positive
#check random slopes
hist(ranef(mod_mean_cvt3.leafP)$Group.2$z.mean_LeafP ) 

# plot random-slope-intercept - adjust to show just the extremes
plot_model(mod_mean_cvt3.leafP,  type="pred",
           terms=c("z.mean_LeafP","Group.2"),pred.type="re",show.legend=F , show.data = T) #
# vars = "c12hour", sample.n = 15
#problem with color palette??

#Seed Mass-----------
names(mean_stability)

lapply(ungroup(mean_stability)%>%filter(woodyness_new=='non-woody')%>%select(mean_SeedMass), hist)
lapply(ungroup(mean_stability)%>%filter(woodyness_new=='non-woody', mean_SeedMass<100)%>%select(mean_SeedMass), function(x){hist(log(x))})


#Outliers Seed Mass >100 9 points (6 spp) are outliers
lapply(ungroup(mean_stability)%>%filter(woodyness_new=='non-woody', mean_SeedMass>100)%>%select(mean_SeedMass), hist)
mean_stability%>%filter(woodyness_new=='non-woody', mean_SeedMass>100)%>%select(Group.1, mean_SeedMass)

lapply(ungroup(mean_stability)%>%filter(woodyness_new=='non-woody')%>%select(z.log.mean_SeedMass), function(x){hist(x)})


mod_mean_cv.seedmass<- lmer(z.sp_cv ~ 
                              z.log.mean_SeedMass + 
                           (1|Group.1) + 
                           ( z.log.mean_SeedMass |Group.2) , 
                         data= mean_stability%>%filter(woodyness_new=='non-woody', z.log.mean_SeedMass>-2), REML = F)
#using random slope it gives a SINGULAR FIT - fixed when using z.log.mean_SeedMass>-2
#Number of obs: 2979, groups:  Group.1, 1211; Group.2, 78


summary(mod_mean_cv.seedmass)
coefplot(mod_mean_cv.seedmass) #significant!
r.squaredGLMM(mod_mean_cv.seedmass) #0.003223209 0.1950927
visreg(mod_mean_cv.seedmass) #negative
#check random slopes
hist(ranef(mod_mean_cv.seedmass)$Group.2$mean_SeedMass ) #one dataset is high
max(ranef(mod_mean_cv.seedmass)$Group.2$mean_SeedMass ) #Rose_ECN


# plot random-slope-intercept
plot_model(mod_mean_cv.seedmass,  type="pred",
           terms=c("z.log.mean_SeedMass","Group.2"),pred.type="re",show.legend=F)
# vars = "c12hour", sample.n = 15
#problem with color palette??

#Caterpillar plot for random slope 
summary(mod_mean_cv.seedmass)
plot_model(mod_mean_cv.seedmass, type='re') 

#do it by hand
ranef(mod_mean_cv.seedmass)$Group.1 #species (651)
ranef(mod_mean_cv.seedmass)$Group.2 #Datasets (77) 

#in the model there are 2 dataset less than in the original data
randoms_mod_mean.sm<-ranef(mod_mean_cv.seedmass, condVar = TRUE)$Group.2
#condVar - an optional logical argument indicating if the conditional variance-covariance matrices 
#of the random effects should be added as an attribute.
#check Dataset names
randoms_mod_mean.sm[is.na(match(row.names(randoms_mod_mean.sm), ids_datasets$Original_dataset)),]
which(is.na(match(row.names(randoms_mod_mean.sm), ids_datasets$Original_dataset)))
row.names(randoms_mod_mean.sm)[51]<-"Penuelas_Garraf_Penuelas" #Peñuelas not working on PC
ids_datasets$Original_dataset[48]<-"Penuelas_Garraf_Penuelas"
match(row.names(randoms_mod_mean.sm), ids_datasets$Original_dataset)
ids_datasets[is.na(match(ids_datasets$Original_dataset, row.names(randoms_mod_mean.sm))),]
randoms_mod_mean.sm$ids<-ids_datasets$`Dataset ID`[na.omit(match(ids_datasets$Original_dataset, row.names(randoms_mod_mean.sm)))]

qq_mod_mean.sm <- attr(ranef(mod_mean_cv.seedmass, condVar = TRUE)[[2]], "postVar")
dim(qq_mod_mean.sm)
qq_mod_mean.sm[1,1,1:78] #for random intercept
qq_mod_mean.sm[2,2,1:78] #for random slope
# rand.interc<-randoms
df_mod_mean.sm<-data.frame(lev.names=randoms_mod_mean.sm$ids, 
                        Intercepts=randoms_mod_mean.sm[,1],
                        sd.interc= 2*sqrt(qq_mod_mean.sm[1,1,1:78]),
                        Slopes=randoms_mod_mean.sm[,2],
                        sd.slopes= 2*sqrt(qq_mod_mean.sm[2,2,1:78]))
# df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)]) #order intercepts when plotting
# df$lev.names<-factor(df$lev.names, levels = rev(levels(df$lev.names)))

#adjust slope values
summary(mod_mean_cv.seedmass) # Slope -0.07691    std error  0.03148
df_mod_mean.sm$new.slope<-c(-0.07691+df_mod_mean.sm$Slopes)
df_mod_mean.sm$min<-c(df_mod_mean.sm$new.slope-df_mod_mean.sm$sd.slopes)
df_mod_mean.sm$max<-c(df_mod_mean.sm$new.slope+df_mod_mean.sm$sd.slopes)
# df_mod_mean[df_mod_mean$Slopes>0&df_mod_mean$min>0|df_mod_mean$Slopes<0&df_mod_mean$max<0,"lev.names"] #34 significative different slopes
# df_mod_mean$direction[c(df_mod_mean$Slopes>0&df_mod_mean$min>0|df_mod_mean$Slopes<0&df_mod_mean$max<0 ) & df_mod_mean$new.slope<0]<-"negative"
# df_mod_mean$direction[c(df_mod_mean$Slopes>0&df_mod_mean$min>0|df_mod_mean$Slopes<0&df_mod_mean$max<0 ) & df_mod_mean$new.slope>0]<-"positive"
# df_mod_mean$direction[is.na(df_mod_mean$direction)]<-"none"
# df_mod_mean$lev.names<-factor(df_mod_mean$lev.names,levels=df_mod_mean$lev.names[order(df_mod_mean$new.slope)]) #order when plotting
# LDMC_positive_slopes<-df_LDMC[c(df_LDMC$Slopes>0&df_LDMC$min>0|df_LDMC$Slopes<0&df_LDMC$max<0 ) & df_LDMC$new.slope>0,]
# write.table(LDMC_positive_slopes, "positive_slopes.txt")

names(df_mod_mean.sm)
df_mod_mean.sm$"lev.names" 
df_mod_mean.sm$new.slope
df_mod_mean.sm<-as_tibble(df_mod_mean.sm)
df_mod_mean.sm<-arrange(df_mod_mean.sm, by=new.slope)
df_mod_mean.sm<- mutate(df_mod_mean.sm,lev.names=factor(lev.names, levels=lev.names))  
summary(mod_mean_cv.seedmass) # Slope -0.23309    std error  0.03184
#random slopes
p_slopes_SM <- ggplot(df_mod_mean.sm,aes(lev.names,new.slope))
#Added horizontal line at y=0, error bars to points and points with size two
p_slopes_SM <- p_slopes_SM + geom_point(size=2) + geom_errorbar(aes(ymin=min, ymax=max), width=0) 
p_slopes_SM <- p_slopes_SM + geom_hline(yintercept=-0.07691, col='red') +  geom_hline(yintercept=c(-0.07691- 0.03148), col='red', alpha=0.5) +  geom_hline(yintercept=c(-0.07691+ 0.03148), col='red', alpha=0.5)
p_slopes_SM <- p_slopes_SM + geom_hline(yintercept=0, col='black',lty=2)
p_slopes_SM <- p_slopes_SM + coord_flip()
#Removed legends 
p_slopes_SM <- p_slopes_SM + theme_bw() + theme(legend.position="none") #axis.text.x = element_text(hjust=1, vjust=1, angle = 45)
p_slopes_SM <- p_slopes_SM + scale_x_discrete(guide = guide_axis(n.dodge=2))
#Changed x and y axis labels
p_slopes_SM <- p_slopes_SM + labs( x="Dataset",  y="Random Effect Slopes (seed mass)")
#change colors
# p_slopes_LDMC <- p_slopes_LDMC+ scale_color_manual(values=c("none"="black", "negative"="red", "positive"="green3")) 


# ggarrange(p_intercept, p_slopes, ncol=2, align = "h",widths =  c(2, 1.6))

#34  datasets with significatively different slope (out of 77) 44% 
#16 resulted in a positive slope 21%

#try with t3
mod_mean_cvt3.seedmass<- lmer(z.sp_cvt3 ~ 
                              z.log.mean_SeedMass + 
                              (1|Group.1) , 
                            data= mean_stability%>%filter(woodyness_new=='non-woody', mean_SeedMass<100), REML = F)

#Number of obs: 2979, groups:  Group.1, 1211; Group.2, 78

summary(mod_mean_cvt3.seedmass) #AIC 8296.8 (random slope)# AIC 8293.7 (just random intercept)
coefplot(mod_mean_cvt3.seedmass) #significant!
r.squaredGLMM(mod_mean_cvt3.seedmass) #0.007383915 0.1950578 (random slope) # 0.006797932 0.1905123 (just random intercept)
visreg(mod_mean_cvt3.seedmass) #negative
#check random slopes
hist(ranef(mod_mean_cvt3.seedmass)$Group.2$z.log.mean_SeedMass ) 

# plot random-slope-intercept
plot_model(mod_mean_cvt3.seedmass,  type="pred",
           terms=c("z.log.mean_SeedMass","Group.2"),pred.type="re",show.legend=F)
# vars = "c12hour", sample.n = 15
#problem with color palette??


#SLA-----------
names(mean_stability)


lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(mean_SLA), hist)
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_SLA<75)%>%select(mean_SLA), hist)

#Outliers Seed Mass >75 6 points (3 spp) REMOVED
lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_SLA>75)%>%select(mean_SLA), hist)
mean_stability%>%filter(woodyness_new=='non-woody', mean_SLA>75)%>%select(Group.1, mean_SLA)

lapply(mean_stability%>%filter(woodyness_new=='non-woody', mean_SLA<75)%>%select(mean_SLA), function(x){hist(log(x))})


mod_mean_cv.sla<- lmer(z.sp_cv ~ 
                         z.log.mean_SLA + 
                              (1|Group.1) + 
                              (1|Group.2) , 
                            data= mean_stability%>%filter(woodyness_new=='non-woody', mean_SLA<75), REML = F)
#using random slope Model failed to converge - transform variable before fitting?

summary(mod_mean_cv.sla)
coefplot(mod_mean_cv.sla) #significant!
r.squaredGLMM(mod_mean_cv.sla) #0.01743962 0.1789547
visreg(mod_mean_cv.sla) #positive
#check random slopes
hist(ranef(mod_mean_cv.sla)$Group.2$`(Intercept)` ) #one dataset is low

# plot random effects
plot_model(mod_mean_cv.sla,  type="pred",
           terms=c("Group.2"),pred.type="re",show.legend=F)
# vars = "c12hour", sample.n = 15
# problem with color palette??

#try with ct3
mod_mean_cvt3.sla<- lmer(z.sp_cvt3 ~ 
                           z.log.mean_SLA + 
                         (1|Group.1),
                       data= mean_stability%>%filter(woodyness_new=='non-woody', mean_SLA<75), REML = F)
# boundary (singular) fit: see ?isSingular - modello con random slope e senza slope
# Number of obs: 2467, groups:  Group.1, 870; Group.2, 76
summary(mod_mean_cvt3.sla) #AIC  6875.9 (with both spp and dataset random intercept) #AIC   6873.9 (with only spp)
coefplot(mod_mean_cvt3.sla) #significant!
r.squaredGLMM(mod_mean_cvt3.sla) #0.01423919 0.1696807 #R squared stays the same without dataset random effect
visreg(mod_mean_cvt3.sla) #positive


#LDMC-----------
names(mean_stability)


lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(mean_LDMC), hist)


mod_mean_cv.ldmc<- lmer(z.sp_cv ~ 
                          z.mean_LDMC + 
                         (1|Group.1) + 
                         (z.mean_LDMC|Group.2) , 
                       data= mean_stability%>%filter(woodyness_new=='non-woody'), REML = F)
#Number of obs: 2090, groups:  Group.1, 651; Group.2, 77

summary(mod_mean_cv.ldmc)
coefplot(mod_mean_cv.ldmc) #significant!
r.squaredGLMM(mod_mean_cv.ldmc) #0.04696485 0.1686547
visreg(mod_mean_cv.ldmc) #negative
#check random slopes
hist(ranef(mod_mean_cv.ldmc)$Group.2$z.mean_LDMC ) #no outliers


# plot random-slope-intercept
plot_model(mod_mean_cv.ldmc,  type="pred",
           terms=c("z.mean_LDMC","Group.2"),pred.type="re",show.legend=F, show.data = T)
# vars = "c12hour", sample.n = 15
#problem with color palette??

#Caterpillar plot for random slope 
summary(mod_mean_cv.ldmc)
plot_model(mod_mean_cv.ldmc, type='re') 

#do it by hand
ranef(mod_mean_cv.ldmc)$Group.1 #species (651)
ranef(mod_mean_cv.ldmc)$Group.2 #Datasets (77) 

#in the model there are 2 dataset less than in the original data
randoms_mod_mean<-ranef(mod_mean_cv.ldmc, condVar = TRUE)$Group.2
#condVar - an optional logical argument indicating if the conditional variance-covariance matrices 
#of the random effects should be added as an attribute.
#check Dataset names
randoms_mod_mean[is.na(match(row.names(randoms_mod_mean), ids_datasets$Original_dataset)),]
which(is.na(match(row.names(randoms_mod_mean), ids_datasets$Original_dataset)))
row.names(randoms_mod_mean)[51]<-"Penuelas_Garraf_Peñuelas" #not working on PC
ids_datasets$Original_dataset[48]<-"Penuelas_Garraf_Penuelas"
match(row.names(randoms_mod_mean), ids_datasets$Original_dataset)
ids_datasets[is.na(match(ids_datasets$Original_dataset, row.names(randoms_mod_mean))),]
randoms_mod_mean$ids<-ids_datasets$`Dataset ID`[na.omit(match(ids_datasets$Original_dataset, row.names(randoms_mod_mean)))]

qq_mod_mean <- attr(ranef(mod_mean_cv.ldmc, condVar = TRUE)[[2]], "postVar")
qq_mod_mean[1,1,1:77] #for random intercept
qq_mod_mean[2,2,1:77] #for random slope
# rand.interc<-randoms
df_mod_mean<-data.frame(lev.names=randoms_mod_mean$ids, 
                        Intercepts=randoms_mod_mean[,1],
                        sd.interc= 2*sqrt(qq_mod_mean[1,1,1:77]),
                        Slopes=randoms_mod_mean[,2],
                        sd.slopes= 2*sqrt(qq_mod_mean[2,2,1:77]))
# df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)]) #order intercepts when plotting
# df$lev.names<-factor(df$lev.names, levels = rev(levels(df$lev.names)))

#adjust slope values
summary(mod_mean_cv.ldmc) # Slope -0.23309    std error  0.03184
df_mod_mean$new.slope<-c(-0.23309+df_mod_mean$Slopes)
df_mod_mean$min<-c(df_mod_mean$new.slope-df_mod_mean$sd.slopes)
df_mod_mean$max<-c(df_mod_mean$new.slope+df_mod_mean$sd.slopes)
# df_mod_mean[df_mod_mean$Slopes>0&df_mod_mean$min>0|df_mod_mean$Slopes<0&df_mod_mean$max<0,"lev.names"] #34 significative different slopes
# df_mod_mean$direction[c(df_mod_mean$Slopes>0&df_mod_mean$min>0|df_mod_mean$Slopes<0&df_mod_mean$max<0 ) & df_mod_mean$new.slope<0]<-"negative"
# df_mod_mean$direction[c(df_mod_mean$Slopes>0&df_mod_mean$min>0|df_mod_mean$Slopes<0&df_mod_mean$max<0 ) & df_mod_mean$new.slope>0]<-"positive"
# df_mod_mean$direction[is.na(df_mod_mean$direction)]<-"none"
# df_mod_mean$lev.names<-factor(df_mod_mean$lev.names,levels=df_mod_mean$lev.names[order(df_mod_mean$new.slope)]) #order when plotting
# LDMC_positive_slopes<-df_LDMC[c(df_LDMC$Slopes>0&df_LDMC$min>0|df_LDMC$Slopes<0&df_LDMC$max<0 ) & df_LDMC$new.slope>0,]
# write.table(LDMC_positive_slopes, "positive_slopes.txt")

names(df_mod_mean)
df_mod_mean$"lev.names" 
df_mod_mean$new.slope
df_mod_mean<-as_tibble(df_mod_mean)
df_mod_mean<-arrange(df_mod_mean, by=new.slope)
df_mod_mean<- mutate(df_mod_mean,lev.names=factor(lev.names, levels=lev.names))  
summary(mod_mean_cv.ldmc) # Slope -0.23309    std error  0.03184
#random slopes
p_slopes_LDMC <- ggplot(df_mod_mean,aes(lev.names,new.slope))
#Added horizontal line at y=0, error bars to points and points with size two
p_slopes_LDMC <- p_slopes_LDMC + geom_point(size=2) + geom_errorbar(aes(ymin=min, ymax=max), width=0) 
p_slopes_LDMC <- p_slopes_LDMC + geom_hline(yintercept=-0.23309, col='red') +  geom_hline(yintercept=c(-0.23309- 0.03184), col='red', alpha=0.5) +  geom_hline(yintercept=c(-0.23309+ 0.02958), col='red', alpha=0.5)
p_slopes_LDMC <- p_slopes_LDMC + geom_hline(yintercept=0, col='black',lty=2)
p_slopes_LDMC <- p_slopes_LDMC + coord_flip()
#Removed legends 
p_slopes_LDMC <- p_slopes_LDMC + theme_bw() + theme(legend.position="none") #axis.text.x = element_text(hjust=1, vjust=1, angle = 45)
p_slopes_LDMC <- p_slopes_LDMC + scale_x_discrete(guide = guide_axis(n.dodge=2))
#Changed x and y axis labels
p_slopes_LDMC <- p_slopes_LDMC + labs( x="Dataset",  y="Random Effect Slopes (LDMC)")
#change colors
# p_slopes_LDMC <- p_slopes_LDMC+ scale_color_manual(values=c("none"="black", "negative"="red", "positive"="green3")) 


# ggarrange(p_intercept, p_slopes, ncol=2, align = "h",widths =  c(2, 1.6))

#34  datasets with significatively different slope (out of 77) 44% 
#16 resulted in a positive slope 21%

# try t3
mod_mean_cvt3.ldmc<- lmer(sp_cv_t3 ~ 
                            z.mean_LDMC + 
                          (1|Group.1) + (z.mean_LDMC|Group.2), 
                        data= mean_stability%>%filter(woodyness_new=='non-woody'), REML = F)
#Number of obs: 2090, groups:  Group.1, 651; Group.2, 77

summary(mod_mean_cvt3.ldmc) #AIC  5763.5 (random slope) #5764.1 (only random intercept)
coefplot(mod_mean_cvt3.ldmc) #significant!
r.squaredGLMM(mod_mean_cvt3.ldmc) #0.03667784 0.1747909 (random slope) # 0.03389146 0.1582597 (only random intercept)
visreg(mod_mean_cvt3.ldmc) #negative
#check random slopes
hist(ranef(mod_mean_cvt3.ldmc)$Group.2$z.mean_LDMC ) #no outliers


# plot random-slope-intercept
plot_model(mod_mean_cvt3.ldmc,  type="pred",
           terms=c("z.mean_LDMC","Group.2"),pred.type="re",show.legend=F, show.data = T)
# vars = "c12hour", sample.n = 15
#problem with color palette??


#put all models together----------
tab<-export_summs(mod_mean_cvt3.height, mod_mean_cvt3.seedmass, 
                  mod_mean_cvt3.sla, mod_mean_cvt3.leafN, mod_mean_cvt3.leafP, mod_mean_cvt3.ldmc,
                  note=NULL, model.names = c("Height", "Seed mass", 
                                             "SLA", "Leaf N content", 
                                             "Leaf P content", "LDMC") ) # results = 'asis'

quick_xlsx(tab, file='tab.xlsx')
# tab$names[1]<-"GEN_MPD_Multi_trait"
number_format(tab)<-gsub("%.2f", "%.4f", number_format(tab))
# tab_final<-tab
# 
# tab_final<- add_rows(tab_final, tab)

#caterpillar plots-----------
ggarrange(p_slopes_LDMC, p_slopes_SM,labels = 'auto')



#coefficent plots-------
p<-plot_summs(
  mod_mean_cvt3.height, mod_mean_cvt3.seedmass, 
  mod_mean_cvt3.sla, mod_mean_cvt3.leafN, 
  mod_mean_cvt3.leafP, mod_mean_cvt3.ldmc,
  point.shape=F, colors = rep("black", 6),
  coefs = c( 'Height'='z.log.mean_height', 
            'Seed Mass'='z.log.mean_SeedMass' ,
             'SLA'='z.log.mean_SLA' ,
             'Leaf N content'='z.mean_LeafN',
             'Leaf P content'='z.mean_LeafP',
             'LDMC'='z.mean_LDMC')) 

p + theme(legend.position = "none") #+ ggtitle('Models ')


#add r.squared values?
r.squared.m<- c(r.squaredGLMM(mod_mean_cvt3.height)[1], r.squaredGLMM(mod_mean_cvt3.seedmass)[1],
r.squaredGLMM(mod_mean_cvt3.sla)[1], r.squaredGLMM(mod_mean_cvt3.leafN)[1],
r.squaredGLMM(mod_mean_cvt3.leafP)[1], r.squaredGLMM(mod_mean_cvt3.ldmc)[1])

r.squared.c<- c(r.squaredGLMM(mod_mean_cvt3.height)[2], r.squaredGLMM(mod_mean_cvt3.seedmass)[2],
                r.squaredGLMM(mod_mean_cvt3.sla)[2], r.squaredGLMM(mod_mean_cvt3.leafN)[2],
                r.squaredGLMM(mod_mean_cvt3.leafP)[2], r.squaredGLMM(mod_mean_cvt3.ldmc)[2])


str(summary(mod_mean_cv.height))

# p+ geom_text(y=, label=paste('R2=', round(p$data$r.squared, 3), sep=""))
# c(3.25,2.25,1.25,3.10,2.10,1.10,2.93,1.93,0.93)

#effect plots---------
#  mod_mean_cvt3.seedmass, 
# mod_mean_cvt3.sla, mod_mean_cvt3.leafN, 
# mod_mean_cvt3.ldmc,

#seed mass
mod_dat_sm<- get_model_data(mod_mean_cvt3.seedmass, type = "pred", terms = c("z.log.mean_SeedMass"), pretty = FALSE)

str(mod_mean_cvt3.seedmass)
dat_sm<-as_tibble(mod_mean_cvt3.seedmass@frame)
# dat_sm$meansppcv<-tapply(dat_sm$z.sp_cv_t3, dat_sm$Group.1, mean)
dat_sm<-mutate(dat_sm%>%group_by(Group.1), mean_cv=mean(z.sp_cvt3))

p1<-ggplot(aes(y=z.sp_cvt3, x=z.log.mean_SeedMass), data = dat_sm)
p1<-p1+  
  geom_point(col=gray(0.5, 0.25)) + #geom_abline(slope = 1.3697, intercept = -4.9615, col="red", lwd=1.5) 
  geom_point(inherit.aes = F, aes(y=mean_cv, x=z.log.mean_SeedMass), col='black') + 
  geom_ribbon(data = mod_dat_sm, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_sm, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, col= "red")+
  theme_classic() + labs(y="CV", x="Seed Mass")
# p1<- p1 + c+ annotate("text", x= 0.55, y = -13 ,  label = "paste(italic(R) ^ 2, \"m= 0.004\")", parse=T, size=6)+ 
  # annotate("text", x= 0.55, y = -14 ,  label = "paste(italic(R) ^ 2, \"c=0.602\")", parse=T, size=6)
  
#mod_mean_cvt3.sla
mod_dat_sla<- get_model_data(mod_mean_cvt3.sla, type = "pred", terms = c("z.log.mean_SLA"), pretty = FALSE)

str(mod_mean_cvt3.sla)
dat_sla<-tibble(mod_mean_cvt3.sla@frame)
# dat_sm$meansppcv<-tapply(dat_sm$z.sp_cv_t3, dat_sm$Group.1, mean)
dat_sla<-mutate(dat_sla%>%group_by(Group.1), mean_cv=mean(z.sp_cvt3))

p2<-ggplot(aes(y=z.sp_cvt3, x=z.log.mean_SLA), data = dat_sla)
p2<-p2+  
  geom_point(col=gray(0.5, 0.25)) + #geom_abline(slope = 1.3697, intercept = -4.9615, col="red", lwd=1.5) 
  geom_point(inherit.aes = F, aes(y=mean_cv, x=z.log.mean_SLA), col='black') + 
  geom_ribbon(data = mod_dat_sla, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_sla, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, col= "red")+
  theme_classic() + labs(y="CV", x="SLA")
# p1<- p1 + c+ annotate("text", x= 0.55, y = -13 ,  label = "paste(italic(R) ^ 2, \"m= 0.004\")", parse=T, size=6)+ 
# annotate("text", x= 0.55, y = -14 ,  label = "paste(italic(R) ^ 2, \"c=0.602\")", parse=T, size=6)

# mod_mean_cvt3.leafN
mod_dat_leafN<- get_model_data(mod_mean_cvt3.leafN, type = "pred", terms = c("z.mean_LeafN"), pretty = FALSE)

str(mod_mean_cvt3.sla)
dat_leafN<-tibble(mod_mean_cvt3.leafN@frame)
# dat_sm$meansppcv<-tapply(dat_sm$z.sp_cv_t3, dat_sm$Group.1, mean)
dat_leafN<-mutate(dat_leafN%>%group_by(Group.1), mean_cv=mean(z.sp_cvt3))

p3<-ggplot(aes(y=z.sp_cvt3, x=z.mean_LeafN), data = dat_leafN)
p3<-p3+  
  geom_point(col=gray(0.5, 0.25)) + #geom_abline(slope = 1.3697, intercept = -4.9615, col="red", lwd=1.5) 
  geom_point(inherit.aes = F, aes(y=mean_cv, x=z.mean_LeafN), col='black') + 
  geom_ribbon(data = mod_dat_leafN, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_leafN, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, col= "red")+
  theme_classic() + labs(y="CV", x="Leaf N content")

# mod_mean_cvt3.ldmc
mod_dat_ldmc<- get_model_data(mod_mean_cvt3.ldmc, type = "pred", terms = c("z.mean_LDMC"), pretty = FALSE)

dat_ldmc<-tibble(mod_mean_cvt3.ldmc@frame)
# dat_sm$meansppcv<-tapply(dat_sm$z.sp_cv_t3, dat_sm$Group.1, mean)
dat_ldmc<-mutate(dat_ldmc%>%group_by(Group.1), mean_cv=mean(z.sp_cvt3))

p4<-ggplot(aes(y=z.sp_cvt3, x=z.mean_LDMC), data = dat_ldmc)
p4<-p4+  
  geom_point(col=gray(0.5, 0.25)) + #geom_abline(slope = 1.3697, intercept = -4.9615, col="red", lwd=1.5) 
  geom_point(inherit.aes = F, aes(y=mean_cv, x=z.mean_LDMC), col='black') + 
  geom_ribbon(data = mod_dat_ldmc, aes(x = x,  ymin = conf.low, ymax=conf.high), 
              inherit.aes = F, show.legend = F, alpha = .25)+
  geom_line(data = mod_dat_ldmc, aes(x = x, y = predicted), 
            size = 1.025, show.legend = F, col= "red")+
  theme_classic() + labs(y="CV", x="LDMC")

#all plots together
ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, labels= 'auto')

visreg(mod_mean_cvt3.seedmass)

#++Qualitative traits---- to change into cvt3!!
mod_mean_cv.growthform<- lmer(z.sp_cv ~ 
                                growthform + 
                                (1|Group.1) + 
                                (1|Group.2) -1 , 
                              data= mean_stability, REML = F)
summary(mod_mean_cv.growthform)
visreg(mod_mean_cv.growthform)
coefplot(mod_mean_cv.growthform)


mod_mean_cv.woodyness<- lmer(z.sp_cv ~ woodyness_new + 
                               (1|Group.1) + (1|Group.2) -1 , 
                             data= mean_stability, REML = F)
coefplot(mod_mean_cv.woodyness)

mod_mean_cv.lifeform<- lmer(z.sp_cv ~ lifeform + 
                              (1|Group.1) + (1|Group.2) -1 , 
                            data= mean_stability, REML = F)
coefplot(mod_mean_cv.lifeform)

mod_mean_cv.lifespan<- lmer(z.sp_cv ~ lifespan_new + 
                              (1|Group.1) + (1|Group.2) -1 , 
                            data= mean_stability, REML = F)
coefplot(mod_mean_cv.lifespan)

#plots-------
est<-coefplot(mod_mean_cv.growthform)$data
# names(est)[2]<-"growthform"
est$Coefficient<- c("fern", "graminoid", "herb", "herb/shrub","shrub", "shrub/tree", "tree" )
m<- ggplot(data = mean_stability[!is.na(mean_stability$growthform)&!is.na(mean_stability$z.sp_cv ),], aes(y=z.sp_cv , x=growthform))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=5)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=2, inherit.aes = F, width = 0)+
  theme_classic(base_size=22) + xlab(label = "") +ylab(label = "species mean CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-0.7, 0.3), zoom.size = 1, zoom.data = z.sp_cv==0)
m

est<-coefplot(mod_mean_cv.woodyness)$data
# names(est)[2]<-"growthform"
est$Coefficient<- c("non-woody", "woody")
m<- ggplot(data = mean_stability[!is.na(mean_stability$woodyness_new)&!is.na(mean_stability$z.sp_cv ),], aes(y=z.sp_cv , x=woodyness_new))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=5)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=2, inherit.aes = F, width = 0)+
  theme_classic(base_size=22) + xlab(label = "") +ylab(label = "species mean CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-0.5, 0.1), zoom.size = 1, zoom.data = z.sp_cv==0)
m

est<-coefplot(mod_mean_cv.lifeform)$data
est$Coefficient<- c("Ch", "Cr", "H", "P", "T")
m<- ggplot(data = mean_stability[!is.na(mean_stability$lifeform)&!is.na(mean_stability$z.sp_cv ),], aes(y=z.sp_cv , x=lifeform))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=5)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=2, inherit.aes = F, width = 0)+
  theme_classic(base_size=22) + xlab(label = "") +ylab(label = "species mean CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-0.5, 0.8), zoom.size = 1, zoom.data = z.sp_cv==0)
m

est<-coefplot(mod_mean_cv.lifespan)$data
est$Coefficient<- c("annual", "not-annual")
m<- ggplot(data = mean_stability[!is.na(mean_stability$lifespan_new)&!is.na(mean_stability$z.sp_cv ),], aes(y=z.sp_cv , x=lifespan_new))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=5)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=2, inherit.aes = F, width = 0)+
  theme_classic(base_size=22) + xlab(label = "") +ylab(label = "species mean CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-0.4, 0.8), zoom.size = 1, zoom.data = z.sp_cv==0)
m

# Eventually fitting a single model for each dataset, 
#considering the traits emerging from the 'global' model 
#and their relationship with the average of species CV within each dataset. 
#Alternatively, present the caterpillar plot presenting the differences 
#between the intercepts of the datasets in the 'global' model that uses them as random effect.

#<---------------->#########
#CV~PCoA AXIS----------
load("pcoa_nonwoody_axis1.2")
head(pcoa_nonwoody_axis1.2)

#add axis to main table - averaged CV at mean_stability

#check spp names match
unique(row.names(pcoa_nonwoody_axis1.2)) #973 spp in the PCoA
unique(mean_stability$Group.1) #1494 single spp in the stability data
setdiff(unique(mean_stability$Group.1), unique(row.names(pcoa_nonwoody_axis1.2)) ) #778/1494 spp in veg database that are not in the PCoA axes
setdiff(unique(row.names(pcoa_nonwoody_axis1.2)), unique(mean_stability$Group.1)) #257/973 spp in PCoA that don't occur in veg database (after selection)

# mean_stability[, PCoA_axis1 := pcoa_nonwoody_axis1.2$Axis.1[match(mean_stability$Group.1, row.names(pcoa_nonwoody_axis1.2))]]
# mean_stability[, PCoA_axis2 := pcoa_nonwoody_axis1.2$Axis.2[match(mean_stability$Group.1, row.names(pcoa_nonwoody_axis1.2))]]

names(mean_stability)

mean_stability.1<-ungroup(mean_stability) %>% mutate(PCoA_axis1 = pcoa_nonwoody_axis1.2$Axis.1[match(mean_stability$Group.1, row.names(pcoa_nonwoody_axis1.2))])
mean_stability.1<-mean_stability.1 %>% mutate(PCoA_axis2 = pcoa_nonwoody_axis1.2$Axis.2[match(mean_stability$Group.1, row.names(pcoa_nonwoody_axis1.2))])
names(mean_stability.1)

mean_stability.1%>%pull(PCoA_axis1)%>%hist()
mean_stability.1%>%pull(PCoA_axis2)%>%hist()

# plot(sp_cv ~ PCoA_axis1 , stability_selected[woodiness_new=="non-woody",]) 
# plot(sp_cv ~ PCoA_axis2, stability_selected[woodiness_new=="non-woody",]) 

mean_stability_all<-ungroup(mean_stability_all) %>% mutate(PCoA_axis1 = pcoa_nonwoody_axis1.2$Axis.1[match(mean_stability_all$Group.1, row.names(pcoa_nonwoody_axis1.2))])
mean_stability_all<-mean_stability_all %>% mutate(PCoA_axis2 = pcoa_nonwoody_axis1.2$Axis.2[match(mean_stability_all$Group.1, row.names(pcoa_nonwoody_axis1.2))])




#Model--------
# Axis 1 - corresponds to foliar traits
# Axis 2 - corresponds to SM (and H)

mod_mean_z.cv.PCoA<- lmer(z.sp_cv ~ PCoA_axis1 +PCoA_axis2 +
                                (1|Group.1) + (1|Group.2), 
                              data= mean_stability.1%>%filter(woodyness_new=='non-woody'), REML = F)

summary(mod_mean_z.cv.PCoA)
coefplot(mod_mean_z.cv.PCoA)
r.squaredGLMM(mod_mean_z.cv.PCoA)
r.squaredGLMM(mod_mean_z.cv.PCoA)
# R2m       R2c
# 0.04603474 0.1655482

c.mean_z.cv.PCoA<-coefplot(mod_mean_z.cv.PCoA, intercept = F, title='', color="black",
                      pointSize = 5, lwdInner = 3, lwdOuter = 1.5, sort = 'magnitude')
c.mean_z.cv.PCoA<-c.mean_z.cv.PCoA+  labs(y="", x="Coefficent") + theme_classic() +
  scale_y_discrete(labels=c("PCoA_axis1"="1st Axis", 
                            "PCoA_axis2"="2nd Axis"))

#single PCoA axis against single traits

#PCoA1----
mod_mean_z.cv.PCoA1<- lmer(z.sp_cv ~ PCoA_axis1 +
                            (1|Group.1) + (1|Group.2), 
                          data= mean_stability.1%>%filter(woodyness_new=='non-woody'), REML = F)

summary(mod_mean_z.cv.PCoA1)
coefplot(mod_mean_z.cv.PCoA1)
r.squaredGLMM(mod_mean_z.cv.PCoA1)
# R2m       R2c
# 0.04353206 0.1641843


#PCoA2------

mod_mean_z.cv.PCoA2<- lmer(z.sp_cv ~ PCoA_axis2 +
                           (1|Group.1) + (1|Group.2), 
                           data= mean_stability.1%>%filter(woodyness_new=='non-woody'), REML = F)
#All together------
mod_mean_z.cv
tab.PCA<-export_summs(mod_mean_z.cv.PCoA, mod_mean_z.cv.PCoA1, mod_mean_z.cv.PCoA2,
                  mod_mean_cv.seedmass, mod_mean_cv.ldmc,
                  note=NULL ) # results = 'asis'

quick_xlsx(tab.PCA, file='tab.PCA.xlsx')
# tab$names[1]<-"GEN_MPD_Multi_trait"
number_format(tab.PCA)<-gsub("%.2f", "%.4f", number_format(tab))
# tab_final<-tab
# 
# tab_final<- add_rows(tab_final, tab)


#<---------------->#########
#ABUNDANCE & SD----------

names(mean_stability.1)
mod_mean_z.cv

#check abundances??
# mean_stability%>%filter(woodyness_new=='non-woody')%>%select(sp_mean_abu)
# lapply(mean_stability%>%filter(woodyness_new=='non-woody')%>%select(sp_mean_abu), hist)
# lapply(mean_stability%>%filter(woodyness_new=='non-woody', sp_mean_abu<100)%>%select(sp_mean_abu), hist)
# 
# hist(sp_mean_abu[sp_mean_abu<100])
# mean_stability[sp_mean_abu>200,]
# 
# mean_stability%>%filter(woodyness_new=='non-woody', sp_mean_abu>500)

# from tempo function
#x is a years x spp matrix
# # x_rel <- x/rowSums(x) single abundance/total abundance per year (across spp)
# sp_mean_abu <- apply(x, 2, mean) mean abundance across years
# sp_mean_abu_rel <- apply(x_rel, 2, mean) mean abundance across years using the relative abundances (considering the sum of abundances in each plot each year)
# sp_sd <- apply(x, 2, sd)
# sp_cv <- sp_sd/sp_mean_abu

# z.scores for abundance and sd-----
mean_stability.2<-mean_stability.1%>%group_by(Group.2)%>%mutate(z.sp_mean_abu=scale(sp_mean_abu),
                                                                z.sp_sd=scale(sp_sd))
pull(ungroup(mean_stability.2), z.sp_mean_abu)%>%hist()
pull(ungroup(mean_stability.2), z.sp_sd)%>%hist()
hist(scale(pull(ungroup(mean_stability.2), sp_sd)))

mean_stability_all<-mean_stability_all%>%group_by(Group.2)%>%mutate(z.sp_mean_abu=scale(sp_mean_abu),
                                                                z.sp_sd=scale(sp_sd))

#models------
mod_mean_z.sp_mean_abu<- lmer(z.sp_mean_abu ~ z.mean_LeafN + z.log.mean_SeedMass + 
                             z.log.mean_SLA +  z.mean_LDMC +
                             (1|Group.1) + (1|Group.2), 
                           data= mean_stability.2%>%filter(woodyness_new=='non-woody'), REML = F)

summary(mod_mean_z.sp_mean_abu)
r.squaredGLMM(mod_mean_z.sp_mean_abu)

mod_mean_sp_sd<- lmer(sp_sd ~ z.mean_LeafN + z.log.mean_SeedMass + 
                                z.log.mean_SLA +  z.mean_LDMC +
                                (1|Group.1) + (1|Group.2), 
                              data= mean_stability.2%>%filter(woodyness_new=='non-woody'), REML = F)
summary(mod_mean_sp_sd)
r.squaredGLMM(mod_mean_sp_sd)
coefplot(mod_mean_z.sp_sd)

#TAYLOR'S POWER LAW
attach(mean_stability.2)
plot(z.sp_mean_abu~z.sp_sd)
plot(z.sp_cv~log(z.sp_mean_abu))
plot(z.sp_cv~z.sp_mean_abu)
plot(log10(sp_mean_abu)~log10(sp_var))
abline(0,1)

cor.test(z.sp_cv, z.sp_mean_abu)

#plots-----------
c.mean_z.sp_mean_abu<-coefplot(mod_mean_z.sp_mean_abu, intercept = F, title='', color="black",
                      pointSize = 5, lwdInner = 3, lwdOuter = 1.5, sort = 'natural')
c.mean_z.sp_mean_abu<-c.mean_z.sp_mean_abu+  labs(y="", x="Coefficent") + theme_classic() +
  scale_y_discrete(labels=c("z.mean_LeafN"="Leaf N content", 
                            "z.log.mean_SeedMass"="Seed Mass",
                            "z.log.mean_SLA"="SLA",
                            "z.mean_LDMC"="LDMC"))
c.mean_z.sp_sd<-coefplot(mod_mean_z.sp_sd, intercept = F, title='', color="black",
                               pointSize = 5, lwdInner = 3, lwdOuter = 1.5, sort = 'natural')
c.mean_z.sp_sd<-c.mean_z.sp_sd+  labs(y="", x="Coefficent") + theme_classic() +
  scale_y_discrete(labels=c("z.mean_LeafN"="Leaf N content", 
                            "z.log.mean_SeedMass"="Seed Mass",
                            "z.log.mean_SLA"="SLA",
                            "z.mean_LDMC"="LDMC"))
# c_sm<-c_sm + xlim(-0.1,0.2)+ theme_classic() 

ggarrange(c.mean_z.sp_mean_abu, c.mean_z.sp_sd,labels = 'auto')

#<---------------->#########
#CATEGORICAL TRAITS----------
names(mean_stability.2)
mod_mean_z.cv

#woodiness, lifeform, lifespan, growthform

#woodiness-----
mod_wood<-lmer(z.sp_cv~woodyness_new  -1 +
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all, REML=F)
summary(mod_wood)

#add points with violin plot
est<-coefplot(mod_wood)$data
est$Coefficient<- as.factor(c("non-woody", "woody"))
m_wood<- ggplot(data = mean_stability_all[!is.na(mean_stability_all$woodyness_new)&!is.na(mean_stability_all$z.sp_cv),], aes(y=z.sp_cv, x=woodyness_new))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=3)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=1.5, inherit.aes = F, width = 0)+
  theme_classic(base_size=15) + xlab(label = "") +ylab(label = "species CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-0.1, 0.2), zoom.size = 1, zoom.data = z.sp_cv==0) #facet zoom on the right
m_wood


#life form-----
mod_lifeform<-lmer(z.sp_cv~lifeform -1 +  
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all, REML=F)
summary(mod_lifeform)

#add points with violin plot
est<-coefplot(mod_lifeform)$data
est$Coefficient<- c("Ch", "Cr", "H", "P", "T")
m_lifeform<- ggplot(data = mean_stability_all[!is.na(mean_stability_all$lifeform)&!is.na(mean_stability_all$z.sp_cv),], aes(y=z.sp_cv, x=lifeform))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=3)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=1.5, inherit.aes = F, width = 0)+
  theme_classic(base_size=15) + xlab(label = "") +ylab(label = "species CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-0.5, 1), zoom.size = 1, zoom.data = sp_cv==0) #facet zoom on the right



#life span-----
mod_lifespan<-lmer(z.sp_cv~lifespan_new -1 +  
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all, REML=F)
summary(mod_lifespan)

#add points with violin plot
est<-coefplot(mod_lifespan)$data
est$Coefficient<- c("annual", "not-annual")
m_lifespan<- ggplot(data = mean_stability_all[!is.na(mean_stability_all$lifespan_new)&!is.na(mean_stability_all$z.sp_cv),], aes(y=z.sp_cv, x=lifespan_new))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=3)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=1.5, inherit.aes = F, width = 0)+
  theme_classic(base_size=15) + xlab(label = "") +ylab(label = "species CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-0.25, 0.75), zoom.size = 1, zoom.data = sp_cv==0) #facet zoom on the right


# growth form-----
mod_growth<-lmer(z.sp_cv~growthform -1 +  
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all, REML=F)
summary(mod_growth)

#add points with violin plot
est<-coefplot(mod_growth)$data
est$Coefficient<-  c("fern", "graminoid", "herb", "herb/shrub", "shrub","shrub/tree", "tree")
m_growth<- ggplot(data = mean_stability_all[!is.na(mean_stability_all$growthform)&!is.na(mean_stability_all$z.sp_cv),], aes(y=z.sp_cv, x=growthform))+
  geom_jitter(colour="gray", size=0.5)+
  geom_violin(colour="black", position="identity", fill=NA ) +
  geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=3)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=1.5, inherit.aes = F, width = 0)+
  theme_classic(base_size=15) + xlab(label = "") +ylab(label = "species CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_zoom(ylim = c(-1, 1), zoom.size = 1, zoom.data = sp_cv==0) #facet zoom on the right


#All together-----
getwd()
tiff('categorical_traits.tiff', units="cm", width=16, height=10, res=300, compression = 'lzw')
ggarrange(m_wood, m_lifespan, m_lifeform, m_growth, labels = 'auto')
dev.off()

#write table with single categorical trait models
tab_cate<-export_summs(mod_lifespan, mod_lifeform, 
                  mod_wood, mod_growth,
                  note=NULL, model.names = c("Life span", "Life form", 
                                             "Woodyness", "Growth form") ) # results = 'asis'
#warning message In summ.merMod(model = new("lmerMod", resp = new("lmerResp", .xData = <environment>),  :
# Could not calculate r-squared. Try removing missing data
# before fitting the model.

r.squaredGLMM(mod_lifespan)
# R2m       R2c
# [1,] 0.04066684 0.2270277
r.squaredGLMM(mod_lifeform)
# R2m       R2c
# [1,] 0.06024877 0.1411521
r.squaredGLMM(mod_wood)
# R2m       R2c
# [1,] 7.035184e-07 0.228913
r.squaredGLMM(mod_growth)
# R2m       R2c
# [1,] 0.01753236 0.2243322


quick_xlsx(tab, file='tab.xlsx')
# tab$names[1]<-"GEN_MPD_Multi_trait"
number_format(tab)<-gsub("%.2f", "%.4f", number_format(tab))
# tab_final<-tab
# 
# tab_final<- add_rows(tab_final, tab)


#<---------------->#########
#WRITE TABELLONE--------
mean_stability_all
names(mean_stability_all)
mean_stability_all<- ungroup(mean_stability_all)%>%mutate(ID_dataset=ids_datasets$`Dataset ID`[match(mean_stability_all$Group.2, ids_datasets$Original_dataset)])
mean_stability_all[is.na(mean_stability_all$ID_dataset),]
mean_stability_all[is.na(mean_stability_all$ID_dataset),'ID_dataset']<- 'D52'
ids_datasets[ids_datasets$`Dataset ID`=="D62",]

write_csv(mean_stability_all, file = 'mean_stability_all.csv')

mean_stability_all[mean_stability_all$ID_dataset=='D52',] #22 entries

mean_stability_all[mean_stability_all$ID_dataset=='D60',] #D62
mean_stability_all[mean_stability_all$ID_dataset=='D61',] #D60
mean_stability_all[mean_stability_all$ID_dataset=='D62',] #D61


#trait correlations---------
names(mean_stability_all)
traits_temp<-mean_stability_all[,c(1,27:33)]
traits_temp<- distinct(traits_temp, Group.1, .keep_all = TRUE)
names(traits_temp)
names(traits_temp)[-1] <- c("Height", "Leaf N", "Leaf P", "Seed Mass", "SLA", "LDMC", "SSD")
# temp_traits<-Traits[!is.na(match(Traits$Species, unique(species))),c(2,5,8,11,14,17,20)]

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,  use = "na.or.complete")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}
pairs(traits_temp[,-1], lower.panel = panel.smooth, upper.panel=panel.cor)




