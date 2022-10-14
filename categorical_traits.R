#Models using separating annuals and perennials
#models using quantitative traits in interaction with lifespan 
#models using quantitative traits in addition to lifespan and growthform

#packages------
library(tidyverse)#managing data
library(lme4)#Mixed Effects models
library(MuMIn) #r.squared + confidence intervals with visreg (replaces predict.merMod?)
library(car)
library(coefplot)

#data-------
#data from original analyses
# load("~/OneDrive - ČZU v Praze/STABILITY Project/Traits/workspace.RData")#nope?
mean_stability_all<-as_tibble(read.csv("mean_stability_all.csv"))
names(mean_stability_all)

#Categorical traits names
# "growthform"  "woodyness" "woodyness_new" "lifeform" "lifespan" "lifespan_new"    

#check categorical trait models


#WOODYNESS------------
summary(as.factor(mean_stability_all$woodyness))
# non-woody non-woody/woody           woody 
# 3271             126             487 
summary(is.na(mean_stability_all$woodyness)) #no NAs
summary(as.factor(mean_stability_all$woodyness_new))


mod_wood<-lmer(z.sp_cv~woodyness_new  -1 +
                 (1|Group.1) + (1|Group.2),
               mean_stability_all, REML=F)
summary(mod_wood)
r.squaredGLMM(mod_wood)
coefplot(mod_wood)


mod_wood1<-lmer(z.sp_cv~woodyness  -1 +
                  (1|Group.1) + (1|Group.2),
                mean_stability_all, REML=F)
summary(mod_wood1)
r.squaredGLMM(mod_wood1)
coefplot(mod_wood1)

#non-woody/woody species have a highly variable CV!


#LIFE FORM------------
summary(as.factor(mean_stability_all$lifeform)) #1381 NAs
mod_lifeform<-lmer(z.sp_cv~lifeform -1 +  
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all, REML=F)
summary(mod_lifeform)
coefplot(mod_lifeform)
r.squaredGLMM(mod_lifeform)


#LIFE SPAN----------
summary(as.factor(mean_stability_all$lifespan)) #1143 NAs (715 spp that were defined as not-annual)
summary(as.factor(mean_stability_all$lifespan_new)) #no Nas?
# annual             not-annual 
# 606  (283 spp)     3278 (1516 spp)

unique(mean_stability_all[is.na(mean_stability_all$lifespan),"Group.1"])
unique(mean_stability_all[mean_stability_all$lifespan_new=="annual","Group.1"]) #283 spp
unique(mean_stability_all[mean_stability_all$lifespan_new=="not-annual","Group.1"]) #1516 spp

mod_lifespan<-lmer(z.sp_cv~lifespan_new -1 +  
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all, REML=F)
summary(mod_lifespan)
coefplot(mod_lifespan)
r.squaredGLMM(mod_lifespan)

#GROWTH FORM------
summary(as.factor(mean_stability_all$growthform)) 
mod_growth<-lmer(z.sp_cv~growthform -1 +  
                   (1|Group.1) + (1|Group.2),
                 mean_stability_all, REML=F)
summary(mod_growth)
coefplot(mod_growth)
r.squaredGLMM(mod_growth)


#////MODELS-------

#Separating by lifespan_new (annuals vs not-annuals)-----
names(mean_stability_all)
#continous traits in the full model z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
# z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD  
#continous traits in the final model? 
mod_annuals<- lmer(z.sp_cv~ z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                     z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD +
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all[mean_stability_all$lifespan_new=='annual',], REML=F)
# boundary (singular) fit: see ?isSingular
summary(mod_annuals)
coefplot(mod_annuals.1)
r.squaredGLMM(mod_annuals) #not reliable
# R2m       R2c
# [1,] 0.3084379 0.3084379
mod_annuals.1<- lmer(z.sp_cv~ z.mean_LeafN  +  z.mean_LDMC + 
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all[mean_stability_all$lifespan_new=='annual'&mean_stability_all$woodyness_new=='non-woody',], REML=F)
r.squaredGLMM(mod_annuals.1)
# R2m        R2c
# [1,] 0.01009974 0.08823534

mod_notannuals<- lmer(z.sp_cv~ z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                     z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD +
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all[mean_stability_all$lifespan_new=='not-annual'&mean_stability_all$woodyness_new=='non-woody',], REML=F)
summary(mod_notannuals)
coefplot(mod_notannuals)
r.squaredGLMM(mod_notannuals)
# R2m        R2c
# [1,] 0.03330625 0.08533796

mod_notannuals.1<- lmer(z.sp_cv~ z.log.mean_SeedMass + z.mean_LDMC + z.log.mean_height +
                        (1|Group.1) + (1|Group.2),
                      mean_stability_all[mean_stability_all$lifespan_new=='not-annual'&mean_stability_all$woodyness_new=='non-woody',], REML=F)
summary(mod_notannuals.1)
coefplot(mod_notannuals.1)
r.squaredGLMM(mod_notannuals.1)
# R2m       R2c
# [1,] 0.03811708 0.1281322


#Quantitative traits in interaction with life span----

mod_interaction<- lmer(z.sp_cv~ (z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                        z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD)*lifespan_new +
                        (1|Group.1) + (1|Group.2),
                      mean_stability_all, REML=F)
summary(mod_interaction)
coefplot(mod_interaction)
r.squaredGLMM(mod_interaction)
# R2m       R2c
# [1,] 0.08069282 0.1188175


mod_interaction.1<- lmer(z.sp_cv~ ( z.mean_LeafN )*lifespan_new +  z.mean_LDMC -1 +
                         (1|Group.1) + (1|Group.2), mean_stability_all, REML=F)
summary(mod_interaction.1)
coefplot(mod_interaction.1)
r.squaredGLMM(mod_interaction.1)
# R2m      R2c
# [1,] 0.09444683 0.176568
mod_interaction.2<- lmer(z.sp_cv~ lifespan_new + z.mean_LeafN  +  z.log.mean_SeedMass +  z.mean_LDMC +z.log.mean_SLA-1 +
                           (1|Group.1) + (1|Group.2), mean_stability_all, REML=F)
summary(mod_interaction.2)
coefplot(mod_interaction.2)
r.squaredGLMM(mod_interaction.2)
# R2m       R2c
# [1,] 0.08979585 0.1685861

#Quant traits in addition to lifespan and growth form--------
mod_addition<- lmer(z.sp_cv~ z.log.mean_height + z.mean_LeafN  +    
                                   z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + 
                                  lifespan_new + graminoid  +
                                   (1|Group.1) + (1|Group.2),
                       mean_stability_all, REML=F)
summary(mod_addition) 
coefplot(mod_addition)
r.squaredGLMM(mod_addition)
# R2m       R2c
# [1,] 0.06947318 0.1157023
levels(as.factor(mod_addition@frame$growthform)) #were are the other growthform categories???
# [1] "graminoid"  "herb"       "herb/shrub"
levels(as.factor(mean_stability_all$growthform)) 
#categories dropped because of missing trait values?

mean_stability_all$graminoid<-'not-graminoid'
mean_stability_all$graminoid[mean_stability_all$growthform=='graminoid']<-'graminoid'


load("dc.phylo.eig.r") 

