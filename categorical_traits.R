#Models using separating annuals and perennials
#models using quantitative traits in interaction with lifespan 
#models using quantitative traits in addition to lifespan and growthform

#packages------
library(tidyverse)#managing data
library(lme4)#Mixed Effects models
library(MuMIn) #r.squared + confidence intervals with visreg (replaces predict.merMod?)
library(car)
library(coefplot)
library(jtools)
library(huxtable)

#data-------
#data from original analyses
# load("~/OneDrive - ČZU v Praze/STABILITY Project/Traits/workspace.RData")#nope?
mean_stability_all<-as_tibble(read.csv("mean_stability_all.csv"))
names(mean_stability_all)
dim(mean_stability_all)
# observations: 3884 (species x plot)
length(unique(mean_stability_all$Group.1)) #1799 species

summary(is.na(mean_stability_all$z.sp_cv)) #15 entries with no standardized CV
mean_stability_all[is.na(mean_stability_all$z.sp_cv),c("Group.1", "Group.2","sp_cv" , "residency_prop", "z.sp_cv")]

#Categorical traits names
# "growthform"  "woodyness" "woodyness_new" "lifeform" "lifespan" "lifespan_new"    

#check categorical trait models

attach(mean_stability_all)

#WOODYNESS------------
summary(as.factor(mean_stability_all$woodyness))
# non-woody non-woody/woody           woody 
# 3271             126             487 
summary(is.na(mean_stability_all$woodyness)) #no NAs
summary(as.factor(mean_stability_all$woodyness_new))
length(unique(mean_stability_all$Group.1[mean_stability_all$woodyness_new=="woody"]))
# non-woody                woody 
# 3397  (1494 species)     487 (305 species)
summary(is.na(mean_stability_all$woodyness_new)) #no NAs

mod_wood<-lmer(z.sp_cv~woodyness_new  -1 +
                 (1|Group.1) + (1|Group.2),
               mean_stability_all, REML=F)
summary(mod_wood)
r.squaredGLMM(mod_wood)
#               R2m      R2c
# [1,] 7.035184e-07 0.228913
coefplot(mod_wood)

# source("/Users/laiouisa/OneDrive - ČZU v Praze/R/R scripts/facetzoom2.R") #for facet zoom on the right 

# est<-coefplot(mod_wood)$data
# est$Coefficient<- as.factor(c("non-woody", "woody"))
m_wood<- ggplot(data = mean_stability_all[!is.na(mean_stability_all$woodyness_new)&!is.na(mean_stability_all$z.sp_cv),], 
                aes(y=z.sp_cv, x=woodyness_new, color=lifespan_new))+
  geom_jitter(size=0.5)+ #colour="gray", 
  geom_violin(colour="black", position="identity", fill=NA ) +
  # geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=3)+
  # geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  # geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=1.5, inherit.aes = F, width = 0)+
  # facet_zoom(ylim = c(-0.1, 0.2), zoom.size = 1, zoom.data = z.sp_cv==0) #facet zoom on the right
  theme_classic(base_size=15) + xlab(label = "") +ylab(label = "species CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
m_wood


mod_wood1<-lmer(z.sp_cv~woodyness  -1 +
                  (1|Group.1) + (1|Group.2),
                mean_stability_all, REML=F)
summary(mod_wood1)
r.squaredGLMM(mod_wood1)
coefplot(mod_wood1)

#non-woody/woody species have a highly variable CV! - focusing on non-woody species


#LIFE FORM------------
summary(as.factor(mean_stability_all$lifeform)) #1381/3884 NAs - 806/1799 species
summary(as.factor(woodyness_new[is.na(lifeform)]))
length(unique(Group.1[is.na(lifeform)]))
mod_lifeform<-lmer(z.sp_cv~lifeform -1 +  
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all, REML=F)
summary(mod_lifeform)
coefplot(mod_lifeform)
r.squaredGLMM(mod_lifeform)


#LIFE SPAN----------
summary(as.factor(mean_stability_all$lifespan)) #1143 NAs (715 spp that were defined as not-annual)
summary(as.factor(mean_stability_all$lifespan_new)) #no Nas?
length(unique(Group.1[lifespan_new=='annual']))
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

#lifespan and other traits
mean_stability_all[mean_stability_all$lifespan_new=='annual',] #606 observations 
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.sp_cv"] ) #No NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.log.mean_height"] ) #9 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.mean_LeafN"] ) #343 NAs !!
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.mean_LeafP"] ) #372 NAs !!
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.log.mean_SeedMass"] ) #23 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.log.mean_SLA"] ) #115 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.mean_LDMC"] ) #199 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='annual',"z.mean_SSD"] ) #562 NAs !!

#lifespan and other traits
mean_stability_all[mean_stability_all$lifespan_new=='not-annual',] #3278 observations 
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.sp_cv"] ) #15 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.log.mean_height"] ) #1050 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.mean_LeafN"] ) #1602 NAs !!
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.mean_LeafP"] ) #1862 NAs !!
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.log.mean_SeedMass"] ) #873 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.log.mean_SLA"] ) #1295 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.mean_LDMC"] ) #1595 NAs
summary(mean_stability_all[mean_stability_all$lifespan_new=='not-annual',"z.mean_SSD"] ) #2618 NAs !!


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
coefplot(mod_annuals)
r.squaredGLMM(mod_annuals) #not reliable
# R2m       R2c
# [1,] 0.3084379 0.3084379
                           

# mod_annuals.1<- lmer(z.sp_cv~ z.mean_LeafN  +  z.mean_LDMC + 
#                      (1|Group.1) + (1|Group.2),
#                    mean_stability_all[mean_stability_all$lifespan_new=='annual'&mean_stability_all$woodyness_new=='non-woody',], REML=F)
# r.squaredGLMM(mod_annuals.1)
# R2m        R2c
# [1,] 0.01009974 0.08823534


mod_annuals_final<- lmer(z.sp_cv~   z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + 
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all[mean_stability_all$lifespan_new=='annual',], REML=F)
summary(mod_annuals_final)
coefplot(mod_annuals_final)
r.squaredGLMM(mod_annuals_final) 

# ---


mod_notannuals<- lmer(z.sp_cv~ z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                     z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD +
                     (1|Group.1) + (1|Group.2),
                   mean_stability_all[mean_stability_all$lifespan_new=='not-annual',], REML=F)
summary(mod_notannuals)
coefplot(mod_notannuals)
r.squaredGLMM(mod_notannuals)
# R2m        R2c
# [1,] 0.03330625 0.08533796

# mod_notannuals.1<- lmer(z.sp_cv~ z.log.mean_SeedMass + z.mean_LDMC + z.log.mean_height +
#                         (1|Group.1) + (1|Group.2),
#                       mean_stability_all[mean_stability_all$lifespan_new=='not-annual'&mean_stability_all$woodyness_new=='non-woody',], REML=F)
# summary(mod_notannuals.1)
# coefplot(mod_notannuals.1)
# r.squaredGLMM(mod_notannuals.1)
# R2m       R2c
# [1,] 0.03811708 0.1281322

mod_notannuals_final<- lmer(z.sp_cv~  z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + 
                        (1|Group.1) + (1|Group.2),
                      mean_stability_all[mean_stability_all$lifespan_new=='not-annual',], REML=F)
summary(mod_notannuals_final)
coefplot(mod_notannuals_final)
r.squaredGLMM(mod_notannuals_final)

#Quantitative traits in interaction with life span----

mod_interaction<- lmer(z.sp_cv~ (z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                        z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD)*lifespan_new -1 +
                        (1|Group.1) + (1|Group.2),
                      mean_stability_all, REML=F)
summary(mod_interaction) #Number of obs: 676, groups:  Group.1, 93; Group.2, 67
coefplot(mod_interaction)
r.squaredGLMM(mod_interaction)
# R2m       R2c
# [1,] 0.08069282 0.1188175

mod_interaction@frame$Group.1 #only 93 levels - spp for which we have all the traits!

mod_interaction.final<-  lmer(z.sp_cv~ ( z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC)*lifespan_new -1 +
                                (1|Group.1) + (1|Group.2),
                              mean_stability_all, REML=F)
summary(mod_interaction.final) #Number of obs: 1630, groups:  Group.1, 395; Group.2, 77
coefplot(mod_interaction.final) #interaction between not-annuals and leafN keeps comming out 
r.squaredGLMM(mod_interaction.final) #R2m 0.1041634 R2c 0.1862911

# annual             not-annual 
# 606  (283 spp)     3278 (1516 spp)
summary(z.mean_LeafN[lifespan_new=='annual']) #343 NAs
summary(z.mean_LeafN[lifespan_new=='not-annual']) #1602 NAs 


p<- ggplot(data = mean_stability_all[!is.na(mean_stability_all$z.sp_cv),], 
                aes(y=z.mean_LeafN, x=lifespan_new))+
  geom_jitter(size=0.5)+ #colour="gray", 
  geom_violin(colour="black", position="identity", fill=NA ) +
  # geom_point(data=est, mapping = aes(y=Value, x=Coefficient), color="red", size=3)+
  # geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowOuter, ymax = HighOuter) , color="red", size=0.5, inherit.aes = F, width = 0)+
  # geom_errorbar(data=est, mapping = aes( x=Coefficient, ymin = LowInner, ymax = HighInner), color="red", size=1.5, inherit.aes = F, width = 0)+
  # facet_zoom(ylim = c(-0.1, 0.2), zoom.size = 1, zoom.data = z.sp_cv==0) #facet zoom on the right
  theme_classic(base_size=15) + xlab(label = "") +ylab(label = "species CV") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p

# mod_interaction.1<- lmer(z.sp_cv~ ( z.mean_LeafN )*lifespan_new +  z.mean_LDMC -1 +
#                          (1|Group.1) + (1|Group.2), mean_stability_all, REML=F)
# summary(mod_interaction.1)
# coefplot(mod_interaction.1)
# r.squaredGLMM(mod_interaction.1)
# # R2m      R2c
# # [1,] 0.09444683 0.176568



mod_additive<- lmer(z.sp_cv~ lifespan_new + z.mean_LeafN  +  z.log.mean_SeedMass +  z.mean_LDMC +z.log.mean_SLA-1 +
                           (1|Group.1) + (1|Group.2), mean_stability_all, REML=F)
summary(mod_additive)
coefplot(mod_additive)
r.squaredGLMM(mod_additive)
# R2m       R2c
# [1,] 0.08979585 0.1685861

#Quant traits in addition to lifespan and graminoid--------

mean_stability_all$graminoid<-'not-graminoid'
mean_stability_all$graminoid[mean_stability_all$growthform=='graminoid']<-'graminoid'

mod_addition<- lmer(z.sp_cv~ z.log.mean_height + z.mean_LeafN  +    
                                   z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + 
                                  lifespan_new + graminoid  -1 +
                                   (1|Group.1) + (1|Group.2),
                       mean_stability_all, REML=F)
summary(mod_addition) 
coefplot(mod_addition) #graminoid category is dropped
r.squaredGLMM(mod_addition)
# R2m       R2c
# [1,] 0.06947318 0.1157023




