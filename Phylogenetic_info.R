#Adding Phylogenetic info into the analyses

#packages--------
getwd()
library(tidyverse)#managing data
library(ape) #plot.phylo
# install.packages("devtools") #not working
# library(devtools)
# devtools::install_github("GuangchuangYu/treeio")
library(phylobase) #for phylo tree manipulation
library(vegan) #for multivariate analyses
library(adephylo) #for phylogenetic distance
library(lme4)
library(MASS)
library(coefplot)
library(MuMIn)
library(lmerTest)
library(ggpubr)


#data from original analyses-----
# load("~/OneDrive - ČZU v Praze/STABILITY Project/Traits/workspace.RData")#nope?
mean_stability_all<-as_tibble(read.csv("mean_stability_all.csv"))


#load phylo data------
load("~/OneDrive - ČZU v Praze/STABILITY Project/Traits/Luisa_tree.RData")
summary(Luisa_tree)
head(Luisa_tree)
class(Luisa_tree) #phylo


#check tree
quartz()
plot.phylo(Luisa_tree, type='radial') #is it too big to render it!
length(Luisa_tree$tip.label)
mean_stability_all$Group.1


#convert tree format-----
Luisa_tree<-as(Luisa_tree, "phylo4") 
class(Luisa_tree) #phylo4
str(Luisa_tree)
Luisa_tree@label

#check differences between species names in main table and in tree
length(Luisa_tree@label)#4710 species in the tree
length(unique(mean_stability_all$Group.1)) #1799 species in the table
setdiff(mean_stability_all$Group.1, Luisa_tree@label) #49 spp not found in tree
setdiff(Luisa_tree@label,unique(mean_stability_all$Group.1)) #2960 spp not found in table

summary(mean_stability_all$Group.1=="mean_stability_all$Group.1")

mean_stability_all$Group.1[grep("Oxalis", mean_stability_all$Group.1)]

as.vector(Luisa_tree@label[grep("Lupinus",Luisa_tree@label)])

phylobase::subset(Luisa_tree, tips.include=as.vector(Luisa_tree@label[grep("Lupinus",Luisa_tree@label)]))
treePlot(phylobase::subset(Luisa_tree, tips.include=as.vector(Luisa_tree@label[grep("Lupinus",Luisa_tree@label)])),
         type="fan")

as.vector(Luisa_tree@label[grep("subsp",Luisa_tree@label)])
unique(mean_stability_all$Group.1)[!unique(mean_stability_all$Group.1)%in%sub_tree@label]
mean_stability_all[grep("Lycopodium", mean_stability_all$Group.1), "growthform"]
#Subspecies were turned to spp in table
#Ferns are missing from tree

#subset tree-----
length(Luisa_tree@label[Luisa_tree@label%in%unique(mean_stability_all$Group.1)])
sub_tree<-phylobase::subset(Luisa_tree, tips.include=Luisa_tree@label[Luisa_tree@label%in%unique(mean_stability_all$Group.1)])
treePlot(sub_tree, show.tip.label = F, type = "fan") #we need a better plot - circular and small label characters

save(sub_tree, file="subset_tree.RData")

#--------Using eigenvalues (PCA axes of phylogenetic distance)------------

load("NEW.dc.phylo.eig.r") 
head(dc.phylo.eig) #??  these are from the PCoA done on the residuals of the db RDA 
dim(dc.phylo.eig) #50 axes for 1202 observations
summary(is.na(dc.phylo.eig))
load("phylo.eig.r") # these are from the PCoA of the phylo
head(phylo.eig)
dim(phylo.eig) #15 axes for 1750 observations
names(phylo.eig)
row.names(phylo.eig)
summary(is.na(phylo.eig))

eigenval<-phylo.eig
eigenval$names<- row.names(eigenval)
summary(is.na(eigenval))
dim(eigenval)

#add all eigenvalues (15 axes selected from the 50 that explain 80% of variability in the phylogeny) to main table 

mean_stability_all_phylo<-left_join(mean_stability_all,eigenval, by=c('Group.1'='names'))
mean_stability_all_phylo$A50
mean_stability_all_phylo[c('Group.1', 'A1')]
names(mean_stability_all_phylo)
unique(mean_stability_all$Group.1) #1799 species

#Option 1 - Effects of traits beyond phylogeny - using residuals of a model of CV explained by eigenvalues-------
names(mean_stability_all_phylo)
mod_phylo<- lm(as.formula(paste('z.sp_cv', paste(names(mean_stability_all_phylo)[43:57], collapse=" + "), sep=" ~ ")), 
                 data=mean_stability_all_phylo) #the model doesn't include random effects
#careful with - 15 rows missing values in z.sp_cv, -71 rows missing values in eigenvalues 
#original dim 3884 observations, actual dim 3798
summary(mod_phylo)
dim(mod_phylo$model)  

#using spp as random effect
mod_phylo$terms
mod_phylo_random<- lmer(z.sp_cv ~ A1+ A4+A7+A10+A14+A19+A23+A24+A25+A31+A33+A36+A38+A43+A50+ 
                                (1|Group.1)+(1|Group.2), data=mean_stability_all_phylo) 
summary(mod_phylo_random) #Number of obs: 3798, groups:  Group.1, 1745; Group.2, 78

#make residual table with species names 
phylo_residuals<-as.tibble(resid(mod_phylo_random))
phylo_residuals$Group.1<-mean_stability_all_phylo$Group.1[!is.na(z.sp_cv)&!is.na(mean_stability_all_phylo['A1'])]
unique(phylo_residuals$Group.1) #1745 species
#add residuals to phylo table
mean_stability_all_phylo #A tibble: 3,884 × 57
unique(mean_stability_all_phylo$Group.1) #1799
phylo_residuals #A tibble: 3,798 × 2
unique(phylo_residuals$Group.1) #1745
# mean_stability_all_phylo$phylo_resid[!is.na(z.sp_cv)&!is.na(mean_stability_all_phylo['A1'])]<-phylo_residuals$value

#model using residuals i.e. portion of the CV variability left unexplained by the phylogeny
names(mean_stability_all_phylo)
mod_cv_phyloresid<- lmer(phylo_resid ~ z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                           z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD +
                           (1|Group.1) + (1|Group.2), data=mean_stability_all_phylo)
#boundary (singular) fit: see ?isSingular
summary(mod_cv_phyloresid) #Number of obs: 675, groups:  Group.1, 92; Group.2, 67
coefplot(mod_cv_phyloresid) #LDMC and Seed mass still negative, LDMC almost sign
r.squaredGLMM(mod_cv_phyloresid) #R2m 0.003629667 R2c 0.01499085

mod_cv_phyloresid_final<- lmer(phylo_resid ~  z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC +
                           (1|Group.1)+ (1|Group.2) , data=mean_stability_all_phylo)
# boundary (singular) fit: see ?isSingular
summary(mod_cv_phyloresid_final) #Number of obs: 1622, groups:  Group.1, 390; Group.2, 77
rePCA(mod_cv_phyloresid_final) #it is beccause of species random factor
coefplot(mod_cv_phyloresid_final)#LDMC and SLA significative 
r.squaredGLMM(mod_cv_phyloresid_final)#R2m 0.009478955 R2c 0.01978616

mod_cv_phyloresid_final.1<- lmer(phylo_resid ~  z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC +
                                   (1|Group.2) , data=mean_stability_all_phylo) #no singularity
summary(mod_cv_phyloresid_final.1) #Number of obs: 1622, groups:  Group.2, 77
coefplot(mod_cv_phyloresid_final.1)#LDMC and SLA significative
r.squaredGLMM(mod_cv_phyloresid_final.1) #R2m 0.009478955 R2c 0.01978616 - it stays the same!! 


#coefplot-------
p.mod_phyloresid.final<-coefplot(mod_cv_phyloresid_final, intercept = F, title='', color="black",
                                  pointSize = 3, lwdInner = 2, lwdOuter = 1, sort = 'natural')
temp.labels<-c( "z.mean_LeafN"  =   "Leaf N content",                       
                "z.log.mean_SeedMass" =  "Seed Mass",                    
                "z.log.mean_SLA" = "SLA" ,                           
                "z.mean_LDMC" = "LDMC" )

p.mod_phyloresid.final<-p.mod_phyloresid.final+  labs(y="", x="Coefficent") + theme_classic(base_size = 15) +
  scale_y_discrete(labels=temp.labels)




#Option 2 - Effects of traits AND phylogeny (excluding overlapping part)-------------
load("~/OneDrive - ČZU v Praze/STABILITY Project/Traits/workspace.RData")#to compare with original model
head(dc.phylo.eig)
decoup.eigenval<-dc.phylo.eig
decoup.eigenval$names<- row.names(dc.phylo.eig)
summary(is.na(decoup.eigenval))
dim(decoup.eigenval) #50 axes (instead of 10?)
match(decoup.eigenval$names, unique(mean_stability_all_phylo$Group.1))
names(decoup.eigenval)[1:50]<-paste(names(decoup.eigenval)[1:50], ".dc", sep = "")

names(mean_stability_all_phylo)
mean_stability_all_phylo[43:56] #selected axes from PCoA using phyl dist

mean_stability_all_phylo<-left_join(mean_stability_all_phylo,decoup.eigenval, by=c('Group.1'='names'))

# mod_phylo.dc <- lm(as.formula(paste('z.sp_cv', paste(names(mean_stability_all_phylo)[57:106], collapse=" + "), sep=" ~ ")), 
#                    data=mean_stability_all_phylo)
# final.mod_phylo.dc<-step(mod_phylo.dc)

#final model in the main text
mod_mean_z.cv
# R2m       R2c
# 0.06957713 0.1824155
#using same spp as in this model
unique(mod_mean_z.cv@frame$Group.1) #395 species

temp_data<- mean_stability_all_phylo%>%filter(Group.1%in%unique(mod_mean_z.cv@frame$Group.1))

#using just the first 10 axes
mod_phylo.dc_random<- lmer(z.sp_cv ~ A1.dc + A2.dc+  A3.dc+A4.dc+A5.dc+
                             A6.dc+A7.dc+A8.dc+A9.dc+ A10.dc+
                          (1|Group.1)+(1|Group.2), data=as.data.frame(temp_data) )
lmerTest::step(mod_phylo.dc_random, reduce.random=F) #only 1 axis maintained
mod_phylo.dc_random.final<- lmer(z.sp_cv ~ A6.dc + (1 | Group.1) + (1 | Group.2), data=as.data.frame(temp_data) )

summary(mod_phylo.dc_random.final)
coefplot(mod_phylo.dc_random.final)
r.squaredGLMM(mod_phylo.dc_random.final)
# R2m       R2c
# [1,] 0.007230586 0.1713662

mod_phylo.dc_traits<- lmer(z.sp_cv ~ z.mean_LeafN + z.log.mean_SeedMass  
                           +    z.log.mean_SLA    +      z.mean_LDMC + scale(A6.dc) + (1 | Group.1) + (1 | Group.2), data=as.data.frame(temp_data) )
summary(mod_phylo.dc_traits)
coefplot(mod_phylo.dc_traits)
r.squaredGLMM(mod_phylo.dc_traits)
# R2m       R2c
# [1,] 0.07437573 0.1886205

#coefplot-----
p.mod_phylo.dc_traits<-coefplot(mod_phylo.dc_traits, intercept = F, title='', color="black",
                                 pointSize = 3, lwdInner = 2, lwdOuter = 1, sort = 'natural')
temp.labels<-c( "z.mean_LeafN"  =   "Leaf N content",                       
                "z.log.mean_SeedMass" =  "Seed Mass",                    
                "z.log.mean_SLA" = "SLA" ,                           
                "z.mean_LDMC" = "LDMC",
                "scale(A6.dc)" = "Decoupled Phylo" )

p.mod_phylo.dc_traits<-p.mod_phylo.dc_traits+  labs(y="", x="Coefficent") + theme_classic(base_size = 15) +
  scale_y_discrete(labels=temp.labels)


#plot together------
# getwd()
# tiff('interaction.life_span.tiff', units="cm", width=33, height=11, res=300, compression = 'lzw')
ggarrange(p.mod_phyloresid.final, p.mod_phylo.dc_traits, labels = 'auto')
# dev.off()


