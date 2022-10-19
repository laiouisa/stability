#Adding Phylogenetic info into the analyses

getwd()
library(tidyverse)#managing data
library(ape) #plot.phylo
# install.packages("devtools") #not working
# library(devtools)
# devtools::install_github("GuangchuangYu/treeio")
library(phylobase)


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

load("dc.phylo.eig.r") 
head(dc.phylo.eig) #?? are these from the db RDA??
dim(dc.phylo.eig) #50 axes for #461 observations??
summary(is.na(dc.phylo.eig))
load("phylo.eig.r")
head(phylo.eig)
dim(phylo.eig) #15 axes for 1750 observations
names(phylo.eig)
row.names(phylo.eig)
summary(is.na(phylo.eig))

eigenval<-phylo.eig
eigenval$names<- row.names(eigenval)
summary(is.na(eigenval))

#add all eigenvalues (15 axes selected from the 50 that explain 80% of variability in the phylogeny) to main table 

mean_stability_all_phylo<-left_join(mean_stability_all,eigenval, by=c('Group.1'='names'))
mean_stability_all_phylo$A50
mean_stability_all_phylo[c('Group.1', 'A1')]
names(mean_stability_all_phylo)

#Option 1 - Effects of traits beyond phylogeny - using residuals of a model of CV explained by eigenvalues
names(mean_stability_all_phylo)
mod_phylo<- lm(as.formula(paste('z.sp_cv', paste(names(mean_stability_all_phylo)[43:57], collapse=" + "), sep=" ~ ")), 
                 data=mean_stability_all_phylo) #the model doesn't include random effects
#careful with - 15 rows missing values in z.sp_cv, -71 rows missing values in eigenvalues 
#original dim 3884 observations, actual dim 3798
summary(mod_phylo)
dim(mod_phylo$model)  

#IT ALL WORKS BETTER IF I USE LESS AXES (StepAIC on the 15 axes..)

#using spp as random effect
mod_phylo$terms
mod_phylo_random<- lmer(z.sp_cv ~ A1+ A4+A7+A10+A14+A19+A23+A24+A25+A31+A33+A36+A38+A43+A50+ 
                                (1|Group.1)+(1|Group.2), data=mean_stability_all_phylo) 
summary(mod_phylo_random) #Number of obs: 3798, groups:  Group.1, 1745; Group.2, 78

#make residual table with species names 
phylo_residuals<-as.tibble(resid(mod_phylo_random))
phylo_residuals$species<-mean_stability_all_phylo$Group.1[!is.na(z.sp_cv)&!is.na(mean_stability_all_phylo['A1'])]
unique(phylo_residuals$species) #1745 species
#add residuals to phylo table
mean_stability_all_phylo<-left_join(mean_stability_all_phylo, phylo_residuals, by=c('Group.1'='species'))
colnames(mean_stability_all_phylo)[58]<-'phylo_residuals'

#model using residuals i.e. portion of the CV variability left unexplained by the phylogeny
names(mean_stability_all_phylo)
mod_cv_phyloresid<- lmer(phylo_residuals ~ z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
                           z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD +
                           (1|Group.1) + (1|Group.2), data=mean_stability_all_phylo)
#boundary (singular) fit: see ?isSingular
summary(mod_cv_phyloresid) #Number of obs: 8595, groups:  Group.1, 92; Group.2, 67
coefplot(mod_cv_phyloresid)
r.squaredGLMM(mod_cv_phyloresid) #R2m 0.001610543 R2c 0.007688418

mod_cv_phyloresid_final<- lmer(phylo_residuals ~  z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC +
                           (1|Group.1)+ (1|Group.2) , data=mean_stability_all_phylo)
# boundary (singular) fit: see ?isSingular
summary(mod_cv_phyloresid_final) #Number of obs: 14124, groups:  Group.1, 390; Group.2, 77
rePCA(mod_cv_phyloresid_final) #it is beccause of datasets 
coefplot(mod_cv_phyloresid_final)
r.squaredGLMM(mod_cv_phyloresid_final)#R2m 0.00423158 R2c 0.01560563

mod_cv_phyloresid_final.1<- lmer(phylo_residuals ~  z.mean_LeafN + z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC +
                                 (1|Group.1) , data=mean_stability_all_phylo) #no singularity
summary(mod_cv_phyloresid_final.1) #Number of obs: 14124, groups:  Group.1, 390
coefplot(mod_cv_phyloresid_final.1)
r.squaredGLMM(mod_cv_phyloresid_final.1) #R2m 0.00423158 R2c 0.01560564 - it stays the same!! 



