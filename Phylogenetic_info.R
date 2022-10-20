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
head(dc.phylo.eig) #?? are these from the db RDA?? detrended correspondence analyses?
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
mean_stability_all_phylo$phylo_resid[!is.na(z.sp_cv)&!is.na(mean_stability_all_phylo['A1'])]<-phylo_residuals$value

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



#Option 2 - Effects of traits AND phylogeny (excluding overlapping part)-------------

mean_stability_all_phylo[43:57] #selected axes from PCoA using phyl dist
sub_tree #tree - from which I can calculate distances - which distance method?
dist_phylo<-distTips(sub_tree)
head(dist_phylo) 
length(dist_phylo)
labels(dist_phylo) #1750 species

#dbRDA Pcoa phylo~traits
rda(dist_phylo~z.log.mean_height + z.mean_LeafN + z.mean_LeafP +    
  z.log.mean_SeedMass + z.log.mean_SLA + z.mean_LDMC + z.mean_SSD, data= mean_stability_all_phylo, na.action=na.omit) 
# Error in qr.fitted(Q, Y) : NA/NaN/Inf in foreign function call (arg 5) because of duplicates??
# Error in X[nas, , drop = FALSE] : subscript out of bounds

names(mean_stability_all_phylo)
mean_stability_all_phylo %>% drop_na(names(mean_stability_all_phylo)[27:33]) #NAs dropped directly on the function L148


#usare file dc.phylo.eig - selezionare primi 10 assi lme(cv~traits+10 assi, random)  da confrontare con lme(cv ~traits) [sulla stessa quantita di species]

#aggiungere 10 assi su mean_stability_all_phylo
dc.phylo.eig$Group.1<-row.names(dc.phylo.eig) #461 spp -- too little 
mean_stability_all_phylo%>%left_join(dc.phylo.eig, by='Group.1')




