#Adding Phylogenetic info into the analyses

getwd()
library(tidyverse)#managing data
library(ape) #plot.phylo
# install.packages("devtools") #not working
# library(devtools)
# devtools::install_github("GuangchuangYu/treeio")
library(phylobase)


#data from original analyses
# load("~/OneDrive - ČZU v Praze/STABILITY Project/Traits/workspace.RData")#nope?
mean_stability_all<-as_tibble(read.csv("mean_stability_all.csv"))


#load phylo data
load("~/OneDrive - ČZU v Praze/STABILITY Project/Traits/Luisa_tree.RData")
summary(Luisa_tree)
head(Luisa_tree)
class(Luisa_tree) #phylo


#check tree
quartz()
plot.phylo(Luisa_tree, type='radial') #is it too big to render it!
length(Luisa_tree$tip.label)
mean_stability_all$Group.1


#convert tree format
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

#subset tree
length(Luisa_tree@label[Luisa_tree@label%in%unique(mean_stability_all$Group.1)])
sub_tree<-phylobase::subset(Luisa_tree, tips.include=Luisa_tree@label[Luisa_tree@label%in%unique(mean_stability_all$Group.1)])
treePlot(sub_tree, show.tip.label = F, type = "fan") #we need a better plot - circular and small label characters

save(sub_tree, file="subset_tree.RData")



