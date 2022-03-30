# following are analyses on samples for which either IgA-seq was performed and both IgA-bound and IgA-unbound bacterial fractions were profiled, or for samples in which autologous paired serum and stool were used and IgG-bound and IgG-unbound bacterial fractions were profiled.

library(labdsv);library(ggplot2);library(lme4);library(lmerTest);library(plyr);library(stringr);library(RColorBrewer);library(Hmisc);library(viridis);library(vegan);library(phyloseq);library(superheat)
library(ggbeeswarm);library(ape);library(reshape2);library(splitstackshape)

# read in metadata mapping file, with rows as sample IDs:
map<-read.table(file="~/Dropbox/post-doc/HIV/bacFACS/170124 seq of BM6-9_Trinch/170216_HIVBM1-9_map_forR.txt",sep="\t",header=T,row.names=1)

# read in dada2-processed, rarefied ASV table, with columns as sample IDs:
tab1<-read.table("~/Dropbox/post-doc/HIV/bacFACS/validation_experiments/180712_MBM1-4_16s_data/seqtab_mc005p_rar23k.txt",sep="\t",header=T,row.names=1)

# filter out ASVs with read abundance <0.01% of total:
tab2<-tab1[rowSums(tab1)>sum(tab1)*0.0001,]

# filter out samples in mapping file not present in ASV table:
map2<-map[rownames(map)%in%colnames(tab2),]

#filter out samples missing either a negative fraction (Ig-unbound) or positive (Ig-bound) fraction; 'expmouse' refers column in metadata that defines ID of the experimental subject, 'negpos' refers to whether sample is negative or positive fraction:
sampsw2sum1<-table(map2$negpos,map2$expmouse)
sampsw2sum<-colSums(sampsw2sum1[c("neg","pos"),])
sampsw2<-names(sampsw2sum [sampsw2sum>1])
map3<-map2[map2$expmouse%in%sampsw2,]



# add pseudocount of 1:
tab<-(tab2+1)

# define which immunoglobulin class to query, if both IgA-seq and IgG-seq was performed, with column named 'Ig' in this case presenting this information and in this case selecting samples for IgG: 
map4igg<-subset(map3, Ig=="IgG")

# create Ig score matrix:
ICmat<-matrix(nrow=nrow(tab2),ncol=unique(map4igg$expmouse))
colnames(ICmat)<-unique(map4igg$expmouse)
rownames(ICmat)<-rownames(tab2)
for(i in 1:length(unique(map4igg$expmouse)))
{
  subjID<-unique(map4igg$expmouse)[i]
  mapsub<-subset(map4igg, expmouse ==subjID)
  possamp<-rownames(subset(mapsub, negpospre =="pos"))
  negsamp<-rownames(subset(mapsub, negpospre =="neg"))
  
#remove ASVs that were not detected in either pos or neg fractions by forcing "NaN" at log transformation calculation:
  zeroesinpos<-(tab[,possamp]==1)
  zeroesinneg<-(tab[,negsamp]==1)
  bothzeroes<-zeroesinpos+ zeroesinneg
  bothzeropositions<-(bothzeroes==2)
  tab[bothzeropositions,possamp]<-0
  tab[bothzeropositions,negsamp]<-0

# calculate log ratios for all taxa for the given pair of pos and neg samples:
  ICmat[,i]<-log(tab[,possamp],10)-log(tab[,negsamp],10)
}
