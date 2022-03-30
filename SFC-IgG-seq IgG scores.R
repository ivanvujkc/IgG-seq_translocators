library(ggplot2);library(lme4);library(lmerTest);library(plyr);library(stringr);library(RColorBrewer);library(Hmisc);library(viridis);library(vegan);library(phyloseq)
library(ggbeeswarm);library(ape)

# read in dada2-generated ASV table, ASV's as rows (with md5 hash of the fasta sequence as ASV identifiers) and samples as columns:
tab1a<-read.table(file="~/home/seqtab_md5.txt",header=T,row.names=1,sep="\t")

# read in metadata, sample names as row names:
metadat<-read.csv("~/home/metadata.csv",row.names=1)

# remove taxa less than 0.001% total read abundance:
tab1b<-tab1a[rowSums(tab1a)>sum(tab1a)*0.00001,]
# get rid of taxa in fewer than 30% of samples:
proppos<-as.data.frame(cbind(proppositive=apply(tab1b, 1,function(x) sum(x>0)/length(x)),blank=NA))
belowthresh<-rownames(subset(proppos, proppositive<0.3))
tab2<-tab1b[!rownames(tab1b)%in%belowthresh,]

# rarefy samples:
ps <- phyloseq(otu_table(t(tab2), taxa_are_rows = FALSE))
tab2rar<-t(rarefy_even_depth(ps,sample.size=400000,replace=F))

# perform 'pseudo-absolute' normalization by calculating ratio to known spike-in control, Staphylococcus hominis, which is absent from SFC and an equivalent quantity was added to each sample. S.hominis md5 hash identifier: f46a7ca244afef522b22a11bd33d27b1
tab3shom<-apply(tab2rar, 2, function(x) {x/x[rownames(tab3)=="f46a7ca244afef522b22a11bd33d27b1"]})

# add pseudocount of minimal non-zero value, log transform:
tab3shom2<-tab3shom+min(tab3shom[tab3shom!=0])
logtab<-log(tab3shom2,base=2)

# perform IgG score calculations, taking ratios with ASV abundance in the pre-enrichment fraction ('no serum' control) as the numerator, and ASV abundance in the negative fractions (IgG-unbound fractions) as the denominator. define the string that is present in the 'no serum' control sample names within the 'grep' argument, in the following case it is "2noAb":
noabscores1<-logtab[,grep("2noAb",colnames(logtab))]
noabscores2<-apply(noabscores1, 1, function(x) mean(as.numeric(x)[is.finite(as.numeric(x))],na.rm=T))
logtab3b<-apply(logtab, 2,function(x) (noabscores2-x))
