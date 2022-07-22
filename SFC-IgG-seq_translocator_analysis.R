
library(ggplot2);library(lme4);library(lmerTest);library(plyr);library(stringr);library(RColorBrewer);library(Hmisc);library(viridis);library(vegan)
library(ggbeeswarm);library(ape)

# translocator culture data:
library(reshape2)
# read in NCBI 16S BLAST results for SFC ASVs:
blastax<-read.table(file="~/Downloads/fig3_SFC_blastax.txt",sep="\t")
blastax[,4]<-gsub("[[]", "", blastax[,4])
blastax[,4]<-gsub("[]]", "", blastax[,4])
blastax2<-colsplit(blastax[,4],pattern=" ",1:10)
blastax3<-cbind(blastax,genspec=paste(blastax2[,1],blastax2[,2], sep=" "))
blastax4<-blastax3[order(blastax3[,3],decreasing=T),]


# read in translocator data for IBD cohort:
xlocs<-read.csv(file="~/Downloads/fig3_Devkota IgG-seq- cultivation and LBP data 1.25.21.csv",header=T)
xlocs2<-matrix(nrow=nrow(xlocs),ncol=ncol(xlocs))
colnames(xlocs2)<-colnames(xlocs)
# rename e. coli to e. marmotae to match V4 16S BLAST, remove stray spaces:
for(i in 1:ncol(xlocs))
{
  xlocs[,i]<-gsub("/Shigella", " marmotae", xlocs[,i])
  xlocs2[,i]<-apply(colsplit(xlocs[,i],pattern=" ",names=1:5)[1:2],1,paste,collapse=" ")
}
xlocsm<-melt(xlocs2)
xlocsm2<-as.data.frame(matrix(nrow=nrow(xlocsm),ncol=10))
xlocsm2[,1:ncol(xlocsm)]<-xlocsm
colnames(xlocsm2)[1:ncol(xlocsm)]<-colnames(xlocsm)


# select only ASVs in SFC with >99% identity to closest 16S sequence in BLAST database:
blastax5<-subset(blastax4,V3>99)
# identify ASVs in SFC that match translocators with >99% identity to type strain:
for(i in 1:length(unique(xlocsm[,3])))
{
  taxon<-unique(xlocsm[,3])[i]
  temp1<-blastax5[blastax5 $genspec%in%taxon,]
  if(nrow(temp1)>0)
    for(j in 1:sum(xlocsm[,3]==taxon))
    {
      string<-which(xlocsm[,3]==taxon)
      xlocsm2[string[j],4:(3+length(unique(temp1[,1])))]<-unique(temp1[,1])	
    }	
}



library(ggplot2);library(reshape2)
# read in mucosal 16S sequencing data from IBD cohort subjects:
tax<-read.csv("~/Downloads/fig3_210128_cell2020_16s/taxaspec_md5_fasta_blast.csv",row.names=1)
muctab<-read.table("~/Downloads/fig3_210128_cell2020_16s/seqtab_md5.xls",sep="\t",header=T,row.names=1)
map<-read.csv("~/Downloads/fig3_210128_cell2020_16s/SraRunTable_amp_map.csv",header=T,row.names=1)
realtax<-rownames(tax)[tax$blastax!=""]
muctab2<-muctab[rownames(muctab)%in%realtax,]
mapmuc<-subset(map,muc=="MUC")
muctab3<-muctab2[,colnames(muctab2)%in%rownames(mapmuc)]

# add up reads among samples from same individual (uninvolved + involved mucosal biopsies):
mucsumtab<-matrix(ncol=length(unique(mapmuc$PID)),nrow=nrow(muctab3))
colnames(mucsumtab)<-unique(mapmuc$PID)
rownames(mucsumtab)<-rownames(muctab3)
for(i in 1:length(unique(mapmuc$PID)))
{
  subj<-unique(mapmuc$PID)[i]
  temp<-subset(mapmuc,PID==subj)
  mucsumtab[,i]<-rowSums(muctab3[,colnames(muctab3)%in%rownames(temp)])
}

# select only samples with >19000 reads:
mucsumtab3<-mucsumtab[,colnames(mucsumtab)%in%colnames(mucsumtab)[colSums(mucsumtab)>19000]]
mucsumtab3b<-apply(mucsumtab3,2,function(x) x/sum(x))
colnames(mucsumtab3b)<-gsub("UC0","UC", colnames(mucsumtab3b))
tax$blastax<-gsub("[[]", "", tax$blastax)
tax$blastax<-gsub("[]]", "", tax$blastax)
tax2<-colsplit(tax$blastax,pattern=" ",1:10)
tax3<-cbind(tax,genspec=paste(tax2[,1],tax2[,2], sep=" "))



# block produces list of all translocators in people who also had 16S data:
allxlocsin16s<-xlocs2[,colnames(xlocs2)%in% colnames(mucsumtab3b)]
allxlocs<-unique(melt(allxlocsin16s)[,3])
allxlocs<-allxlocs[!allxlocs%in%c(" ","NA NA")]
allxlocs[allxlocs%in%"Klebsiella "]<-blastax3$genspec[grep("Klebsiella",blastax3$genspec)[1]]
sumtabxloc<-mucsumtab3b[,colnames(mucsumtab3b)%in%colnames(allxlocsin16s)]


# read in SFC-IgG-seq data on IBD cohort:
m4<-read.csv("~/Downloads/fig3_IgG-scores.csv",row.names=1)


# create list of all possible translocator ASV IDs to examine:
taxlist<-unique(unlist(xlocsm2[,4:ncol(xlocsm2)]))
taxlist<-taxlist[!is.na(taxlist)]
# remove internal spike-in control ASV:
taxlist<-taxlist[!taxlist%in%"f46a7ca244afef522b22a11bd33d27b1"]
resmat<-matrix(nrow=length(taxlist),ncol=10)
resmat[,1]<-taxlist
resmat<-as.data.frame(resmat)
colnames(resmat)<-c("tax","meanwithoverwo","numwith","numwo","zhealthy","znoxloc","zxloc","lmp")
# loop that iteratively goes through each ASV in the translocator sequences data:
for(k in 1:length(taxlist))
{
  temp<-apply(xlocsm2, 2, function(x) grep(taxlist[k],x))
  xlocsm2$temp<-NA
  xlocsm2$temp[unlist(temp)]<-TRUE	
  subjswithxloc<-unique(xlocsm2[unlist(temp),2])
  # remove non-mucosal sample "CD11MLN":
  subjswithxloc<-subjswithxloc[!subjswithxloc%in%"CD11MLN"]
  paste(subjswithxloc,collapse=" ")
  # select all remaining subjects subjected to culture (i.e. subjects for whom bacterium was not found to have translocated):
  subjswithout<-unique(xlocsm2[,2])[!unique(xlocsm2[,2])%in%subjswithxloc]
  subjswithout<-subjswithout[!subjswithout%in%"CD11MLN"]
  paste(subjswithout,collapse=" ")
  taxgenspec<-xlocsm2[unlist(temp)[1],3]
  m4a<-m4
  # check if ASV is in both IgG score data and translocator sequences data:
  if(sum(colnames(m4a)%in%paste("X",taxlist[k],sep=""))>0)
  {
    # check if ASV is in both translocator sequences data and mucosal 16S data from same subjects:
    if(taxlist[k]%in%rownames(sumtabxloc))
    {
      # select only taxa for which over 4 individuals exhibited translocation:
      if(length(subjswithxloc)>4)
      {
        colnames(m4a)[colnames(m4a)%in%paste("X",taxlist[k],sep="")]<-"taxon"
        m4a$xlocstatus<-as.vector("")
        m4a$xlocstatus<-as.vector(m4a$xlocstatus)
        m4a[m4a[,"XsampIDnodup"]%in%subjswithxloc,"xlocstatus"]<-"xloc"
        m4a[m4a[,"XsampIDnodup"]%in% subjswithout,"xlocstatus"]<-"no_xloc"
        m4a$muc16s<-as.vector("")
        m4a$muc16s <-as.vector(m4a$muc16s)
        m4a[m4a[,"XsampIDnodup"]%in%colnames(sumtabxloc),"muc16s"]<-"not present"
        m4a[m4a[,"XsampIDnodup"]%in%colnames(sumtabxloc)[sumtabxloc[taxlist[k],]>0],"muc16s"]<-"present"
        # calculate ratio of mean IgG score in subjects with translocation over mean IgG score in subjects without translocation:
        resmat[k,2]<-mean(as.numeric(m4a[m4a[,"XsampIDnodup"]%in%subjswithxloc,"taxon"]))/mean(as.numeric(m4a[m4a[,"XsampIDnodup"]%in%subjswithout,"taxon"]))
        resmat[k,3]<-sum(subjswithxloc%in%m4a[,"XsampIDnodup"])
        resmat[k,4]<-sum(subjswithout%in%m4a[,"XsampIDnodup"])
        resmat[k,1]<-paste(taxgenspec,"_",substr(taxlist[k],0,4),sep="")
        m4a$xlocstatus[m4a$xlocstatus%in%""]<-"Healthy"
        m4a$xlocstatus[m4a$xlocstatus%in%"xloc"]<-"IBD, translocation +"
        m4a$xlocstatus[m4a$xlocstatus%in%"no_xloc"]<-"IBD, translocation -"
        m4a[,"taxon"]<-scale(as.numeric(as.vector(m4a[,"taxon"])))
        # mean IgG scores per subject grouping:
        resmat[k,5]<-mean(as.numeric(m4a[m4a[,"XIBD_stat"]=="Healthy","taxon"]))
        resmat[k,6]<-mean(as.numeric(m4a[m4a[,"XsampIDnodup"]%in% subjswithout,"taxon"]))
        resmat[k,7]<-mean(as.numeric(m4a[m4a[,"XsampIDnodup"]%in%subjswithxloc,"taxon"]))
        # set ordinal variable of non-IBD, IBD with no translocation, IBD with translocation (in that order):
        m4a$IBD_statnum<-NA
        m4a$IBD_statnum[m4a[,"XIBD_stat"]=="Healthy"]<-1
        m4a$IBD_statnum[m4a[,"XsampIDnodup"]%in% subjswithout]<-2
        m4a$IBD_statnum[m4a[,"XsampIDnodup"]%in%subjswithxloc]<-3
        # linear model testing step-wise difference in IgG score from non-IBD, to IBD with no translocation, to IBD with translocation for the given ASV tested in this iteration of loop:
        resmat[k,8]<-summary(lm(scale(as.numeric(taxon))~as.numeric(IBD_statnum),data=m4a))$coefficients[2,4]
        }
    }
  }
}
resmat2<-resmat[!is.na(resmat[,2]),]
resmat3<-melt((resmat2),id=colnames(resmat)[!colnames(resmat)%in%c("zhealthy","znoxloc","zxloc")])
resmat3$tax <- factor(resmat3$tax, levels = unique(resmat3$tax[order(resmat3$lmp,decreasing=T)]))
resmat3$variable<-gsub("zhealthy","Healthy",resmat3$variable)
resmat3$variable<-gsub("znoxloc","IBD, translocation -",resmat3$variable)
resmat3$variable<-gsub("zxloc","IBD, translocation +",resmat3$variable)
ggplot(resmat3,aes(y=tax,x=variable,fill=as.numeric(as.vector(value))))+geom_tile()+scale_fill_distiller(palette="Blues",direction=1)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("")+ylab("")

resmat3 $tax <- factor(resmat3 $tax, levels =rev(levels(resmat3 $tax)))
ggplot(resmat3,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.65,fatten=5,outlier.shape=NA)+geom_line(aes(group=tax),color="black",size=2)+geom_line(aes(group=tax,colour=tax),size=1.5)+geom_point()+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("mean z-IgG score")+xlab("")+scale_color_viridis(discrete=T,option="C",direction=-1) +geom_boxplot(width=0.65,fatten=5,alpha=0.8,outlier.shape=NA)



