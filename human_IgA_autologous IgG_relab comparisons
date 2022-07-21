#the following was used to investigate relationships between mucosal IgA scores, systemic IgG scores, and fecal relative abundances, with results shown in Figure 1

library(labdsv);library(ggplot2);library(lme4);library(lmerTest);library(plyr);library(stringr);library(RColorBrewer);library(Hmisc);library(viridis);library(vegan)
library(ggbeeswarm);library(ape)

  mock2<-read.csv("~/fig1_rar11k.csv",row.names=1)
  map3<-read.csv("~/fig1_map_forR.csv",row.names=1)
  
  #calculate IgA scores:
  ICmat<-matrix(nrow=nrow(mock2),ncol=1)
  rownames(ICmat)<-rownames(mock2)
  mock<-(mock2+1) #pseudocount of 1
  for(i in 1:length(unique(map3$SubjectID)))
  {
    subjID<-unique(map3$SubjectID)[i]
    mapsub<-subset(map3, SubjectID ==subjID & Ig=="IgA"|  SubjectID ==subjID & Ig=="pre")
    possamp<-rownames(subset(mapsub,posneg=="pos"))
    negsamp<-rownames(subset(mapsub,posneg=="neg"))
    # which ASVs in possamp are 0's?
    zeroesinpos<-(mock[,possamp]==1)
    # which ASVs in negsamp are 0's?	
    zeroesinneg<-(mock[,negsamp]==1)
    # 2 means both are "True" which means both are 0 abundance
    bothzeroes<-zeroesinpos+ zeroesinneg
    # "TRUE" are ASVs that are either 0 or 1, so have abundance in either pos or neg
    bothzeropositions<-(bothzeroes==2)
    mock[bothzeropositions,possamp]<-0
    mock[bothzeropositions,negsamp]<-0
    ICmat<-cbind(ICmat,log(mock[,possamp],10)-log(mock[,negsamp],10))
    colnames(ICmat)[ncol(ICmat)]<-as.character(subjID)
  }
  ICmatIgA<-ICmat
  ICmatIgApruned<-ICmatIgA[apply(ICmatIgA,1,function(x)sum(is.finite(x)))>(.1*length(colnames(ICmatIgA))),]

  
  #calculate IgG scores:
  ICmat<-matrix(nrow=nrow(mock2),ncol=1)
  rownames(ICmat)<-rownames(mock2)
  mock<-(mock2+1)
  for(i in 1:length(unique(map3$SubjectID)))
  {
    subjID<-unique(map3$SubjectID)[i]
    mapsub<-subset(map3, SubjectID ==subjID & Ig=="IgG"|  SubjectID ==subjID & Ig=="pre")
    possamp<-rownames(subset(mapsub,posneg=="pos"))
    negsamp<-rownames(subset(mapsub,posneg=="neg"))
    # which ASVs in possamp are 0's?
    zeroesinpos<-(mock[,possamp]==1)
    # which ASVs in negsamp are 0's?	
    zeroesinneg<-(mock[,negsamp]==1)
    # 2 means both are "True" which means both are 0 abundance
    bothzeroes<-zeroesinpos+ zeroesinneg
    # "TRUE" are ASVs that are either 0 or 1, so have abundance in either pos or neg
    bothzeropositions<-(bothzeroes==2)
    mock[bothzeropositions,possamp]<-0
    mock[bothzeropositions,negsamp]<-0
    ICmat<-cbind(ICmat,log(mock[,possamp],10)-log(mock[,negsamp],10))
    colnames(ICmat)[ncol(ICmat)]<-as.character(subjID)
  }
  ICmatIgG<-ICmat
  ICmatIgGpruned<-ICmatIgG[apply(ICmatIgG,1,function(x)sum(is.finite(x)))>(.1*length(colnames(ICmatIgG))),]
  library(Hmisc)
  ICmatIgG2<-ICmatIgGpruned[,!is.na(colnames(ICmatIgGpruned))]
  IgGandIgA<-merge(ICmatIgApruned, ICmatIgG2,by=0)
  library(ggplot2)
  mat1<-matrix(nrow=length(unique(colnames(ICmatIgG2))),ncol=2)
  igahigher2 <-NA
  igghigher2 <-NA
  
  # iteratively examining each human participant, collate ASV names of ASVs with 10 fold differences in IgG vs IgA score (1 unit in log10 space): 
  for(i in 1:length(unique(colnames(ICmatIgG2)))){
    IgGandIgAcomp<-as.data.frame(	IgGandIgA[,grep(unique(colnames(ICmatIgG2))[i],colnames(IgGandIgA))])
    igghigher2 <-c(igghigher2,IgGandIgA$Row.names[which((IgGandIgAcomp[,1]-IgGandIgAcomp[,2])<(-1))])
    igahigher2 <-c(igahigher2,IgGandIgA$Row.names[which((IgGandIgAcomp[,1]-IgGandIgAcomp[,2])>1)])
  }
igahigher<-igahigher2[2:length(igahigher2)]
igghigher<-igghigher2[2:length(igghigher2)]



# generate relative abundance table by selecting "pre" fractionation samples (stool with no selection of Ig+/- bacteria):
preonly<-mock[,grep("pre",colnames(mock))]
colnames(preonly)<-map3[match(colnames(preonly),rownames(map3)),"SubjectID"]
preonly<-preonly[,!is.na(colnames(preonly))]
preonly2<-preonly[,colnames(preonly)%in%colnames(ICmatIgG2)]
preonly3<-log(preonly2,10)
preonly3<-preonly3[,order(colnames(preonly3))]



# calculate correlation statistics comparing IgG scores to fecal relative abundance:
mpre<-melt(cbind(preonly3,arbtax=rownames(preonly3)))
migg<-melt(ICmatIgG2)
colnames(migg)<-c("arbtax","variable","iggvalue")
migg2<-cbind(migg,arbvar=paste(migg$arbtax,migg$variable))
mpre2<-cbind(mpre,arbvar=paste(mpre$arbtax,mpre$variable))
mboth<-merge(migg2,mpre2,by="arbvar",suffixes=c(".igg",".pre"))
compmat<-matrix(nrow=nrow(preonly3),ncol=3)
colnames(compmat)<-c("otu","rho","pvalue")
compmat[,1]<-rownames(preonly3)
for(i in 1:length(unique(rownames(preonly3))))
{
  temp1<-subset(mboth,arbtax.igg==unique(rownames(preonly3))[i])
  temp<-cbind(as.numeric(as.vector(temp1[,"iggvalue"])),as.numeric(as.vector(temp1[,"value"])))
  temp2<-temp[!is.na(rowSums(temp)),]
  if (length(temp2)[1]>10)
  {
    if(mean(as.numeric(as.vector(temp2[,1])),na.rm=T)>0){
      compmat[i,2] = rcorr( temp2, type="spearman")$r[2,1]
      compmat[i,3] = rcorr(temp, type="spearman")$P[2,1]
    }
  }
}
compmat2<-compmat[!is.na(compmat[,2]),]
colnames(mboth)[2]<-"otu"
mboth2<-merge(mboth,compmat2,by="otu")
mboth3<-mboth2[!is.na(rowSums(mboth2[c("iggvalue","value")])),]
mboth3iggsub<-mboth3



# calculate correlation statistics comparing IgA scores to fecal relative abundance:
mpre<-melt(cbind(preonly3,arbtax=rownames(preonly3)))
miga<-melt(ICmatIgApruned)
colnames(miga)<-c("arbtax","variable","igavalue")
miga2<-cbind(miga,arbvar=paste(miga$arbtax,miga$variable))
mpre2<-cbind(mpre,arbvar=paste(mpre$arbtax,mpre$variable))
mboth<-merge(miga2,mpre2,by="arbvar",suffixes=c(".iga",".pre"))
compmat<-matrix(nrow=nrow(preonly3),ncol=3)
colnames(compmat)<-c("otu","rho","pvalue")
compmat[,1]<-rownames(preonly3)
for(i in 1:length(unique(rownames(preonly3))))
{
  temp1<-subset(mboth,arbtax.iga==unique(rownames(preonly3))[i])
  temp<-cbind(as.numeric(as.vector(temp1[,"igavalue"])),as.numeric(as.vector(temp1[,"value"])))
  temp2<-temp[!is.na(rowSums(temp)),]
  if (length(temp2)[1]>10)
  {
    if(mean(as.numeric(as.vector(temp2[,1])),na.rm=T)>0){
      compmat[i,2] = rcorr( temp2, type="spearman")$r[2,1]
      compmat[i,3] = rcorr(temp, type="spearman")$P[2,1]
    }
  }
}
compmat2<-compmat[!is.na(compmat[,2]),]
colnames(mboth)[2]<-"otu"
mboth2<-merge(mboth,compmat2,by="otu")
mboth3<-mboth2[!is.na(rowSums(mboth2[c("igavalue","value")])),]
mboth3igasub<-mboth3


# merge results:
colnames(mboth3iggsub)<-gsub("igg","ig",colnames(mboth3iggsub))
colnames(mboth3igasub)<-gsub("iga","ig",colnames(mboth3igasub))
bothbo<-rbind(cbind(mboth3iggsub,ig="igg"),cbind(mboth3igasub,ig="iga"))
bothbo2<-bothbo[match(unique(paste(bothbo$otu,bothbo$ig)),paste(bothbo$otu,bothbo$ig)),]
bothbo2$logp<-(-log(as.numeric(as.vector(bothbo2$pvalue)),10))
bothbo3<-bothbo2[bothbo2$otu%in%names(table(bothbo2$otu))[which(table(bothbo2$otu)==2)],]
bothbo4<-bothbo3[bothbo3$otu%in% unique(c(igahigher, igghigher)),]


library(mratios)
bothbo4$ig<-factor(bothbo4$ig,levels=c("iga","igg"))
ggplot(bothbo4,aes(x=ig,y=abs(as.numeric(as.vector(logp)))))+geom_point(aes(colour= ig),alpha=0.4)+geom_line(aes(group=otu),alpha=0.081)+theme_bw()+ scale_color_viridis(end=0.5,na.value = "grey50",discrete=T)+ylab("ASV correlation to fecal\nrelative abundance, P value")+xlab("")+geom_boxplot(alpha=0.4,outlier.shape=NA,width=0.3,fatten=4)+ylim(c(0,3))
ttestratio(abs(as.numeric(logp))~as.factor(ig), data=bothbo4)

ggplot(bothbo4,aes(x=ig,y=abs(as.numeric(as.vector(rho))*10)))+geom_point(aes(colour= ig),alpha=0.4)+geom_line(aes(group=otu),alpha=0.081)+theme_bw()+ scale_color_viridis(end=0.5,na.value = "grey50",discrete=T)+ylab("ASV correlation to fecal\nrelative abundance, rho value")+xlab("")+geom_boxplot(alpha=0.4,outlier.shape=NA,width=0.3,fatten=4)
ttestratio(abs(as.numeric(rho))~as.factor(ig), data=bothbo4)
