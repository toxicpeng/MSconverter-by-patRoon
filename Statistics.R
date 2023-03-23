rm(list = ls())

library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggrepel)

####set working path
path_data<-'E:/48NRs project/Data_analysis'
setwd(path_data)

####read data and set target protein
polarity<--1 
Rawdata<-read.csv("48NRs_2_Neg.csv")
Batch<-"Batch2"
sampleID<-read_excel("SampleID.xlsx",sheet=Batch)
Targetp<-"PXR"
Targetex<-"Dust"

####Calculate fc and pvalue
####Dust extract
mydata<-as_tibble(Rawdata)

mydata$fold1<-rep(1,nrow(mydata))
mydata$fold2<-rep(1,nrow(mydata))
mydata$fold3<-rep(1,nrow(mydata))

mydata$pvalue1<-rep(0,nrow(mydata))
mydata$pvalue2<-rep(0,nrow(mydata))
mydata$pvalue3<-rep(0,nrow(mydata))

temp<-sampleID[which(sampleID$Protein==Targetp),]
temp_WT<-sampleID[which(sampleID$Protein=="WT"),]

index_test<-temp[which(temp$ExtractID==Targetex),]
index_crtl1<-temp_WT[which(temp$ExtractID==Targetex),]
index_crtl2<-temp[which(temp$ExtractID=="DMSO"),]
index_crtl3<-temp_WT[which(temp$ExtractID=="DMSO"),]

FileID_test<-index_test$ID
FileID_crtl1<-index_crtl1$ID
FileID_crtl2<-index_crtl2$ID
FileID_crtl3<-index_crtl3$ID

testdata<-mydata[, colnames(mydata)%in%FileID_test]
data.crtl1<-mydata[, colnames(mydata)%in%FileID_crtl1]
data.crtl2<-mydata[, colnames(mydata)%in%FileID_crtl2]
data.crtl3<-mydata[, colnames(mydata)%in%FileID_crtl3]

for (i in 1:nrow(mydata)){

    print(i)
  
    test<-testdata[i,]
    ctrl1<-data.crtl1[i,]
    ctrl2<-data.crtl2[i,]
    ctrl3<-data.crtl3[i,]
    
    fold1<-(sum(test)*3)/(sum(ctrl1)*3)
    fold2<-(sum(test)*3)/(sum(ctrl2)*3)
    fold3<-(sum(test)*3)/(sum(ctrl3)*3)
    
    mydata$fold1[i]<-fold1
    mydata$fold2[i]<-fold2
    mydata$fold3[i]<-fold3
   
    if (sd(test)>0){
      ttest1<-t.test(test,ctrl1)
      pvalue1<-ttest1$p.value
      mydata$pvalue1[i]<-pvalue1
    } else {mydata$pvalue1[i]<-1}
    
    if (sd(test)>0){
      ttest2<-t.test(test,ctrl2)
      pvalue2<-ttest2$p.value
      mydata$pvalue2[i]<-pvalue2
    } else {mydata$pvalue2[i]<-1}
    
    if (sd(test)>0){
      ttest3<-t.test(test,ctrl3)
      pvalue3<-ttest3$p.value
      mydata$pvalue3[i]<-pvalue3
    } else {mydata$pvalue3[i]<-1}   

  }

mydata<-as.data.frame(lapply(mydata,as.numeric))

####filter out significantly increased features
sig_temp<-mydata%>%filter(fold1>5&pvalue1<0.05)%>%
  filter(fold2>5&pvalue2<0.05)%>%
  filter(fold3>5&pvalue3<0.05)
ID_Targetex<-sampleID[which(sampleID$Group==Targetex),]
index_Targetex<-which(colnames(sig_temp)%in%ID_Targetex$ID)
sig<-filter(sig_temp, sig_temp[,index_Targetex[[1]]]!=100 & sig_temp[,index_Targetex[[2]]]!=100)

####calculate average intensity
sig$avgex<-rep(1,nrow(sig))
sig$avgsam<-rep(1,nrow(sig))

data_avgsam<-sig[, colnames(sig)%in%FileID_test]

index_avgex<-sampleID[which(sampleID$Group==Targetex),]
FileID_avgex<-index_avgex$ID
data_avgex<-sig[, colnames(sig)%in%FileID_avgex]

for (j in 1:nrow(sig)){
     
     print(j)
     print(FileID_avgex)
     print(FileID_test)

     ex<-data_avgex[j,]
     avgex<-sum(ex)/2
     sig$avgex[j]<-avgex
     
     sam<-data_avgsam[j,]
     avgsam<-sum(sam)/3
     sig$avgsam[j]<-avgsam

}

####save the sig data
write.table(sig, file='PXR_Dust_Neg.csv', sep=',', row.names = FALSE)

####prepare dataset
dataset<-mydata

dataset$avgex<-rep(100,nrow(dataset))
dataset$avgsam<-rep(100,nrow(dataset))
dataset$change<-rep("stable",nrow(dataset))
sig$change<-rep("up",nrow(sig))

dataset<-dataset[-which(dataset$pvalue2<0.05&dataset$fold2<0.2),]
dataset<-dataset[-which(dataset$pvalue2<0.05&dataset$fold2>5),]
dataset<-dataset[-which(dataset$pvalue1<0.05&dataset$fold1>5),]
dataset<-dataset[-which(dataset$pvalue3<0.05&dataset$fold3>5),]

dataset<-bind_rows(dataset, sig)

####volcano plot
cut_off_pvalue = -log10(0.05)
cut_off_logFC = log10(5)
temp_p<-arrange(sig,desc(avgsam))

p <- ggplot(
  dataset, aes(x = log10(fold1), y = -log10(pvalue1), colour=change, size=avgsam)) +
  scale_color_manual(values=c("#231f20", "#be1e2d")) +
  geom_point(alpha = 0.7, shape = 16) +
  geom_hline(yintercept = cut_off_pvalue, col="#be1e2d", linetype = 2) +
  geom_vline(xintercept = cut_off_logFC, col="#be1e2d",linetype = 2) +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA)) +
  theme(panel.grid.minor=element_line(colour=NA)) +
  geom_text_repel(
    data=temp_p[1:5,],
    aes(label=round(mz,digits=4)),
    size=4,
    color="#be1e2d",
    segment.color="black",
    show.legend = FALSE) +
  theme(legend.position = "none") +
  xlab("Log(Fold change)") +
  ylab("-Log P value") +
  theme(axis.text.x = element_text(size = rel(1.5), color = "black")) +
  theme(axis.text.y = element_text(size = rel(1.5), color = "black")) +
  theme(axis.title.x = element_text(size = rel(1.2), color = "black")) +
  theme(axis.title.y = element_text(size = rel(1.2), color = "black"))
  
p
 
ggsave(p, filename = "PXR pulldown_Neg.tiff", width = 6, height = 5) 







