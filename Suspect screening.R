rm(list = ls())

library(ChemmineR)
library(rcdk)
library(tidyr)
library(dplyr)

####build the search function
InitialSearch<-function(Library,ppm,polarity,Database){
  if (polarity==1){
    Adducts<-c('[M+H]+','[M+H-H2O]+','[M+NH4]+','[M+CH3OH+H]+')
    MW.adducts<-c(1.007825,-17.00274,18.03437,33.03404)-0.00054
  }
  
  if (polarity==-1){
    Adducts<-c('[M-H]-','[M-H-H2O]-','[M+Cl]-','[M+CH2O2-H]-')
    MW.adducts<-c(-1.007825,-19.01839,34.96885,44.99765)+0.00054
  }
  
  indexresults<-NULL
  Library$ID<-rep(0,nrow(Library))
  Library$SMILES<-rep(0,nrow(Library))
  Library$ADDUCT<-rep(0,nrow(Library))
  Library$LOGP<-rep(0,nrow(Library))
  for (i in 1:nrow(Library)){
    mz<-Library$mz[i]
    for (j in 1:length(MW.adducts)){
      mserror<-(mz-Database$MONOISOTOPIC_MASS-MW.adducts[j])/mz
      index1<-which(abs(mserror)<(ppm*10^(-6)))
      if (length(index1)>0){#save the id with match
        indexresults<-c(indexresults,i)
        if (Library$ID[i]==0){
          Library$ID[i]<-paste(c(Database$Order[index1]),collapse = ';')
        }else{
          Library$ID[i]<-paste(c(Library$ID[i],Database$Order[index1]),collapse = ';')
        }
        if (Library$SMILES[i]==0){
          Library$SMILES[i]<-paste(c(Database$SMILES[index1]),collapse = ';')
        }else{
          Library$SMILES[i]<-paste(c(Library$SMILES[i],Database$SMILES[index1]),collapse = ';')
        }
        Adduct.paste<-rep(Adducts[j],length(index1))
        if (Library$ADDUCT[i]==0){
          Library$ADDUCT[i]<-paste(Adduct.paste,collapse = ';')
        }else{
          Library$ADDUCT[i]<-paste(c(Library$ADDUCT[i],Adduct.paste),collapse = ';')
        }
        for (k in 1:length(index1)){
          if (Library$LOGP[i]==0){
            mol<-Database$SMILES[index1[k]]
            mol<-parse.smiles(mol)[[1]]
            convert.implicit.to.explicit(mol)
            Library$LOGP[i]<-get.xlogp(mol)
          }else{
            mol<-Database$SMILES[index1[k]]
            mol<-parse.smiles(mol)[[1]]
            convert.implicit.to.explicit(mol)
            Library$LOGP[i]<-paste(c(Library$LOGP[i],get.xlogp(mol)),collapse = ';')
          }
          
        }
      }
    }##findout matching, ppm
  }
  return(Library[indexresults,])
}

####set working path
path_data <- 'E:/48NRs project/Data_analysis'
setwd(path_data)

####load databases
####Remember to change "MONOISOTOPIC_MASS" and "SMILES" as these will be used by the function

sig <- read.csv("SigResult/Sig_PXR_Dust_Neg_filtered.csv")
sig <- sig %>% select(mz, rt, fold1, fold2, pvalue1, pvalue2, avgex, avgsam)

####load databases and set parameters
Database_TSCA <- read.csv("Database/TSCAMS_edited.csv")
Database_Norman <- read.csv("Database/NORMAN_edited.csv")

####Set the parameters
polarity<--1
Library<-sig
ppm<-3

####Search against the databases
MatchRes_tsca <- InitialSearch(sig, ppm, polarity, Database_TSCA)
MatchRes_nor <- InitialSearch(sig, ppm, polarity, Database_Norman)

####Drop duplicates and label the source
MatchRes_tsca <- MatchRes_tsca[!duplicated(MatchRes_tsca), ]
MatchRes_tsca$database <- "TSCA"
MatchRes_tsca$ID <- NULL
MatchRes_nor <- MatchRes_nor[!duplicated(MatchRes_nor), ]
MatchRes_nor$database <- "Norman"
MatchRes_nor$ID <- NULL

####Combine them together and save
MatchRes <- rbind(MatchRes_tsca, MatchRes_nor)

####Drop duplicate matches from different databases
MatchRes <- separate_rows(MatchRes, SMILES, ADDUCT, LOGP, sep = ";")

####Save the final data
write.table(MatchRes, file = "MatchResults/MatchRes_PXR_Dust_neg.csv", sep = ",", row.names = FALSE)

