library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(patRoon)
library(gtools)
library(MSnbase)

options(patRoon.path.SIRIUS = "C:/Program Files/sirius") # location where SIRIUS was extracted
options(patRoon.path.pngquant = "D:/Software/pngquant") # directory containing pngquant binary
options(patRoon.path.MetFragCL = "D:/Software/MetFragCommandLine-2.5.0.jar") # full location to the jar file
options(patRoon.path.MetFragPubChemLite = "D:/Software/PubChemLite_exposomics_20220429.csv") # full location to desired PubChemLite CSV file
options(patRoon.path.MetFragPubChemLite = "D:/Software/PubChem_OECDPFAS_largerPFASparts_20220324.csv") # full location to PFAS DB (NOTE: configured like PubChemLite)
options(patRoon.path.BioTransformer = "D:/Software/biotransformer-3.0.0.jar")
options(patRoon.path.OpenMS = "C:/Program Files/OpenMS-2.7.0/bin")

patRoon::verifyDependencies()

rm(list = ls())
datapath <- "E:/NTA workflow/data"
filters <- c("msLevel [1,1]","polarity -")
outpath <- "E:/NTA workflow/negative"

setwd(datapath)
convertMSFiles(files = datapath, 
               algorithm = "pwiz", 
               from = "thermo", 
               to = "mzXML", 
               outPath = outpath, 
               filters = filters)
