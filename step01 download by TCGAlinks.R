
#  Screening the potential prognostic long-non coding RNA signatures based on competing endogenous RNA networks in metastasis melanoma###
#  TCGAbiolinks is able to access The National Cancer Institute (NCI) Genomic Data Commons (GDC) thorough its
#  GDC Application Programming Interface (API) to search, download and prepare relevant data for analysis in R.
## TCGAbiolinks main steps: GDCquery;  GDCdownload;  GDCprepare

getwd()
rm(list=ls())
options(stringsAsFactors = F)
#install and load.package
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
#download RNA-seq data
#Genome of reference: hg38  #Accessing GDC.
cancer_type <- "TCGA-SKCM"
query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query, method = "api", files.per.chunk = 100) 
expdat <- GDCprepare(query = query,save = TRUE, save.filename = "exp.rda") #expdat<-load("exp.rda")
RNA_count_matrix <- assay(expdat)
RNA_count_matrix[1:4,1:4]
write.csv(RNA_count_matrix,file = paste(cancer_type,"RNA-seq-Counts.csv",sep = "-"),row.names = F)


#download miRNA-seq data
mi_query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "miRNA Expression Quantification")
GDCdownload(mi_query, method = "api", files.per.chunk = 20)
mi_expdat<- GDCprepare(query = mi_query,save = TRUE, save.filename = "mi_exp.rda")
mi_expdat[1:4,1:4]
del2 <- seq(3, nrow(mi_expdat), by = 3)
del3 <-  seq(3, nrow(mi_expdat), by = 3)
mi_expdat <- mi_expdat[,-c(del2,del3)]
colnames(mi_expdat) <- substr(colnames(mi_expdat),12,27)
write.csv(mi_expdat,file = paste(cancer_type,"miRNA-seq.csv",sep = "-"),row.names = F)

#download clinical data
#Parse XML clinical data
#Genome of reference: hg38
cli_query <- GDCquery(project = "TCGA-SKCM", 
                  data.category = "Clinical", 
                  file.type = "xml" )
GDCdownload(cli_query)
clinical <- GDCprepare_clinic(cli_query, clinical.info = "patient")
nrow(clinical) #493
write.csv(clinical,file = paste(cancer_type,"clinical.csv",sep="-"),row.names = F)



