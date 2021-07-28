# Differential expression analysis of GSE148729
# SARS-CoV (S1) and SARS-CoV2 (S2) and Mock as control

setwd("~/projects/comp-sys-pharm/SARS-CoV-2-cholesterol/code/")

library('BiocParallel')
library('DESeq2')
library('tidyverse')
library('stringr')

n=10
register(MulticoreParam(n))

#read data
data=read.csv('../data/gse148729/GSE148729_raw.csv',sep=',',header=T,row.names = 1)
meta=read.csv('../data/gse148729/meta.csv',sep=',',header=T,row.names = 1)
meta=meta[colnames(data),]

#remove all 0 genes
fil=apply(data,1,sum)!=0
data=data[fil,]

# Columns as factors
meta$Treatment <- as.factor(meta$Treatment)
meta$Time <- as.factor(meta$Time)
data = sapply(data, as.integer)

# Variance Stabilizing Transformation 
dds=DESeqDataSetFromMatrix(countData = data,
                           colData = meta,
                           design = ~ Treatment + Time)
dds=DESeq(dds,parallel = T)
vsd <- vst(dds, blind=TRUE)
write.csv(assay(vsd),'../results/gse148729/gse148729_vst.csv')

# Function of differential expression analysis
# Parameters: selected control (mock/untr), virus (S1/S2) and time (4H/24H)
de_analysis <- function(control, virus, time){

  # Filter
  meta_temp=meta[(meta$Treatment==control | meta$Treatment==virus) & meta$Time==time,]
  print(meta_temp)
  
  # New column Treatment_Time as design of DE design formula
  meta_temp$TreatmentTime = paste0(meta_temp$Treatment,'_',meta_temp$Time)
  meta_temp = within(meta_temp, rm(Treatment, Time))
  meta_temp$TreatmentTime <- as.factor(meta_temp$TreatmentTime)
  data_temp=data[,rownames(meta_temp)]

  # DE analysis
  dds=DESeqDataSetFromMatrix(countData = data_temp,
                             colData = meta_temp,
                             design = ~ TreatmentTime)
  dds=DESeq(dds,parallel = T)
  # Save results
  result = results(dds)[!is.na(results(dds)$stat),]
  write.csv(result,paste0('../data/gse148729/gse148729_DE_',control,'_',virus,'_',time,'.csv'))
}

# Analysis
de_analysis("mock", "S1", "24H")
de_analysis("mock", "S2", "24H")
de_analysis("mock", "S1", "4H")
de_analysis("mock", "S2", "4H")


