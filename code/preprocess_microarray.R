setwd("~/projects/comp-sys-pharm/SARS-CoV-2-cholesterol/code")
library(limma)
library("hgug4112a.db")

for (dname in c('28166','33267','37571','45042','56677')){
  

  dictname=paste0('../data/gse',dname,'/GSE',dname,'_RAW/')
  
  fnames=list.files(dictname)
  fnames=fnames[grep('GSM',fnames)]

  data=read.maimages(path = dictname,files = fnames,source="agilent", 
                     green.only=TRUE, other.columns="gIsWellAboveBG")
  data_bg=backgroundCorrect(data, method="normexp")
  data_norm=normalizeBetweenArrays(data_bg, method="quantile")
  data_norm$genes$SYMBOL=mapIds(hgug4112a.db,data_norm$genes$ProbeName,keytype="PROBEID", column="SYMBOL")
  
  expression=data_norm$E
  rownames(expression)=data_norm$genes$SYMBOL
  data=limma::avereps(expression[!is.na(rownames(expression)),])

  dictname=paste0('../data/gse',dname,'/gse',dname,'_norm_mean.csv')
  write.csv(data,file = dictname)
} 

print("SUCCESS")