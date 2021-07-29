library(limma)
setwd("~/projects/comp-sys-pharm/SARS-CoV-2-cholesterol/code/")

for (dname in c('28166','33267','37571','45042','56677')){
  print(dname)
  fname=paste0('../data/gse',dname,'/gse',dname,'_gex.csv')
  mname=paste0('../data/gse',dname,'/gse',dname,'_meta.csv')
  
  data=read.csv(fname,sep=',',header = T,row.names = 1)
  meta=read.csv(mname,sep=',',header=T, row.names = 1)
  
  meta$Virus = as.factor(meta$Virus)
  meta$Virus=relevel(meta$Virus,ref='Mock')

  meta$FINAL_ID <- gsub(" ", "", paste(meta$Virus,"_",meta$Time))
  
  
  f=meta$FINAL_ID
  design=model.matrix(~0 + f)
  fit <- lmFit(data, design)

  cont=paste0('fGSE',dname,'_Virus_24H-fGSE',dname,'_Virus_0H+fMock_0H-fMock_24H')

  if (dname=='56677'){
    cont=paste0('fGSE',dname,'_Virus_18H-fGSE',dname,'_Virus_0H+fMock_0H-fMock_18H')
  }
  contrast=makeContrasts(cont,levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  fname=paste0('../data/gse',dname,'/gse',dname,'_DE.csv')
  write.csv(topTable(fit2, adjust="BH",number = 1000000),fname)
}
