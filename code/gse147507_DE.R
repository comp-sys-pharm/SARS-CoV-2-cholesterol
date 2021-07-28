setwd("~/projects/comp-sys-pharm/SARS-CoV-2-cholesterol/code/")

library('BiocParallel')
library('DESeq2')
library('tidyverse')
library('stringr')

n=10
register(MulticoreParam(n))

# read in data
file = '../data/gse147507/GSE147507_raw.csv'
data=read.csv(file,sep=',',header=T,row.names = 1)
meta=read.csv('../data/gse147507/meta.csv',sep=',',header=T,row.names = 1)
meta=meta[colnames(data),]

#remove 0 genes
fil=apply(data,1,sum)!=0
data=data[fil,]

# Columns as factors
meta[,2:ncol(meta)] <- lapply(meta[,2:ncol(meta)], as.factor)

# Variance Stabilizing Transformation 
dds=DESeqDataSetFromMatrix(countData = data,
                           colData = meta,
                           design = ~ Cell)
dds=DESeq(dds,parallel = T)
vsd <- vst(dds, blind=TRUE)
write.csv(assay(vsd),paste0('../results/gse147507/gse147507_vst.csv'))

# Function of differential expression analysis
# Params: data: RNA seq raw data 
# Params: meta: dataframe of selected conditions (Treatment)
# Params: name: string, save file as 'meta+series_number_DE.csv'
de_analysis <- function(data, meta, name){
  data=data[,rownames(meta)]
  fil=apply(data,1,sum)!=0
  data=data[fil,]
  meta$Virus=as.factor(meta$Virus)
  dds=DESeqDataSetFromMatrix(countData = data,
                             colData = meta,
                             design = ~ Virus)
  dds=DESeq(dds,parallel = T)
  write.csv(results(dds),paste0('../data/gse147507/',name,'_DE.csv'))
}

for (series in unique(meta$Series)){
  print(series)
  series_number = as.numeric(substr(series, 7, nchar(series)))
  if (series_number == 8){
    # In Series 8 more virus were tested (RSV and HPIV3)
    meta1 = meta[(meta$Series == series) & ((meta$Treatment == 'Mock')| (meta$Treatment == 'RSV')),]
    meta2 = meta[(meta$Series == series) & ((meta$Treatment == 'Mock')| (meta$Treatment == 'HPIV3')),]
    print(meta1)
    de_analysis(data, meta1, paste0('series',series_number,'a'))
    print(meta2)
    de_analysis(data, meta2, paste0('series',series_number,'b'))
    #write.csv(data1, paste0('../data/',accession,'/meta',series_number',a.csv'))
    #write.csv(data2, paste0('../data/',accession,'/meta',series_number',b.csv'))
  } else if (series_number == 9){
    # In Series 9 more virus were tested (IAV selected)
    meta1 = meta[(meta$Series == series) & ((meta$Treatment == 'Mock')| (meta$Treatment == 'IAV')),]
    print(meta1)
    de_analysis(data, meta1, paste0('series',series_number))
  } else {
    meta1 = meta[meta$Series==series, ]
    print(meta1)
    #write.csv(data1, paste0('../data/',accession,'/meta',series_number',.csv'))
    de_analysis(data, meta1, paste0('series',series_number))
  }
}