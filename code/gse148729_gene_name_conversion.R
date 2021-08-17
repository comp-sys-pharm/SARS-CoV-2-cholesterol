
setwd("~/projects/comp-sys-pharm/SARS-CoV-2-cholesterol/code/")

library(biomaRt)

OUTFILE='../data/gse148729/gse148729_name_conversion.csv'
data=read.csv('../data/gse148729/GSE148729_Calu3_totalRNA_readcounts.txt',
              header=T,stringsAsFactors=FALSE,sep='\t',row.names = 1)
genelist=sapply(strsplit(rownames(data), "[.]"), head, 1)

#possibilities: "hsapiens_gene_ensembl", ="mmusculus_gene_ensembl", "rnorvegicus_gene_ensembl"
species='hsapiens_gene_ensembl'
#possibilities: 'hgnc_symbol', "mgi_symbol","rgd_symbol"
#possibilities: 'entrezgene','ensembl_gene_id','uniprotswissprot','uniprotsptrembl'
from='ensembl_gene_id'
to='hgnc_symbol'

mart = useMart("ensembl", dataset = species)

genelist=getBM(values=genelist,attributes = c(from,to), 
                filters = from,mart = mart)

write.csv(genelist,file=OUTFILE)
