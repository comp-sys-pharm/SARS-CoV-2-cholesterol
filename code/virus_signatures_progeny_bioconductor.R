library(progeny)

# Progeny pathway activity analysis
setwd("~/projects/comp-sys-pharm/SARS-CoV-2-cholesterol/code/")


set.seed(2020)
# Read in data
data_exp=read.csv('../results/microarray_and_rnaseq_expression.csv',
                  sep = ',',header = T,row.names = 1, check.names = FALSE)


pathways_zscore <- t(progeny(as.matrix(data_exp), 
                             scale=TRUE, organism="Human", 
                             top = 100, perm = 10000, z_scores = TRUE))
rownames(pathways_zscore)[rownames(pathways_zscore) == 'JAK.STAT'] <- 'JAK_STAT'

write.csv(pathways_zscore,'../results/microarray_and_rnaseq_progeny_z_bioconductor.csv')


pathways_inputCarnival <- 
  t(progeny(as.matrix(data_exp), 
            scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE))

write.csv(pathways_inputCarnival,'../results/microarray_and_rnaseq_progeny_bioconductor_inputCarnival.csv')

