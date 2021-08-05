setwd("~/projects/comp-sys-pharm/SARS-CoV-2-cholesterol/code/")

library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
library(readr)
library(progeny)
library('BiocParallel')

n=10
register(MulticoreParam(n))


set.seed(2020)

# Dorothea
data(dorothea_hs, package = "dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

calculate_tf_activities <- function(filename) {
  exp=read.csv(paste0('../data/drug_signatures/signatures_',filename,'.csv'),sep = ',',header = T,row.names = 1, check.names = FALSE)
  
  tf_activities <- run_viper(exp, regulons, 
                             options =  list(method = 'none', minsize = 4, 
                                             eset.filter = FALSE, cores = 1, 
                                             verbose = FALSE))
  tf_activities <- t(tf_activities)
  write.csv(tf_activities,paste0('../results/drug_signatures/signatures_',filename,'_dorothea.csv'))
  
}


# Progeny
calculate_pathway_activities <- function(filename){
  
  exp=read.csv(paste0('../data/drug_signatures/signatures_',filename,'.csv'),sep = ',',header = T,row.names = 1, check.names = FALSE)
  print('READ  IN OK')
  
  pathways_zscore <- t(progeny(as.matrix(exp),
                               scale=TRUE, organism="Human",
                               top = 100, perm = 10000, z_scores = TRUE))
  colnames(pathways_zscore) = colnames(exp)
  pathways_zscore = data.frame(pathways_zscore, check.names = FALSE)
  
  print('Z scores ok')
  write.csv(pathways_zscore,paste0('../results/drug_signatures/signatures_',filename,'_progeny_z.csv'))
  print('Saved')
  
  print('Calculate carnival input')
  if (nrow(exp) > 1000){
    
      for(i in seq(1, ncol(exp), 1)) {

        if(i %% 300 != 0 && i !=ncol(exp)) {
          next
        }

        length = if(i==ncol(exp)) ncol(exp)%%300 else 300
        print(i+1-length)
        print(i)
        exp_part = exp[,(i+1-length):i]
        
        pathways_inputCarnival_part <-
          t(progeny(as.matrix(exp_part),
                    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE))
        
        colnames(pathways_inputCarnival_part) = colnames(exp_part)
        pathways_inputCarnival_part = data.frame(pathways_inputCarnival_part, check.names = FALSE)
        
        if (exists("pathways_inputCarnival")){
          pathways_inputCarnival <- cbind(pathways_inputCarnival, pathways_inputCarnival_part)
          print('next concatenated')
          print(ncol(pathways_inputCarnival))
        } else {
          print('first part saved')
          pathways_inputCarnival <- pathways_inputCarnival_part
          
        }
       }
       
  } else {
    pathways_inputCarnival <- 
      t(progeny(as.matrix(exp), 
                scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE))
    
    colnames(pathways_inputCarnival) = colnames(exp)
    pathways_inputCarnival = data.frame(pathways_inputCarnival, check.names = FALSE)
    
    print('progeny ok')    
  }

  write.csv(pathways_inputCarnival,paste0('../results/drug_signatures/signatures_',filename,'_progeny_carnival.csv'))
  print('saved')
  rm(pathways_inputCarnival)
}



calculate_tf_activities('lm_gene')
calculate_tf_activities('bing_gene')
calculate_pathway_activities('lm_gene')
calculate_pathway_activities('bing_gene')