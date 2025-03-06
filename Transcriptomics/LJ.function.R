suppressMessages({
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(tibble)
  
  library(miloR)
  library(SingleCellExperiment)
  
  library(scater)
  library(scran)
  library(igraph)
  
  library(batchelor)
  library(uwot)
  library(Seurat)
  library(DoubletFinder)
  library(GGally)
  library(viridis)
})

########### TPM function ##############################################################################
count2tpm <- function(countMat, idType = "Ensembl"){                                                 ##
  if(!sum(class(countMat)%in%"matrix"))  countMat<-as.matrix(countMat)                               ##
  ensembl<-read.csv(paste0(Gpath,"Ensembl_length.csv"))                                              ##
  if(toupper(idType) == "ENSEMBL"){                                                                  ##
    len <- ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), "Length"]                      ##
  }else if(toupper(idType) == "SYMBOL"){                                                             ##
    len <- ensembl[match(rownames(countMat), ensembl$hgnc_symbol), "Length"]                         ##
  }else if(toupper(idType) == "ENTREZ"){                                                             ##
    len <- ensembl[match(rownames(countMat), ensembl$entrezgene_id), "Length"]                       ##
  }else{                                                                                             ##
    stop("Please input right type of gene name, such as Ensembl or gene Symbol ...")                 ##
  }                                                                                                  ##
  na_idx <- which(is.na(len))                                                                        ##
  if(length(na_idx)>0){                                                                              ##
    warning(paste0("Omit ", length(na_idx), " genes of which length is not available !"))            ##
    cat(paste0("There are ", length(na_idx), " genes that did not match any gene length.\n"))        ##
    countMat <- countMat[!is.na(len),]                                                               ##
    len = len[!is.na(len)]                                                                           ##
  }                                                                                                  ##
  tmp <- countMat / len                                                                              ##
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))                                                              ##
  TPM <- TPM[!duplicated(rownames(TPM)),]                                                            ##
  TPM <- TPM[!is.na(rownames(TPM)),]                                                                 ##
  TPM <- TPM[!rownames(TPM)==" ",]                                                                   ##
  TPM <- TPM[,!is.na(colnames(TPM))]                                                                 ##
  TPM <- TPM[,!colnames(TPM)==" "]                                                                   ##
  return(TPM)                                                                                        ##
}                                                                                                    ##
#######################################################################################################
