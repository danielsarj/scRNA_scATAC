library(Signac)
library(Seurat)
library(data.table)
library(argparse)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project2/gilad/daraujo/scRNA_scATAC')
parser <- ArgumentParser()
parser$add_argument('-c')
args <- parser$parse_args()

# load dataset
assay <- LoadSeuratRds('compiled_atac_processed.rds')

# find DA peaks
da_peaks <- FindMarkers(object=assay, ident.1=as.character(args$c), test.use='LR', 
                           latent.vars='nCount_ATAC')

# save result
fwrite(da_peaks, file='DA_peaks_c' %&% as.character(args$c) %&% '.txt',
      sep=' ', row.names=T, col.names=T)