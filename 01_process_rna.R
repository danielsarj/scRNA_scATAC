library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project2/gilad/daraujo/scRNA_scATAC')
  
# check if compiled file already exists. run code if not.
if ('compiled_rna.rds' %in% list.files() == F){
  # load the RNA data, remove unassigned cells based in vireo's results,
  # remove the scATAC-seq layer, add batch/lane IDs
  cells_to_keep <- fread('/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/vireo/donor_ids.tsv.gz') %>% 
    filter(donor_id %in% c('doublet','unassigned') == F) %>% select(cell) %>% pull()
  counts <- Read10X_h5('/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/filtered_feature_bc_matrix.h5') %>%
    CreateSeuratObject(project='1_1') %>% subset(cells=cells_to_keep) %>% RenameCells(add.cell.id='1_1')
  counts@assays$RNA$counts.Peaks <- NULL
  
  # combine all batches/lanes
  for (j in seq(2,6)){
    cells_to_keep <- fread('/project/gilad/kenneth/caQTL/highpass/batch1/lane'%&%j%&%'/outs/vireo/donor_ids.tsv.gz') %>% 
      filter(donor_id %in% c('doublet','unassigned') == F) %>% select(cell) %>% pull()
    f <- Read10X_h5('/project/gilad/kenneth/caQTL/highpass/batch1/lane'%&%j%&%'/outs/filtered_feature_bc_matrix.h5') %>%
      CreateSeuratObject(project='1_'%&%j) %>% subset(cells=cells_to_keep) %>% RenameCells(add.cell.id='1_'%&%j)
    f@assays$RNA$counts.Peaks <- NULL
    counts <- merge(counts, f)
  }
  for (j in seq(1,7)){
    cells_to_keep <- fread('/project/gilad/kenneth/caQTL/highpass/batch2/lane'%&%j%&%'/outs/vireo/donor_ids.tsv.gz') %>% 
      filter(donor_id %in% c('doublet','unassigned') == F) %>% select(cell) %>% pull()
    f <- Read10X_h5('/project/gilad/kenneth/caQTL/highpass/batch2/lane'%&%j%&%'/outs/filtered_feature_bc_matrix.h5') %>%
      CreateSeuratObject(project='2_'%&%j) %>% subset(cells=cells_to_keep) %>% RenameCells(add.cell.id='2_'%&%j)
    f@assays$RNA$counts.Peaks <- NULL
    counts <- merge(counts, f)
  }
  rm(f, cells_to_keep)
  SaveSeuratRds(counts, file='compiled_rna.rds')
}

# perform standard analysis and save Seurat object
if (exists('counts')==F){
  counts <- LoadSeuratRds('compiled_rna.rds')}
counts <- NormalizeData(counts)
counts <- FindVariableFeatures(counts)
counts <- ScaleData(counts)
counts <- RunPCA(counts)
counts <- RunUMAP(counts, dims=1:30)
SaveSeuratRds(counts, file='compiled_rna_processed.rds')

DimPlot(counts)
ggsave('umap_scRNA_processed_bybatch.lane.pdf', height=5, width=8)
