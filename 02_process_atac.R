library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project2/gilad/daraujo/scRNA_scATAC')

# check if compiled file already exists. run code if not.
if ('compiled_atac.rds' %in% list.files() == F){
  # get ATAC peaks from 10x
  atac.peaks <- Read10X_h5('/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/filtered_feature_bc_matrix.h5')
  atac.peaks <- atac.peaks$Peaks 
  # get ATAC cell IDs
  atac.cells <- CountFragments('/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/atac_fragments.tsv.gz') %>% 
    filter(frequency_count>2000) %>% select(CB) %>% pull()
  # load all fragments
  atac.frags <- CreateFragmentObject(path='/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/atac_fragments.tsv.gz',
                                   cells=atac.cells)
  # create featurematrix object
  counts <- FeatureMatrix(fragments=atac.frags, features=GRanges(atac.peaks@Dimnames[[1]]), cells=atac.cells)
  # get cell IDs with donors assigned by vireo
  cells_to_keep <- fread('/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/vireo/donor_ids.tsv.gz') %>% 
    filter(donor_id %in% c('doublet','unassigned') == F) %>% select(cell) %>% pull()
  # create ATAC Seurat object
  atac.assay <- CreateChromatinAssay(counts=counts, min.features=1000, fragments=atac.frags) %>%
    CreateSeuratObject(project='1_1', assay='ATAC') %>% subset(cells=cells_to_keep) %>% RenameCells(add.cell.id='1_1')

  # combine all batches/lanes
    for (j in seq(2,6)){
      # get ATAC peaks from 10x
      atac.peaks <- Read10X_h5('/project/gilad/kenneth/caQTL/highpass/batch1/lane'%&%j%&%'/outs/filtered_feature_bc_matrix.h5')
      atac.peaks <- atac.peaks$Peaks 
      # get ATAC cell IDs
      atac.cells <- CountFragments('/project/gilad/kenneth/caQTL/highpass/batch1/lane'%&%j%&%'/outs/atac_fragments.tsv.gz') %>% 
        filter(frequency_count>2000) %>% select(CB) %>% pull()
      # load all fragments
      atac.frags <- CreateFragmentObject(path='/project/gilad/kenneth/caQTL/highpass/batch1/lane'%&%j%&%'/outs/atac_fragments.tsv.gz',
                                         cells=atac.cells)
      # create featurematrix object
      counts <- FeatureMatrix(fragments=atac.frags, features=GRanges(atac.peaks@Dimnames[[1]]), cells=atac.cells)
      # get cell IDs with donors assigned by vireo
      cells_to_keep <- fread('/project/gilad/kenneth/caQTL/highpass/batch1/lane'%&%j%&%'/outs/vireo/donor_ids.tsv.gz') %>% 
        filter(donor_id %in% c('doublet','unassigned') == F) %>% select(cell) %>% pull()
      # create ATAC Seurat object
      atac.assay.f <- CreateChromatinAssay(counts=counts, min.features=1000, fragments=atac.frags) %>%
        CreateSeuratObject(project='1_'%&%j, assay='ATAC') %>% subset(cells=cells_to_keep) %>% RenameCells(add.cell.id='1_'%&%j)
      atac.assay <- merge(atac.assay, atac.assay.f)
    }
  for (j in seq(1,7)){
    # get ATAC peaks from 10x
    atac.peaks <- Read10X_h5('/project/gilad/kenneth/caQTL/highpass/batch2/lane'%&%j%&%'/outs/filtered_feature_bc_matrix.h5')
    atac.peaks <- atac.peaks$Peaks 
    # get ATAC cell IDs
    atac.cells <- CountFragments('/project/gilad/kenneth/caQTL/highpass/batch2/lane'%&%j%&%'/outs/atac_fragments.tsv.gz') %>% 
      filter(frequency_count>2000) %>% select(CB) %>% pull()
    # load all fragments
    atac.frags <- CreateFragmentObject(path='/project/gilad/kenneth/caQTL/highpass/batch2/lane'%&%j%&%'/outs/atac_fragments.tsv.gz',
                                       cells=atac.cells)
    # create featurematrix object
    counts <- FeatureMatrix(fragments=atac.frags, features=GRanges(atac.peaks@Dimnames[[1]]), cells=atac.cells)
    # get cell IDs with donors assigned by vireo
    cells_to_keep <- fread('/project/gilad/kenneth/caQTL/highpass/batch2/lane'%&%j%&%'/outs/vireo/donor_ids.tsv.gz') %>% 
      filter(donor_id %in% c('doublet','unassigned') == F) %>% select(cell) %>% pull()
    # create ATAC Seurat object
    atac.assay.f <- CreateChromatinAssay(counts=counts, min.features=1000, fragments=atac.frags) %>%
      CreateSeuratObject(project='2_'%&%j, assay='ATAC') %>% subset(cells=cells_to_keep) %>% RenameCells(add.cell.id='2_'%&%j)
    atac.assay <- merge(atac.assay, atac.assay.f)
  }
  
  rm(atac.assay.f, cells_to_keep, atac.peaks, atac.cells, atac.frags, counts)
  SaveSeuratRds(atac.assay, file='compiled_atac.rds')
}

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# perform standard QC/analysis and save Seurat object
atac.assay <- LoadSeuratRds('compiled_atac.rds') #skip this if running script for the first time
Annotation(atac.assay) <- annotation
atac.assay <- NucleosomeSignal(atac.assay)
atac.assay <- TSSEnrichment(atac.assay)

DensityScatter(atac.assay, x='nCount_ATAC', y='TSS.enrichment', log_x=T, quantiles=T)
ggsave('scatter_scATAC_processed_ncount.tssenrichment.pdf', height=5, width=8)

SaveSeuratRds(atac.assay, file='compiled_atac_processed.rds')