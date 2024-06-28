library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project2/gilad/daraujo/scRNA_scATAC')

cells_to_keep <- fread('/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/vireo/donor_ids.tsv.gz') %>% 
  filter(donor_id %in% c('doublet','unassigned') == F) %>% select(cell) %>% pull()

frags <- fread('zgrep -v "^#" /project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/atac_fragments.tsv.gz') %>%
  select(V1, V2, V3) %>% makeGRangesFromDataFrame(seqnames.field='V1', start.field='V2', end.field='V3', ignore.strand=T)

fragcounts <- CountFragments('/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/atac_fragments.tsv.gz')

atac.cells <- fragcounts %>% filter(frequency_count>2000) %>% select(CB) %>% pull()

atac.frags <- CreateFragmentObject(path='/project/gilad/kenneth/caQTL/highpass/batch1/lane1/outs/atac_fragments.tsv.gz',
                                   cells=atac.cells)

counts <- FeatureMatrix(fragments=atac.frags, features=Signac::granges(frags), cells=atac.cells)

