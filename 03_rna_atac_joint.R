library(Signac)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep = "")
setwd('/project2/gilad/daraujo/scRNA_scATAC')

# load datasets
rna.assay <- LoadSeuratRds('compiled_rna_processed.rds')
atac.assay <- LoadSeuratRds('compiled_atac_processed.rds') 

# add dataset-identifying metadata
rna.assay$dataset <- 'RNA'
atac.assay$dataset <- 'ATAC'

# umaps
p1 <- DimPlot(rna.assay) + ggtitle('RNA')
p2 <- DimPlot(atac.assay) + ggtitle('ATAC')
p1 + p2
ggsave('umap_scRNA.scATAC_processed_bybatch.lane.pdf', height=5, width=14)

# quantify gene activity
gene.activities <- GeneActivity(atac.assay, features=VariableFeatures(rna.assay))
atac.assay[['ACTIVITY']] <- CreateAssayObject(counts=gene.activities)

# identify anchors
transfer.anchors <- FindTransferAnchors(reference=rna.assay, query=atac.assay, 
                                        features=VariableFeatures(object=rna.assay),
                                        reference.assay='RNA', 
                                        query.assay='ACTIVITY', 
                                        reduction='cca')

# co-embed scRNA-seq and scATAC-seq datasets 
var.genes <- VariableFeatures(rna.assay)
refdata <- GetAssayData(JoinLayers(rna.assay), assay='RNA', slot='data')[var.genes,]
####need to load Seurat4.4.0 and run next code
atac.assay[['RNA']] <- TransferData(anchorset=transfer.anchors, refdata=refdata, 
                                    weight.reduction=atac.assay[['lsi']],
                                    dims=2:30) 
####load Seurat5.1.0 again
coembed <- merge(x=rna.assay, y=atac.assay)
#coembed[['RNA']] <- as(object=coembed[['RNA']], Class='Assay')
coembed <- ScaleData(coembed, features=var.genes, do.scale=F)
coembed <- RunPCA(coembed, features=var.genes, verbose=F)
coembed <- RunUMAP(coembed, dims=1:30)
DimPlot(coembed, group.by=c('dataset'))
ggsave('umap_scRNA.scATAC.joint_by.dataset.pdf', height=5, width=7)
SaveSeuratRds(coembed, file='rna_atac_joint.rds')