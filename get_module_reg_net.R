# Run as ./get_module_reg_net.R tpm_matrix TF_list

# This script reads a TPM normalized GE matrix and a list of TF gene IDs (no header)
# It applies WGCNA to find modules of coexpressed genes and GENIE3 to find
# TF targets.

# Set the desired minimun module size in the minModuleSize parameter

# An output file called 'tpm_wgcna.modules.genesets' in written to the current directory
# This output file has genes in the first column and module IDs in the second.
# The second output file is 'tpm.TF-gene.genie3.top1perc.txt'
# Each line is an edge with third column depicting edge weight

library(WGCNA)
library(dplyr)
library(edgeR)
library(flashClust)
library(GENIE3)


rm(list=ls())


args = commandArgs(trailingOnly=TRUE) #Genes X Samples (TPM matrix)

tpm=read.table(args[1], header=T, sep="\t",stringsAsFactors=FALSE)
TF=read.table(args[2], header=F, sep="\t")

regulators.expmat=tpm[rownames(tpm) %in% TF$V1,]

#run wgcna to find modules
datExpr0=t(tpm)
net = blockwiseModules(datExpr0, power = 6, minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE, saveTOMFileBase = "tpm", verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)

genesets=as.data.frame(moduleLabels)

write.table(genesets, file="tpm_wgcna.modules.genesets", sep="\t",quote=F, col.names=FALSE)


#run genie
weightMatrix = GENIE3(as.matrix(expmat), regulators=rownames(regulators.expmat),nCores=32)
linkList = getLinkList(weightMatrix)

#Select top 1% edges
nedges=nrow(linkList)*0.01

write.table(linkList[1:nedges,],"tpm.TF-gene.genie3.top1perc.txt",col.names=TRUE, row.names=TRUE, sep="\t", quote=F)
