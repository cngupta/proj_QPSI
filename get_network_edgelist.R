library(data.table)
library(WGCNA)
library(readxl)
library(minet)
library(igraph)
library(dplyr)
library(ggplot2)

rm(list=ls())

source('../scripts/exportNetToCytoscape.R')
#functions
filterGexpr=function(mat,tpmTh,sampleTh)
{
  th=round(ncol(mat)*sampleTh)
  print(paste0("Threshold:",th))
  keep <- rowSums(mat > tpmTh) >= th
  mat.filt=mat[keep,]
  mat.filt
}

get_modules=function(mat)
{
  datExpr0=t(mat)
  net = blockwiseModules(datExpr0, power = 6, minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 3)
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  genesets=as.data.frame(moduleLabels)
  print(dim(genesets))
  print(length(unique(rownames(genesets))))
  genesets
}

#get edge list based on mutual information usign package minet
get_edgeList=function(mat,th=0.2)#report top 20% edges
{
  datExpr0=t(mat)
  net = build.mim(datExpr0, estimator = "pearson", nbins = nrow(datExpr0)/3)
  clr=clr(net, skipDiagonal=1) #play with one witohout CLR
  print(dim(clr))
  g=graph_from_adjacency_matrix( clr,mode="undirected",weighted=TRUE,diag=FALSE)
  g
  #df=adjTodf(clr,edgeFile=NULL,nodeFile = NULL,threshold=5)
  #edges=df$edgeData
  #edges=edges[order(-edges$weight),]
  #edges=edges[1:(nrow(edges)*th),]
  #edges
}


#metadata
meta=read_xlsx("/Users/chiraggupta/Desktop/QPSI/phase1/20400Tim_metadata_comparison matrix_V2_MB.xlsx", col_names=T, sheet=2)[1:9]
colnames(meta)=c("No.","2","Position","Sample_ID","Tissue","Treatment","Timepoint","Species","Extended")
meta=as.data.frame(meta)
meta$Extended=gsub("-",".",meta$Extended)

meta.pop.leaf.fe=meta[(meta$Species %like% "Populus" & meta$Treatment %like% "Fe" & meta$Tissue %like% "Leaf"),]
meta.pop.root.fe=meta[(meta$Species %like% "Populus" & meta$Treatment %like% "Fe" & meta$Tissue %like% "Root"),]
meta.pop.leaf.zn=meta[(meta$Species %like% "Populus" & meta$Treatment %like% "Zn" & meta$Tissue %like% "Leaf"),]
meta.pop.root.zn=meta[(meta$Species %like% "Populus" & meta$Treatment %like% "Zn" & meta$Tissue %like% "Root"),]

meta.sor.leaf.fe=meta[(meta$Species %like% "Sorghum" & meta$Treatment %like% "Fe" & meta$Tissue %like% "Leaf"),]
meta.sor.root.fe=meta[(meta$Species %like% "Sorghum" & meta$Treatment %like% "Fe" & meta$Tissue %like% "Root"),]
meta.sor.leaf.zn=meta[(meta$Species %like% "Sorghum" & meta$Treatment %like% "Zn" & meta$Tissue %like% "Leaf"),]
meta.sor.root.zn=meta[(meta$Species %like% "Sorghum" & meta$Treatment %like% "Zn" & meta$Tissue %like% "Root"),]


pop=read.table("/Users/chiraggupta/Desktop/QPSI/data/poplar/20400Tim-raw_genes_tpm.csv", sep=",", header=T, row.names=1)
sor=read.table("/Users/chiraggupta/Desktop/QPSI/data/sorghum/20400Tim-raw_genes_tpm.csv", sep=",", header=T, row.names=1)

colnames(pop)=gsub("X","",colnames(pop))
colnames(sor)=gsub("X","",colnames(sor))


pop.leaf.fe.expr=pop[,colnames(pop) %in% meta.pop.leaf.fe$Extended]
pop.root.fe.expr=pop[,colnames(pop) %in% meta.pop.root.fe$Extended]
pop.leaf.zn.expr=pop[,colnames(pop) %in% meta.pop.leaf.zn$Extended]
pop.root.zn.expr=pop[,colnames(pop) %in% meta.pop.root.zn$Extended]


sor.leaf.fe.expr=sor[,colnames(sor) %in% meta.sor.leaf.fe$Extended]
#remove one bad sample
sor.leaf.fe.expr=sor.leaf.fe.expr[,!(colnames(sor.leaf.fe.expr) %like% "20400Tim_S56L.4")]

sor.root.fe.expr=sor[,colnames(sor) %in% meta.sor.root.fe$Extended]
#remove one bad sample
sor.root.fe.expr=sor.root.fe.expr[,!(colnames(sor.root.fe.expr) %like% "20400Tim_S73R.4")]

sor.leaf.zn.expr=sor[,colnames(sor) %in% meta.sor.leaf.zn$Extended]
sor.root.zn.expr=sor[,colnames(sor) %in% meta.sor.root.zn$Extended]


#filtering:
#keep genes that have TPM values above 1 in at least 20% libraries:

pop.leaf.fe.expr=filterGexpr(pop.leaf.fe.expr,1,0.2)
pop.root.fe.expr=filterGexpr(pop.root.fe.expr,1,0.2)
pop.leaf.zn.expr=filterGexpr(pop.leaf.zn.expr,1,0.2)
pop.root.zn.expr=filterGexpr(pop.root.zn.expr,1,0.2)

sor.leaf.fe.expr=filterGexpr(sor.leaf.fe.expr,1,0.2)
sor.root.fe.expr=filterGexpr(sor.root.fe.expr,1,0.2)
sor.leaf.zn.expr=filterGexpr(sor.leaf.zn.expr,1,0.2)
sor.root.zn.expr=filterGexpr(sor.root.zn.expr,1,0.2)


#Network inference
#Poplar
pop.leaf.fe.edgeList=get_edgeList(pop.leaf.fe.expr) #try louvain clustering on these networks using igraph in R
pop.root.fe.edgeList=get_edgeList(pop.root.fe.expr)
pop.leaf.zn.edgeList=get_edgeList(pop.leaf.zn.expr)
pop.root.zn.edgeList=get_edgeList(pop.root.zn.expr)


#sorghum
sor.leaf.fe.edgeList=get_edgeList(sor.leaf.fe.expr)
sor.root.fe.edgeList=get_edgeList(sor.root.fe.expr)
sor.leaf.zn.edgeList=get_edgeList(sor.leaf.zn.expr)
sor.root.zn.edgeList=get_edgeList(sor.root.zn.expr)


g=pop.leaf.fe.edgeList
pop.leaf.fe.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
pop.leaf.fe.edgeList.strength=as.data.frame(pop.leaf.fe.edgeList.strength)
pop.leaf.fe.edgeList.strength$gene=rownames(pop.leaf.fe.edgeList.strength)

g=pop.root.fe.edgeList
pop.root.fe.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
pop.root.fe.edgeList.strength=as.data.frame(pop.root.fe.edgeList.strength)
pop.root.fe.edgeList.strength$gene=rownames(pop.root.fe.edgeList.strength)


g=pop.leaf.zn.edgeList
pop.leaf.zn.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
pop.leaf.zn.edgeList.strength=as.data.frame(pop.leaf.zn.edgeList.strength)
pop.leaf.zn.edgeList.strength$gene=rownames(pop.leaf.zn.edgeList.strength)


g=pop.root.zn.edgeList
pop.root.zn.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
pop.root.zn.edgeList.strength=as.data.frame(pop.root.zn.edgeList.strength)
pop.root.zn.edgeList.strength$gene=rownames(pop.root.zn.edgeList.strength)


df.pop.root=left_join(pop.root.zn.edgeList.strength,pop.root.fe.edgeList.strength)
df=df.pop.root
p=ggplot(df,aes(x=df[,1],y=df[,3]))+geom_point()+labs(x = "degree zn", y = "degree fe")+
ggtitle("Poplar root")
ggsave(p,file="Figures/df.pop.root.scatter.pdf")


df.pop.leaf=left_join(pop.leaf.zn.edgeList.strength,pop.leaf.fe.edgeList.strength)
df=df.pop.leaf
p=ggplot(df,aes(x=df[,1],y=df[,3]))+geom_point()+labs(x = "degree zn", y = "degree fe")+
ggtitle("Poplar leaf")
ggsave(p,file="Figures/df.pop.leaf.scatter.pdf")


df.pop.zn=left_join(pop.leaf.zn.edgeList.strength,pop.root.zn.edgeList.strength)
df=df.pop.zn
p=ggplot(df,aes(x=df[,1],y=df[,3]))+geom_point()+labs(x = "degree leaf", y = "degree root")+
ggtitle("Poplar Zinc")
ggsave(p,file="Figures/df.pop.Zn.scatter.pdf")


df.pop.fe=left_join(pop.leaf.fe.edgeList.strength,pop.root.fe.edgeList.strength)
df=df.pop.fe
p=ggplot(df,aes(x=df[,1],y=df[,3]))+geom_point()+labs(x = "degree leaf", y = "degree root")+
ggtitle("Poplar Iron")
ggsave(p,file="Figures/df.pop.Fe.scatter.pdf")



#############
g=sor.leaf.fe.edgeList
sor.leaf.fe.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
sor.leaf.fe.edgeList.strength=as.data.frame(sor.leaf.fe.edgeList.strength)
sor.leaf.fe.edgeList.strength$gene=rownames(sor.leaf.fe.edgeList.strength)

g=sor.root.fe.edgeList
sor.root.fe.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
sor.root.fe.edgeList.strength=as.data.frame(sor.root.fe.edgeList.strength)
sor.root.fe.edgeList.strength$gene=rownames(sor.root.fe.edgeList.strength)


g=sor.leaf.zn.edgeList
sor.leaf.zn.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
sor.leaf.zn.edgeList.strength=as.data.frame(sor.leaf.zn.edgeList.strength)
sor.leaf.zn.edgeList.strength$gene=rownames(sor.leaf.zn.edgeList.strength)


g=sor.root.zn.edgeList
sor.root.zn.edgeList.strength=strength(g,vids = V(g),mode="all",loops = FALSE)
sor.root.zn.edgeList.strength=as.data.frame(sor.root.zn.edgeList.strength)
sor.root.zn.edgeList.strength$gene=rownames(sor.root.zn.edgeList.strength)


df.sor.root=left_join(sor.root.zn.edgeList.strength,sor.root.fe.edgeList.strength)
df=df.sor.root
p=ggplot(df,aes(x=df[,1],y=df[,3]))+geom_point()
ggsave(p,file="Figures/df.sor.root.scatter.pdf")


df.sor.leaf=left_join(sor.leaf.zn.edgeList.strength,sor.leaf.fe.edgeList.strength)
df=df.sor.leaf
p=ggplot(df,aes(x=df[,1],y=df[,3]))+geom_point()
ggsave(p,file="Figures/df.sor.leaf.scatter.pdf")


#write files
#write.table(count_df, file=paste("clr_outs","count_table.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)

#write.table(pop.leaf.fe.edgeList[,1:3], file=paste("clr_outs","pop.leaf.fe.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)
#write.table(pop.root.fe.edgeList[,1:3], file=paste("clr_outs","pop.root.fe.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)
#write.table(pop.leaf.zn.edgeList[,1:3], file=paste("clr_outs","pop.leaf.zn.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)
#write.table(pop.root.zn.edgeList[,1:3], file=paste("clr_outs","pop.root.zn.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)

#write.table(sor.leaf.fe.edgeList[,1:3], file=paste("clr_outs","sor.leaf.fe.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)
#write.table(sor.root.fe.edgeList[,1:3], file=paste("clr_outs","sor.root.fe.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)
#write.table(sor.leaf.zn.edgeList[,1:3], file=paste("clr_outs","sor.leaf.zn.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)
#write.table(sor.root.zn.edgeList[,1:3], file=paste("clr_outs","sor.root.zn.edgeList.clr.txt",sep="/"), sep="\t",quote=F, col.names=TRUE)







#just an example
tmp=pop.root.zn.expr[c("Potri.001G000400","Potri.001G001600"),]
tmp=t(tmp)
pdf(p,file="Figures/exampleMI.pdf")
scatter.smooth(tmp[,1],tmp[,2])
dev.off()
cor.test(tmp[,1],tmp[,2])
