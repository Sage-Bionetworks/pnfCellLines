# Quality Control Figures for RNA Seq data
# Xindi Guo

library(ggplot2)
library(dplyr)
library(synapseClient)
source("../../bin/RNASeqData.R")

minCount <- 2
minTpm <- 0.1

#### grch38 data
mat.tpm <- rnaKallistoMatrix()
mat.count <- rnaKallistoMatrix(metric = "est_counts")

annotes.grch38 <- samp.mappings[,c("Sample Name","RNA-Seq Data", "Sample Genotype")]
annotes.grch38 <- annotes.grch38[complete.cases(annotes.grch38),]
colnames(annotes.grch38) <- c("sampleName","synapseId","genotype")

# tpm
res <- tidyr::gather(as.data.frame(mat.tpm),key=synapseId,value=TPM)
res <- merge(res,annotes.grch38,by="synapseId")
res$Gene <- rep(rownames(mat.tpm),ncol(mat.tpm))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.tpm <- subset(res, TPM >= minTpm)
p<-ggplot(res.tpm,aes(y=TPM+1,x=genotype))+geom_boxplot(aes(fill=sampleName)) + scale_y_log10()
pdf('boxPlotOfTpmMin0.1.pdf')
print(p)
dev.off()

## count
res <- tidyr::gather(as.data.frame(mat.count),key=synapseId,value=count)
res <- merge(res,annotes.grch38,by="synapseId")
res$Gene <- rep(rownames(mat.count),ncol(mat.count))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.count <- subset(res, count >= minCount)
p<-ggplot(res.count,aes(y=count,x=genotype))+geom_boxplot(aes(fill=sampleName))+ scale_y_log10()
pdf('boxPlotOfCountMin2.pdf')
print(p)
dev.off()

#### gencode data
mat.tpm.gencode <- rnaGencodeKallistoMatrix()
mat.count.gencode <- rnaGencodeKallistoMatrix(metric = "est_counts")

annotes.gencode <- samp.mappings[,c("Sample Name","RNA-Seq Data (Gencode)", "Sample Genotype")]
annotes.gencode <- annotes.gencode[complete.cases(annotes.gencode),]
colnames(annotes.gencode) <- c("sampleName","synapseId","genotype")

# tpm
res <- tidyr::gather(as.data.frame(mat.tpm.gencode),key=synapseId,value=TPM)
res <- merge(res,annotes.gencode,by="synapseId")
res$Gene <- rep(rownames(mat.tpm.gencode),ncol(mat.tpm.gencode))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.tpm.gencode <- subset(res, TPM >= minTpm)
p<-ggplot(res.tpm.gencode,aes(y=TPM+1,x=genotype))+geom_boxplot(aes(fill=sampleName)) + scale_y_log10()
pdf('boxPlotOfTpmMin0.1_gencode.pdf')
print(p)
dev.off()

## count
res <- tidyr::gather(as.data.frame(mat.count.gencode),key=synapseId,value=count)
res <- merge(res,annotes.gencode,by="synapseId")
res$Gene <- rep(rownames(mat.count.gencode),ncol(mat.count.gencode))
res$Gene <- sapply(res$Gene,function(x) unlist(strsplit(x,split='.',fixed=T))[1])
res.count.gencode <- subset(res, count >= minCount)
p<-ggplot(res.count.gencode,aes(y=count,x=genotype))+geom_boxplot(aes(fill=sampleName))+ scale_y_log10()
pdf('boxPlotOfCountMin2_gencode.pdf')
print(p)
dev.off() 

#### upload filesw
rnaqc <- 'syn8360472'
scripturl <- 'https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2017-03-01/qc_figures_rnaSeq.R'

uploadFile2Synapse <- function(file,parentId,usedEnt,usedScript){
  synStore(File(path = file,parentId = parentId),
           used=list(list(entity=usedEnt,wasExecuted=FALSE),
                     list(url=usedScript,wasExecuted=TRUE)))
}
#the count files
uploadFile2Synapse(file = "boxPlotOfCountMin2.pdf",parentId = rnaqc,usedEnt = "syn5562376", usedScript = scripturl)
uploadFile2Synapse(file = "boxPlotOfCountMin2_gencode.pdf", parentId = rnaqc, usedEnt = "syn5580378", usedScript = scripturl)

#the tpm files
uploadFile2Synapse(file = "boxPlotOfTpmMin0.1.pdf", parentId = rnaqc, usedEnt = "syn5562378", usedScript = scripturl)
uploadFile2Synapse(file = "boxPlotOfTpmMin0.1_gencode.pdf", parentId = rnaqc, usedEnt = "syn5580347", usedScript = scripturl)
