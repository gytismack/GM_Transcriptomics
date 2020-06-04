#!/usr/bin/env Rscript
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)

args = commandArgs(trailingOnly=TRUE)
FIRST = TRUE
for (arg in args){
    if (FIRST==TRUE){
        countdata = read.table(arg, header=TRUE)
        countdata = countdata[,c(1,7)]
        FIRST = FALSE
    }
    else{
        tmp = read.table(arg, header=TRUE)
        tmp = tmp[,c(1,7)]
        countdata = inner_join(countdata, tmp, by='Geneid')
    }
}
colnames(countdata) = gsub("\\.bam$", "", colnames(countdata))
colnames(countdata) = gsub("_NameSorted", "", colnames(countdata))
colnames(countdata) = gsub("out\\.BAMs\\.MAPPED_", "", colnames(countdata))
rownames(countdata) = countdata$Geneid
countdata = countdata[,2:ncol(countdata)]
cond = c("HBR", "UHRR")
types = c("Collibri", "KAPA")

for (t in types){
    tmp_matrix = select(countdata,contains(t))
    condition_v = c()
    for (name in colnames(tmp_matrix)){
        if (grepl(cond[1], name, fixed=T)){
            condition_v = append(condition_v, cond[1])
        }
        else if(grepl(cond[2], name, fixed=T)){
            condition_v = append(condition_v, cond[2])
        }
        else{
            print("cond NOT FOUND")
        }
    }
    condition = factor(condition_v)
    coldata = data.frame(row.names=colnames(tmp_matrix), condition)
    dds = DESeqDataSetFromMatrix(countData=tmp_matrix, colData=coldata, design=~condition)
    dds = DESeq(dds)
    res = results(dds)
    res = res[order(res$padj), ]
    
    plot = EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', xlim = c(-5, 8))
    ggsave(paste(t,"_vulcano.png", sep=""), plot=plot, device="png")

}

