#loading of installed in R environment
library(stringr)
library(DESeq2)
library(ggplot2)
library(ashr)
library(tidyverse)
library(dplyr)
library(apeglm)
library(vsn)
library(gplots)


#get working directory
getwd()


#list files in working directory
list.files()


#set working directory
setwd("C:/Users/rishi/Documents/BI342_Adv_Prog")


#read the csv file
counts=read.csv("MBI18001_counts.csv", sep=",", head = T, skip =1, row.names = "Geneid")


#extracting the required columns:diseased and controlled
countsnew <- counts [(-c(1:5))]
countsnew
colnames(countsnew)


#alternative code for extracting the required columns:diseased and controlled
colnames(countsnew)=str_split_fixed(colnames(countsnew),"\\.",6)[,1]
colnames(countsnew)


#Define conditions for the samples
mycols = data.frame(row.names = (colnames(countsnew)))
coldata <- data.frame(mycols, condition = factor(c(rep("combined", 3), rep("control",3))))
coldata

#check if row and column names from sample and count matrix matches.
all(rownames(coldata) %in% colnames(countsnew))
all(rownames(coldata) %in% colnames(countsnew))


#loading libraries for DESeq2
library(DESeq2)
library(gplots)
library(ggplot2)


# Build the data frame to be used by DESeq2 package
dds=DESeqDataSetFromMatrix(countData = countsnew,colData = coldata,design = ~ condition)
dds


# remove rows with zero
dds <- dds[rowSums(counts(dds)) > 10,]
dds <-DESeq(dds)


#define VST table
vst <- vst(dds, blind=FALSE)


#Plot PCA plot
plotPCA(vst, intgroup="condition", ntop=nrow(counts(dds)))


#Plot Box Plot
a <- DESeq2::plotPCA(vst, intgroup="condition")
a + geom_label(aes(label =coldata$condition),)
nudge <- position_nudge(y =1)
a + geom_label(aes(label = coldata$condition), position = nudge)
a + geom_text(aes(label = coldata$condition), position = nudge, size=3)
boxplot(assay(vst), col= c("Red", "Red", "Red", "Green", "Green", "Green"), pch=".",
        vertical=TRUE, cex.axis=0.5, main = "Boxplot of heat & drought affected rice using vst method",
        las=2, ylab="assay(vst)", xlab="Samples", ylin=c(-10,30),
        font.main=5, font.axis=0.5, font.lab=2)


#Plot correlation heatmap
cU <-cor( as.matrix(assay(vst)))
cols <- c("dodgerblue3", "firebrick3")[coldata$condition]
heatmap.2(cU, symm=TRUE, col= colorRampPalette(c("darkblue","white"))(100),
          labCol=colnames(cU), labRow=colnames(cU),
          distfun=function(c) as.dist(1 - c),
          trace="none",
          Colv=TRUE, cexRow=0.9, cexCol=0.9, key=F,
          font=2,
          RowSideColors=cols, ColSideColors=cols)
          

#Plot dispersion plot
plotDispEsts(dds)


# define filename and condition reprensentation for output file generation
res <- results(dds, contrast=c("condition", "combined","control"))
summary(res)
grp.mean <- sapply(levels(dds$condition),
                   function(lvl)
                     rowMeans(counts(dds,normalized=TRUE)[,dds$condition== lvl]))
norm.counts <- counts(dds, normalized=TRUE)
all <- data.frame(res, assay(vst))
nrow(all)
write.table(all, file="rice_combined_stress_main.csv",sep=",")


#write.table(assay(vst), file="tea_heat_stress1.csv",sep=",")



