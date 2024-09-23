# DESeq2 analysis
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

# package library
library(DESeq2)
library(ggplot2)
library(readr)
library(dplyr)
library(pheatmap)

setwd("/Users/yongbin/projects/proj_deg_scz/")
figurepath = "./figures/"

# Load count data
df = read_csv("./data/data_count.csv")

# ================= QC ================= 
# Check the distribution of count data
dftmp = filter(df, SZ40_count < 10000)
gcf = ggplot(dftmp) +
  geom_histogram(aes(x = SZ40_count), stat = "bin", bins = 100) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
gcf
ggsave(paste(figurepath, "Histogram_count.png"), dpi = 300, device = 'png', gcf)

# Check mean vs. variance
mean_counts <- apply(df[,2:8], 1, mean)        
variance_counts <- apply(df[,2:8], 1, var)
df1 <- data.frame(mean_counts, variance_counts)
gcf1 = ggplot(df1) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")
gcf1
ggsave(paste(figurepath, "Mean_variance.png"), dpi = 300, device = 'png', gcf1)

# construct coldata
coldata <- read.csv("./data/data_meta.csv", row.names=1)
coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)
coldata$age <- factor(coldata$age)
coldata$edu <- factor(coldata$edu)

# data normalization
cts <- as.matrix(read.csv("./data/data_count.csv", row.names=1))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="./data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# PCA
rld <- vst(dds, blind=TRUE)
plotPCA(rld, intgroup="condition")

# heatmap
rld_mat <- assay(rld)    
rld_cor <- cor(rld_mat)    ## cor() is a base R function
head(rld_cor)
pheatmap(rld_cor, annotation = coldata)


# ================= RUN DESeq2 ====================
# DEG test
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition + age + sex)
dds <- DESeq(dds)
plotDispEsts(dds)

# show results
res <- results(dds, contrast=c("condition", "SZ", "HS"), independentFiltering = F)
summary(res)

# lfcShrink
res_shrink <- lfcShrink(dds, coef = "condition_SZ_vs_HS", type = "apeglm")
resValid <- res_shrink[which(!is.na(res$padj)),]
plotMA(res, ylim=c(-2,2))
plotMA(res_shrink, ylim=c(-2,2))


# ============= Visualization (volcano plot) ===============
res_table <- as.data.frame(resValid) 
resdata <- merge(as.data.frame(resValid), 
                 as.data.frame(counts(dds, normalized=TRUE)), 
                 by="row.names", sort=FALSE)
nonkeep <- rowSums(resdata[,7:109]==0) >= 10
res_table[nonkeep, 'padj']=NA

res_table2 <- res_table %>% 
  mutate(threshold_OE = padj < 0.05)
ggplot(res_table2) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Differential expression in SCZ") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  


# ============= Merge results ============= 
resdata <- merge(res_table, 
                 as.data.frame(counts(dds, normalized=TRUE)), 
                 by="row.names", sort=FALSE)
resdata$significant <- "unchanged"
resdata$significant[resdata$padj <= 0.05 & resdata$log2FoldChange >= 0 ] <- "upregulated"
resdata$significant[resdata$padj <= 0.05 & resdata$log2FoldChange < 0 ] <- "downregulated"
write.csv(as.data.frame(resdata),
          file="./results/DESeq2_deg_results_all.csv")
