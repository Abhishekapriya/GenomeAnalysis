

#Installing DESeq required updation of R version and then the following command
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("apeglm")

install.packages ("plyr")
install.packages ("pheatmap")
install.packages ("RColorBrewer")
install.packages("ashr")

library("DESeq2")
library("plyr")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("ashr")

directory <- "E:/MS/Sem_2/Period_4/Genome_Analysis/Report/HtSeq"
sampleFiles<-c('RNA_BH_72_counts.txt', 'RNA_BH_73_counts.txt', 'RNA_BH_74_counts.txt', 'RNA_Serum_69_counts.txt', 'RNA_Serum_70_counts.txt', 'RNA_Serum_71_counts.txt')
sampleCondition<- c('BHI', 'BHI', 'BHI', 'Serum', 'Serum', 'Serum')
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition)

# Differential Expression analysis
dds <- DESeq(dds)

# Filtering out reads with less than 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Setting condition level to the control
dds$condition <- relevel(dds$condition, ref = "BHI")
res <- results(dds)

# Obtaining genes with adjusted P-value cutoff (alpha) set to 0.05: 
res1 <- results(dds, alpha=0.05)

summary(res)
summary(res1)
res <- results(dds, name=resultsNames(dds)[2])
res

# Transforming the result data to show Log Fold Change for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_Serum_vs_BHI", type="apeglm")

# Ordering results by their p-value
resOrdered <- res[order(res$pvalue),]

# Filtering and identifying number of genes with P-adjusted value less than 0.001: 
filterL2FCG2L05 <- (abs(resOrdered$log2FoldChange) > 2) | (abs(resOrdered$log2FoldChange) < 0.5)
filtered_res <- resOrdered[filterL2FCG2L05,]
sum(filtered_res$padj < 0.001, na.rm=TRUE)


# Identifying and Ordering the top 100 genes genes with most significant expression by adjusted p-value for each sequencing run of RNA from human serum
n = 100
topResults <- rbind( resOrdered[ abs(resOrdered[,'log2FoldChange']) > 2,][1:n,])
topResults[c(60:100), c('baseMean','log2FoldChange','padj')]


# Obtaining MA plot for normalized DEseq results
plotMA(results(dds),main = "MA Plot of dds",ylim=c(-4,4), colNonSig = "gray", colSig = "purple")
plotMA(resLFC, main = "MA Plot of LFC of result data", ylim=c(-2,2), colNonSig = "gray", colSig = "purple")
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

# Obtaining MA plots for different shrinkage factors
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
rld <- rlog(dds, blind=TRUE)
head(assay(vsd), 3)

# Obtaining log2(n + 1)
ntd <- normTransform(dds)

# Obtaning heat maps for different transformation factors
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","sizeFactor")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

# Obatning PCA plot of condition vs counts
data <- plotPCA(rld, intgroup=c("condition", "sizeFactor"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

# Obtaining heat map of top results
hmcol <- brewer.pal(11,'RdBu')
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))

# Loading merged counts data for a full heatmap display
count.table <- read.table(file="merged.tmp", sep="\t", header=FALSE, row.names = 1)
epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(count.table + epsilon)), breaks=100, col="lightblue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)