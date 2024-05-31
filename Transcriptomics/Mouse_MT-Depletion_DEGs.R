# Clear all variables in the environment
rm(list=ls(all=T))

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationHub)
library(org.Mm.eg.db)
library(clusterProfiler)
library(plyr)

# Set working directory
setwd("~/Mouse_Mito-depletion/Rresult")

# Define sample numbers
controlnum <- 2
casenum <- 2
num <- casenum + controlnum

# Load data
datad4 <- as.matrix(read.csv("../RNAdata/genematrix/gene_count_matrix.csv", row.names="gene_id"))
colnames(datad4) <- c("WT_iPS1", "WT_iPS2", "MD_D4_iPS1", "MD_D4_iPS2")
datad4 <- datad4[, c(3, 4, 1, 2)]

# Compare Day4 with WT
sampletitle <- "Mouse MT_D4 vs WT_iPS"
countData <- datad4
condition <- factor(c(rep("treat", casenum), rep("untreat", controlnum)))
colData <- data.frame(row.names=colnames(countData), condition)

# Differential expression analysis
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize=TRUE)), by="row.names", sort=FALSE)
resdata <- resdata[!is.na(resdata$padj),]
resdata <- resdata[order(resdata$padj),]

# Gene annotation
PluripotencyGene <- read.table("Pluripotency_Genes.txt", header = F)
MT_coding <- read.table("Human mitochondria genome coding genes.csv", sep = ",", header = T)
MT_asso <- read.table("Mouse_chromosome_related_mito_genelist.txt", sep = "\t", header = F)
resdata$genename <- sapply(resdata$Row.names, function(x) strsplit(x, "\\|")[[1]][2])
resdata$genename <- toupper(resdata$genename)
resdata$genetype <- NA
resdata$genetype[resdata$genename %in% toupper(c(MT_coding[, 2], "MT-CYTB"))] <- "MT_coding"
resdata$genetype[resdata$genename %in% toupper(MT_asso[, 1])] <- "MT_associated"
resdata$genetype[resdata$genename %in% toupper(PluripotencyGene[, 1])] <- "PluripotencyGenes"

# PCA plot for RNA-seq
library(ggbiplot)
pcamat <- as.data.frame(t(resdata[, 8:11]))
rownames(pcamat) <- c("MT_Depletion_Day4_S1", "MT_Depletion_Day4_S2", "WT_Day0_S1", "WT_Day0_S2")
pcamatgroup <- c(rep("MT_Depletion", 2), rep("WT", 2))
pca_result <- prcomp(pcamat, scale = TRUE, center = FALSE)
pca_result$x <- scale(pca_result$x)

pdf(file = "Mouse_MT_Depletion_pcaplot.pdf", width = 5, height = 5)
ggbiplot(pca_result,
         var.axes = FALSE,
         obs.scale = 1,
         groups = pcamatgroup,
         ellipse = FALSE,
         circle = FALSE) +
  geom_text(aes(label = rownames(pcamat)), vjust = 1.5, size = 2) +
  scale_x_continuous(limits = c(-1.1, 1.1)) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  theme_bw()
dev.off()

# Volcano plot and heatmap for differentially expressed genes
diff_gene_p <- resdata[resdata$padj < 0.05, ]
diff_gene_p$stat[diff_gene_p$log2FoldChange < 0] <- "Up"
diff_gene_p$stat[diff_gene_p$log2FoldChange > 0] <- "Down"
diff_gene_p$log2FoldChange <- -diff_gene_p$log2FoldChange
diff_gene_p <- diff_gene_p[abs(diff_gene_p$log2FoldChange) > 1.5, ]
write.table(diff_gene_p, sep = "\t", row.names = F, quote = F, "Mouse_MD_DEGs.txt")

# Volcano plot
resdata$log2FoldChange <- -resdata$log2FoldChange
data <- resdata
data$threshold <- as.factor(ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1.5, ifelse(data$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
data$pos <- rowMeans(data[, 8:(8 + casenum - 1)])
data$neg <- rowMeans(data[, (8 + casenum):(8 + num - 1)])
data$genetype[!data$Row.names %in% diff_gene_p$Row.names] <- NA
data$genename[is.na(data$genetype)] <- NA
p <- ggplot(data, aes(x = log2(pos), y = log2(neg), colour = threshold)) + geom_point()
p1 <- p + geom_point(aes(color = data$threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha = 0.4, size = 1)
p3 <- p2 + theme_bw() + labs(x = "Case_genes", y = "Control_genes", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 0.6)
# Calculate and annotate Pearson correlation
a <- cor.test(data$pos, data$neg, method = "pearson")
p6 <- p5 + annotate("text", x = 14, y = -0.5, label = paste("Pearson's correlation ", round(a$estimate, 3)), size = 3)
ggsave(p6, file = paste(sampletitle, "_correlation2.pdf", sep = ""), width = 8, height = 8)

# Prepare data for volcano plot
vdata <- arrange(data, padj)
vdata$threshold <- as.factor(ifelse(vdata$pvalue < 0.05 & abs(vdata$log2FoldChange) >= 1.5, ifelse(vdata$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
vdata$genetype[vdata$threshold == "Not"] <- "Not"
vdata$genetype <- factor(vdata$genetype, levels = c("MT_associated", "Not", "MT_coding", "PluripotencyGenes"))
# General volcano plot
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% c("MT_coding", "MT_associated", "PluripotencyGenes")),
             aes(shape = threshold, color = genetype), size = 1) +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#f8766d", "lightgray", "green", "#0033A0"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-max(abs(vdata$log2FoldChange)), max(abs(vdata$log2FoldChange))))
ggsave(v4, file = paste(sampletitle, "_volcano3.pdf", sep = ""), width = 5, height = 4)
# MT_associated Volcano Plot
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(data = subset(vdata, !genetype %in% "MT_associated"), aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% "MT_associated"), aes(shape = threshold, color = genetype), size = 2, color = "#f8766d") +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("lightgray", "#7F7F7F", "#7F7F7F"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(-4, 160) + xlim(-11, 11)
ggsave(v4, file = "Mouse_RNA_volcano_MT_associated.pdf", width = 5, height = 4)
# Pluripotency Genes Volcano Plot
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(data = subset(vdata, !genetype %in% "PluripotencyGenes"), aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% "PluripotencyGenes"), aes(shape = threshold, color = genetype), size = 2, color = "#0033A0") +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#7F7F7F", "lightgray", "#7F7F7F"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(-4, 160) + xlim(-11, 11)
ggsave(v4, file = "Mouse_RNA_volcano_P.pdf", width = 5, height = 4)
# MT_coding Genes Volcano Plot
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(data = subset(vdata, !genetype %in% "MT_coding"), aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% "MT_coding"), aes(shape = threshold, color = genetype), size = 2, color = "green") +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#7F7F7F", "lightgray", "#7F7F7F"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(-4, 160) + xlim(-11, 11)
ggsave(v4, file = "Mouse_RNA_volcano_MT_coding.pdf", width = 5, height = 4)

## Heatmap
# Prepare the data for heatmap
datamat <- as.matrix(diff_gene_p[, c(8:11)])
rownames(datamat) <- diff_gene_p$genename
anno_gene2 <- diff_gene_p[!is.na(diff_gene_p$genetype), c(12, 13)]
rownames(datamat)[rownames(datamat) == "NA"] <- ""
# Create a complex heatmap using ComplexHeatmap package
library(ComplexHeatmap)
# Convert the expression matrix to log scale and normalize
A <- as.matrix(datamat)
samples <- rep(c('MD_D4', 'WT_iPS'), c(2, 2))
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))
# Define the color palette and annotations
B <- Heatmap(A, col = colorRampPalette(c("blue", "white", "red"))(100),
             show_row_names = FALSE,
             top_annotation = HeatmapAnnotation(Group = samples,
                                                simple_anno_size = unit(2, 'mm'),
                                                col = list(Group = c('MD_D4' = '#00DAE0', 'WT_iPS' = '#FF9289')),
                                                show_annotation_name = FALSE))
# Add row annotations for specific gene types
C <- B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "MT_associated"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "MT_associated"], labels_gp = gpar(fontsize = 2.5, col = "#f8766d")))
D <- C + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "PluripotencyGenes"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "PluripotencyGenes"], labels_gp = gpar(fontsize = 4, col = "blue")))
E <- D + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "MT_coding"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "MT_coding"], labels_gp = gpar(fontsize = 4, col = "green")))
# Save the complex heatmap to a PDF file
pdf(file = "Mouse_DEGs6.pdf", width = 7, height = 15)
print(E)
dev.off()

# Calculate and plot the percentage of genes in each set
resdata$threshold <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1.5, ifelse(resdata$log2FoldChange > 0, 'Up', 'Down'), 'Not'))
merge_bar_all <- cbind(as.matrix(resdata$threshold), as.matrix(resdata$genetype))
# Initialize an empty data frame for storing table information
table_type <- data.frame()
# Loop through gene sets and calculate statistics
for (i in c("PluripotencyGenes", "MT_associated", "MT_coding")) {
  sign <- table(merge_bar_all[merge_bar_all[, 2] == i, 1])
  singletype_table <- cbind(rep(i, length(sign)), names(sign), as.matrix(sign))
  rownames(singletype_table) <- NULL
  table_type <- rbind(table_type, singletype_table)
}
# Define column names and convert data types
colnames(table_type) <- c("Celltype", "DEGtype", "Genenumber")
table_sample_type <- as.data.frame(table_type)
table_sample_type$DEGtype <- as.factor(table_sample_type$DEGtype)
table_sample_type$Genenumber <- as.numeric(as.matrix(table_sample_type$Genenumber))
# Calculate percentage weights and prepare for plotting
ce <- ddply(table_sample_type, "Celltype", transform, percent_weight = Genenumber / sum(Genenumber) * 100)
ce$DEGtype <- factor(ce$DEGtype, levels = c("Up", "Not", "Down"))
ce$Celltype <- factor(ce$Celltype, levels = c("MT_associated", "PluripotencyGenes", "MT_coding"))
ce <- ddply(ce, "Celltype", transform, labely = cumsum(percent_weight))
# Create the bar plot
p <- ggplot(ce, aes(x = Celltype, y = percent_weight, fill = DEGtype)) + geom_bar(stat = "identity") + theme_bw()
p_bar2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#FF9289", "lightgray", "#00DAE0")) +
  scale_y_continuous(name = "Gene Percent", expand = c(0, 0)) + labs(x = "Gene Sets", colour = "") +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(y = labely, label = Genenumber), vjust = 1.5, colour = "black", size = 6)
# Save the bar plot to a PDF file
ggsave("Mouse_DEGs_Geneset_percent313.pdf", p_bar2, width = 6, height = 6)

## Gene Function Annotation for Mouse
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
# Define the function for GO enrichment analysis and plotting
perform_GO_enrichment <- function(gene_list, filename_prefix, title_prefix, color_fill) {
  geneLists <- data.frame(symbol = gene_list)
  results <- merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
  id_list <- na.omit(results$gene_id)
  
  ego <- enrichGO(OrgDb = "org.Mm.eg.db", gene = id_list, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
  ego_result <- ego@result[ego@result$Count >= 10 & ego@result$pvalue < 0.05, ]
  write.csv(ego_result, paste0(filename_prefix, "_GO.csv"))
  
  ego_draw <- ego_result[1:10, ]
  p <- ggplot(ego_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
    geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
    xlab("") + scale_y_continuous(name = "Gene Count") + coord_flip() +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient(low = color_fill, high = color_fill, name = "Genes") +
    theme(text = element_text(size = 16, family = "serif")) + labs(title = paste0(title_prefix, " GO Annotation"))
  
  ggsave(p, filename = paste0(filename_prefix, "_GO.pdf"), height = 5, width = 15)
}
# Define the function for KEGG enrichment analysis and plotting
perform_KEGG_enrichment <- function(gene_list, filename_prefix, title_prefix, color_fill) {
  geneLists <- data.frame(symbol = gene_list)
  results <- merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
  id_list <- na.omit(results$gene_id)
  
  kegg <- enrichKEGG(organism = "mmu", gene = id_list, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
  kegg <- setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kegg_result <- kegg@result[kegg@result$Count >= 10 & kegg@result$pvalue < 0.05, ]
  write.csv(kegg_result, paste0(filename_prefix, "_KEGG.csv"))
  
  kegg_draw <- kegg_result[1:10, ]
  p <- ggplot(kegg_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
    geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
    xlab("") + scale_y_continuous(name = "Gene Count") + coord_flip() +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient(low = color_fill, high = color_fill, name = "Genes") +
    theme(text = element_text(size = 16, family = "serif")) + labs(title = paste0(title_prefix, " KEGG Annotation"))
  
  ggsave(p, filename = paste0(filename_prefix, "_KEGG.pdf"), height = 5, width = 10)
}
# Get the list of upregulated genes
diff_gene_p_up <- diff_gene_p$genename[diff_gene_p$stat == "Up"]
perform_GO_enrichment(diff_gene_p_up, "Mouse_mito_depletion_deg_UP", "Mouse Mitochondria Depletion Up", "#FF9289")
perform_KEGG_enrichment(diff_gene_p_up, "Mouse_mito_depletion_deg_UP", "Mouse Mitochondria Depletion Up", "#FF9289")
# Get the list of downregulated genes
diff_gene_p_down <- na.omit(diff_gene_p$genename[diff_gene_p$stat == "Down"])
perform_GO_enrichment(diff_gene_p_down, "Mouse_mito_depletion_deg_down", "Mouse Mitochondria Depletion Down", "#00DAE0")
perform_KEGG_enrichment(diff_gene_p_down, "Mouse_mito_depletion_deg_down", "Mouse Mitochondria Depletion Down", "#00DAE0")
# Save the workspace image
save.image(file = "Mouse_Mito_Depletion.RDATA")
