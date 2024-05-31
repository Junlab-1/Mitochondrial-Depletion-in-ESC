# Clear all objects from the workspace
rm(list = ls(all = TRUE))

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(plyr)

# Set working directory
setwd("~/Human_Mito-depletion/12.5update")

# Define the number of control and case samples
controlnum = 2
casenum = 2
num = casenum + controlnum

# Load gene count matrix data
dataall <- as.matrix(read.csv("../RNAdata/genematrix/gene_count_matrix.csv", row.names = "gene_id"))
datad4 <- dataall[, c(5, 6, 3, 4)]
colnames(datad4) <- c("MD_D4_iPS1", "MD_D4_iPS2", "WT_iPS1", "WT_iPS2")

# Compare Day4 with WT
sampletitle = "Human MT_D4 vs WT_iPS"
countData <- datad4
condition <- factor(c(rep("treat", casenum), rep("untreat", controlnum)))
colData <- data.frame(row.names = colnames(countData), condition)

# Create DESeq2 dataset and run differential expression analysis
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Merge results with normalized counts
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize = TRUE)), by = "row.names", sort = FALSE)
resdata <- resdata[!is.na(resdata$padj), ]
resdata <- resdata[order(resdata$padj), ]

# Gene annotation
PluripotencyGene <- read.table("Pluripotency_Genes.txt", header = FALSE)
MT_coding <- read.table("Human mitochondria genome coding genes.csv", sep = ",", header = TRUE)
MT_asso <- read.table("Human_chromosome_mito_genelist.csv", sep = ",", header = TRUE)

# Annotate gene names and types
resdata$genename <- sapply(resdata$Row.names, function(x) strsplit(x, "\\|")[[1]][2])
resdata$genetype <- NA
resdata$genetype[resdata$genename %in% MT_coding[, 2]] <- "MT_coding"
resdata$genetype[resdata$genename %in% MT_asso[, 2]] <- "MT_associated"
resdata$genetype[resdata$genename %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"

# PCA plot for RNA-seq data using Day4
library(ggbiplot)
pcamat <- as.data.frame(t(resdata[, 8:11]))
rownames(pcamat) <- c("MT_Depletion_Day4_S1", "MT_Depletion_Day4_S2", "WT_Day0_S1", "WT_Day0_S2")
pcamatgroup <- c(rep("MT_Depletion4", 2), rep("WT", 2))
pca_result <- prcomp(pcamat, scale = TRUE, center = FALSE)
pca_result$x <- scale(pca_result$x)
pdf(file = "Human_MT_Depletion_pcaplot.pdf", width = 6, height = 6)
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

diff_gene_p <- resdata[resdata$padj < 0.05, ]
diff_gene_p$stat[diff_gene_p$log2FoldChange < 0] <- "Up"
diff_gene_p$stat[diff_gene_p$log2FoldChange > 0] <- "Down"
diff_gene_p$log2FoldChange <- -diff_gene_p$log2FoldChange
diff_gene_p <- diff_gene_p[abs(diff_gene_p$log2FoldChange) > 1.5, ]
resdata$log2FoldChange=-resdata$log2FoldChange
write.table(diff_gene_p, sep = "\t", row.names = FALSE, quote = FALSE, file = "Human_MD_D4.txt")
write.table(resdata, sep = "\t", row.names = FALSE, quote = FALSE, file = "Human_MD_D4_allgenes.txt")

# Volcano plot for DEGs
vdata <- arrange(resdata, padj)
vdata$threshold <- as.factor(ifelse(vdata$padj < 0.05 & abs(vdata$log2FoldChange) >= 1.5, ifelse(vdata$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
vdata$genetype[is.na(vdata$genetype)] <- "Not"
vdata$genetype <- factor(vdata$genetype, levels = c("MT_associated", "Not", "MT_coding", "PluripotencyGenes"))
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(padj)))
p1 <- p + geom_point(aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% c("MT_coding", "MT_associated", "PluripotencyGenes")),
             aes(shape = threshold, color = genetype), size = 1) +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#f8766d", "lightgray", "green", "#0033A0"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(-5, 160) + xlim(-10, 10)
ggsave(v4, file = "Human_RNA_volcano2.pdf", width = 5, height = 4)
# MT_associated Volcano plot
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(data = subset(vdata, !genetype %in% "MT_associated"), aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% "MT_associated"),
             aes(shape = threshold, color = genetype), size = 2, color = "#f8766d") +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("lightgray", "#7F7F7F", "#7F7F7F"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(v4, file = "Human_RNA_volcano_MT_associated.pdf", width = 8, height = 7)

# Pluripotency Volcano plot
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(data = subset(vdata, !genetype %in% "PluripotencyGenes"), aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% "PluripotencyGenes"),
             aes(shape = threshold, color = genetype), size = 2, color = "#0033A0") +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#7F7F7F", "lightgray", "#7F7F7F"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(v4, file = "Human_RNA_volcano_P.pdf", width = 8, height = 7)

# MT_coding Volcano plot
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(data = subset(vdata, !genetype %in% "MT_coding"), aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% "MT_coding"),
             aes(shape = threshold, color = genetype), size = 2, color = "green") +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#7F7F7F", "lightgray", "#7F7F7F"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(v4, file = "Human_RNA_volcano_MT_coding.pdf", width = 8, height = 7)

# Heatmap for differentially expressed genes
datamat <- as.matrix(diff_gene_p[, c(8:11)])
rownames(datamat) <- diff_gene_p$genename
anno_gene2 <- diff_gene_p[!is.na(diff_gene_p$genetype), c(12, 13)]
anno_gene <- as.factor(anno_gene2$genetype)
anno_gene <- as.matrix(anno_gene)
anno_gene <- as.data.frame(anno_gene)
rownames(anno_gene) <- anno_gene2$genename
rownames(datamat)[rownames(datamat) == "NA"] <- ""
# Create heatmap using ComplexHeatmap
library(ComplexHeatmap)
A <- as.matrix(datamat)  # Convert data matrix to matrix format
samples <- rep(c('MD_D4', 'WT_iPS'), c(2, 2))  # Define sample group information
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))  # Standardize data

# Create heatmap with annotations
B <- Heatmap(A,
             col = colorRampPalette(c("blue", "white", "red"))(100),
             show_row_names = FALSE,
             top_annotation = HeatmapAnnotation(Group = samples,
                                                simple_anno_size = unit(2, 'mm'),
                                                col = list(Group = c('MD_D4' = '#FF9289', 'WT_iPS' = '#00DAE0')),
                                                show_annotation_name = FALSE))

C <- B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "MT_associated"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "MT_associated"], labels_gp = gpar(fontsize = 4, col = "#f8766d")))
D <- C + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype != "MT_associated"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype != "MT_associated"], labels_gp = gpar(fontsize = 4, col = "green")))

# Save heatmap to PDF
pdf(file = "Human_DEGs7.pdf", width = 8, height = 10)
print(D)
dev.off()

# Set column names and convert to data frame
colnames(table_type) <- c("Celltype", "Genetype", "Genenumber")
table_sample_type <- as.data.frame(table_type)
table_sample_type$Genetype <- as.factor(table_sample_type$Genetype)
table_sample_type$Genenumber <- as.numeric(as.matrix(table_sample_type$Genenumber))

# Calculate percent weight for each cell type
ce <- ddply(table_sample_type, "Celltype", transform, percent_weight = Genenumber / sum(Genenumber) * 100)
ce$Genetype <- factor(ce$Genetype, levels = c("Up", "Not", "Down"))
ce$Celltype <- factor(ce$Celltype, levels = c("MT_associated", "PluripotencyGenes", "MT_coding"))
ce <- ddply(ce, "Celltype", transform, labely = cumsum(percent_weight))

# Create bar plot
p <- ggplot(ce, aes(x = Celltype, y = percent_weight, fill = Genetype)) + geom_bar(stat = "identity") + theme_bw()
p_bar2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#FF9289", "lightgray", "#00DAE0")) +
  scale_y_continuous(name = "Gene Percent", expand = c(0, 0)) + labs(x = "Gene Sets", colour = "") +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(y = labely, label = Genenumber), vjust = 1.5, colour = "black", size = 6) + coord_flip()
ggsave("Human_DEGs_Geneset_percent313_2.pdf", p_bar2, width = 6, height = 2)

# Venn diagram for MD and WT
library(eulerr)
MD = rownames(datad4[rowSums(datad4[, 1:2]) > 0, ])
WT = rownames(datad4[rowSums(datad4[, 3:4]) > 0, ])
vd <- euler(c(Mito = length(setdiff(MD, WT)), WT = length(setdiff(WT, MD)), "Mito&WT" = length(MD[MD %in% WT])))
pdf(file = "Human_venn.pdf", width = 5, height = 5)
plot(vd, fills = list(fill = c("#00DAE0", "#FF9289", "#7fb6b4"), alpha = 0.6), labels = list(col = "black", font = 10), edges = FALSE, quantities = TRUE)
dev.off()

# Venn diagram for DEGs
UP = rownames(resdata[resdata$threshold == "Up", ])
DOWN = rownames(resdata[resdata$threshold == "Down", ])
NOT = rownames(resdata[resdata$threshold == "Not", ])
vd2 <- euler(c(UP = length(UP), Down = length(DOWN), "UP&Down" = length(NOT)))
pdf(file = "Human_DEGs_venn.pdf", width = 5, height = 5)
plot(vd2, fills = list(fill = c("#FF9289", "#00DAE0", "lightgray"), alpha = 0.6), labels = list(col = "black", font = 15), edges = FALSE, quantities = TRUE)
dev.off()

# Gene Ontology (GO) and KEGG pathway enrichment for human DEGs
# GO enrichment for upregulated genes
diff_gene_p_up <- diff_gene_p$genename[diff_gene_p$stat == "Up"]
EG2symbol = toTable(org.Hs.egSYMBOL)
fid = as.character(diff_gene_p_up)
geneLists = data.frame(symbol = fid)
results_up = merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_up = na.omit(results_up$gene_id)
ego_up <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_up, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_up_result <- ego_up@result[ego_up@result$Count >= 10 & ego_up@result$pvalue < 0.05, ]
write.csv(ego_up_result, "Human_mito_depletion_deg_UPGO.csv")

# Draw GO enrichment bar plot for upregulated genes
ego_up_draw <- ego_up_result[1:10, ]
p <- ggplot(ego_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count") + coord_flip() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) +
  labs(title = "Human Mitochondria Depletion Up GO Annotation")
ggsave(p, filename = "Human_mito_depletion_deg_UPGO.pdf", height = 5, width = 10)

# KEGG pathway enrichment for upregulated genes
kegg_up <- enrichKEGG(organism = "hsa", gene = id_up, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_up_result <- kegg_up@result[kegg_up@result$Count >= 10 & kegg_up@result$pvalue < 0.05, ]
write.csv(kegg_up_result, "Human_mito_depletion_deg_UPKEGG.csv")

# Draw KEGG pathway bar plot for upregulated genes
kegg_up_draw <- kegg_up_result[1:10, ]
p <- ggplot(kegg_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count") + coord_flip() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) +
  labs(title = "Human Mitochondria Depletion Up KEGG Annotation")
ggsave(p, filename = "Human_mito_depletion_deg_UPKEGG.pdf", height = 5, width = 10)

# GO enrichment for downregulated genes
diff_gene_p_down <- diff_gene_p$genename[diff_gene_p$stat == "Down"]
diff_gene_p_down <- na.omit(diff_gene_p_down)
fid = as.character(diff_gene_p_down)
geneLists = data.frame(symbol = fid)
results_down = merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_down = na.omit(results_down$gene_id)
ego_down <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_down, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_down_result <- ego_down@result[ego_down@result$Count >= 10 & ego_down@result$pvalue < 0.05, ]
write.csv(ego_down_result, "Human_mito_depletion_deg_downGO.csv")

# Draw GO enrichment bar plot for downregulated genes
ego_down_draw <- ego_down_result[1:10, ]
p <- ggplot(ego_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count") + coord_flip() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) +
  labs(title = "Human Mitochondria Depletion Down GO Annotation")
ggsave(p, filename = "Human_mito_depletion_deg_downGO.pdf", height = 5, width = 10)

# KEGG pathway enrichment for downregulated genes
kegg_down <- enrichKEGG(organism = "hsa", gene = id_down, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down_result <- kegg_down@result[kegg_down@result$Count >= 10 & kegg_down@result$pvalue < 0.05, ]
write.csv(kegg_down_result, "Human_mito_depletion_deg_downKEGG.csv")

# Draw KEGG pathway bar plot for downregulated genes
kegg_down_draw <- kegg_down_result[1:10, ]
p <- ggplot(kegg_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count") + coord_flip() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) +
  labs(title = "Human Mitochondria Depletion Down KEGG Annotation")
ggsave(p, filename = "Human_mito_depletion_deg_downKEGG.pdf", height = 5, width = 10)
