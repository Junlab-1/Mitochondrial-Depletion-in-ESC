# Remove all objects
rm(list = ls(all = TRUE))

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(plyr)
library(eulerr)

# Set working directory
setwd("~/Human_MD_differentiation/Rresult")

# Define sample numbers
controlnum = 2
casenum = 2
num = casenum + controlnum

# Load data
dataall <- as.matrix(read.csv("../RNAdata/genematrix/gene_count_matrix.csv", row.names = "gene_id"))
dataall <- dataall[, c(4, 5, 2, 3)]

#-------------------------- Compare D7 with D14 -------------------------
sampletitle = "Human MT_neuro_D7 vs MT_neuro_D14"
countData <- dataall
condition <- factor(c(rep("treat", casenum), rep("untreat", controlnum)))
colData <- data.frame(row.names = colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize = TRUE)), by = "row.names", sort = FALSE)
resdata <- resdata[!is.na(resdata$padj), ]
resdata <- resdata[order(resdata$padj), ]

# Gene annotation
PluripotencyGene <- read.table("Pluripotency_Genes.txt", header = FALSE)
MT_coding <- read.table("Human mitochondria genome coding genes.csv", sep = ",", header = TRUE)
MT_asso <- read.table("human_chromosome_mito_genelist.csv", sep = ",", header = TRUE)
neurogene <- read.table("neuron_marker2.csv", sep = "\t", header = TRUE)
resdata$genename <- sapply(resdata$Row.names, function(x) strsplit(x, "\\|")[[1]][2])
resdata$genetype <- NA
resdata$genetype[resdata$genename %in% MT_coding[, 2]] <- "MT_coding"
resdata$genetype[resdata$genename %in% MT_asso[, 2]] <- "MT_associated"
resdata$genetype[resdata$genename %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"
resdata$genetype[resdata$genename %in% neurogene[, 1]] <- "NeuroDevelopment"
resdata$genename[is.na(resdata$genename)] <- resdata$Row.names[is.na(resdata$genename)]

# Volcano plot and heatmap for significant genes
diff_gene_p <- resdata[resdata$padj < 0.05, ]
diff_gene_p$stat[diff_gene_p$log2FoldChange < 0] <- "Up"
diff_gene_p$stat[diff_gene_p$log2FoldChange > 0] <- "Down"
diff_gene_p$log2FoldChange <- -diff_gene_p$log2FoldChange
diff_gene_p <- diff_gene_p[abs(diff_gene_p$log2FoldChange) > 1.5, ]
write.table(diff_gene_p, sep = "\t", row.names = FALSE, quote = FALSE, "Human_MT_neuro_D7vsD14.txt")

# Volcano plot
resdata$log2FoldChange <- -resdata$log2FoldChange
data <- resdata
data$threshold <- as.factor(ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1.5, ifelse(data$log2FoldChange >= 1.5, 'Up', 'Down'), 'Not'))
data$genetype[!data$Row.names %in% diff_gene_p$Row.names] <- NA
data$genename[is.na(data$genetype)] <- NA
data$pos <- rowMeans(data[, 8:(8 + casenum - 1)])
data$neg <- rowMeans(data[, (8 + casenum):(8 + num - 1)])
p <- ggplot(data, aes(x = log2(pos), y = log2(neg), colour = threshold)) + geom_point()
p1 <- p + geom_point(aes(color = threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha = 0.4, size = 1)
p3 <- p2 + theme_bw() + labs(x = "Case_genes", y = "Control_genes", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 0.6)
a <- cor.test(data$pos, data$neg, method = "pearson")
p6 <- p5 + annotate("text", x = 14, y = -0.5, label = paste("Pearson's correlation ", round(a$estimate, 3)), size = 3)
ggsave(p6, file = "Human_neurod7_vs_d14_RNA-seq_correlation2.pdf", width = 8, height = 8)
vdata <- arrange(data, padj)
vdata$threshold <- as.factor(ifelse(vdata$pvalue < 0.05 & abs(vdata$log2FoldChange) >= 1.5, ifelse(vdata$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
vdata$genetype[vdata$threshold == "Not"] <- "Not"
vdata$genetype <- factor(vdata$genetype, levels = c("MT_associated", "Not", "MT_coding", "PluripotencyGenes", "NeuroDevelopment"))
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% c("MT_coding", "MT_associated", "PluripotencyGenes", "NeuroDevelopment")),
             aes(shape = threshold, color = genetype), size = 1) +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) + geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#f8766d", "lightgray", "#7F7F7F", "#0033A0", "#7F7F7F"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlim(-12, 12) + ylim(-5, 105)
ggsave(v4, file = paste(sampletitle, "_volcano3.pdf", sep = ""), width = 5, height = 4)

# Heatmap
datamat <- as.matrix(diff_gene_p[, c(8:11)])
rownames(datamat) <- diff_gene_p$genename
anno_gene2 <- diff_gene_p[!is.na(diff_gene_p$genetype), c(12, 13)]
rownames(datamat)[rownames(datamat) == "NA"] <- ""
library(ComplexHeatmap)
A <- as.matrix(datamat)
samples <- rep(c('MD_neuro_D7', 'MD_neuro_D14'), c(2, 2))
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))
B <- Heatmap(A,
             col = colorRampPalette(c("blue", "white", "red"))(100),
             show_row_names = FALSE,
             top_annotation = HeatmapAnnotation(Group = samples,
                                                simple_anno_size = unit(2, 'mm'),
                                                col = list(Group = c('MD_neuro_D7' = '#F2AFEF', 'MD_neuro_D14' = '#9F70FD')),
                                                show_annotation_name = FALSE))
C <- B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "MT_associated"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "MT_associated"], labels_gp = gpar(fontsize = 4, col = "#f8766d")))
D <- C + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "MT_coding"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "MT_coding"], labels_gp = gpar(fontsize = 4, col = "green")))
E <- D + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "PluripotencyGenes"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "PluripotencyGenes"], labels_gp = gpar(fontsize = 4, col = "#0033A0")))
FF <- E + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "NeuroDevelopment"]),
                                         labels = anno_gene2$genename[anno_gene2$genetype == "NeuroDevelopment"], labels_gp = gpar(fontsize = 4, col = "#7F27FF")))

pdf(file = "Human_neuro7_14_2.pdf", width = 8, height = 10)
print(FF)
dev.off()

# Map gene percent
resdata$threshold <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1.5, ifelse(resdata$log2FoldChange > 0, 'Up', 'Down'), 'Not'))
write.table(resdata, sep = "\t", row.names = FALSE, quote = FALSE, "Human_MT_neuro_D7vsD14_all.txt")
merge_bar_all <- cbind(as.matrix(resdata$threshold), as.matrix(resdata$genetype))
for (i in c("PluripotencyGenes", "MT_associated", "MT_coding", "NeuroDevelopment")) {
  sign <- table(merge_bar_all[merge_bar_all[, 2] == i, 1])
  singletype_table <- cbind(rep(i, length(sign)), names(sign), as.matrix(sign))
  rownames(singletype_table) <- NULL
  if (i == "PluripotencyGenes") {
    table_type <- singletype_table
  } else {
    table_type <- rbind(table_type, singletype_table)
  }
}
colnames(table_type) <- c("Celltype", "Cancertype", "Cellnumber")
table_sample_type <- as.data.frame(table_type)
table_sample_type$Cancertype <- as.factor(table_sample_type$Cancertype)
table_sample_type$Cellnumber <- as.numeric(as.matrix(table_sample_type$Cellnumber))
ce <- ddply(table_sample_type, "Celltype", transform, percent_weight = Cellnumber / sum(Cellnumber) * 100)
ce$Cancertype <- factor(ce$Cancertype, levels = c("Up", "Not", "Down"))
ce$Celltype <- factor(ce$Celltype, levels = c("PluripotencyGenes", "MT_associated", "NeuroDevelopment", "MT_coding"))
ce <- ddply(ce, "Celltype", transform, labely = cumsum(percent_weight))
p <- ggplot(ce, aes(x = Celltype, y = percent_weight, fill = Cancertype)) + geom_bar(stat = "identity") + theme_bw()
p_bar2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#FF9289", "lightgray", "#00DAE0")) +
  scale_y_continuous(name = "Gene Percent", expand = c(0, 0)) + labs(x = "Gene Sets", colour = "") +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(y = labely, label = Cellnumber), vjust = 1.5, colour = "black", size = 6)
ggsave("Human_neuro_Geneset_percent313.pdf", p_bar2, width = 6, height = 6)


# Venn Diagram for MD_D7 and MD_D14
MD_D7 = rownames(dataall[rowSums(dataall[,1:2]) > 0, ])
MD_D14 = rownames(dataall[rowSums(dataall[,3:4]) > 0, ])
vd <- euler(c(MDD7 = length(setdiff(MD_D7, MD_D14)), MDD14 = length(setdiff(MD_D14, MD_D7)),
              "MDD7&MDD14" = length(MD_D7[MD_D7 %in% MD_D14])))
pdf(file = "Human_neuro_venn.pdf", width = 5, height = 5)
plot(vd,
     fills = list(fill = c("#F2AFEF", "#9F70FD", "lightgray"), alpha = 0.6),
     labels = list(col = "black", font = 10), 
     edges = FALSE,
     quantities = TRUE)
dev.off()

# Venn Diagram for DEGs
UP = rownames(resdata[resdata$threshold == "Up", ])
DOWN = rownames(resdata[resdata$threshold == "Down", ])
NOT = rownames(resdata[resdata$threshold == "Not", ])
vd2 <- euler(c(UP = length(UP), Down = length(DOWN),
               "UP&Down" = length(NOT)))
pdf(file = "Human_DEGs_venn.pdf", width = 5, height = 5)
plot(vd2,
     fills = list(fill = c("#FF9289", "#00DAE0", "lightgray"), alpha = 0.6),
     labels = list(col = "black", font = 15), 
     edges = FALSE,
     quantities = TRUE)
dev.off()

# Gene Function Annotation for DEGs
# Upregulated genes
diff_gene_p_up <- diff_gene_p$genename[diff_gene_p$stat == "Up"]
EG2symbol <- toTable(org.Hs.egSYMBOL)
fid <- as.character(diff_gene_p_up)
geneLists <- data.frame(symbol = fid)
results_up <- merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_up <- na.omit(results_up$gene_id)
ego_up <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_up, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_up_result <- ego_up@result[ego_up@result$Count >= 10 & ego_up@result$pvalue < 0.05, ]
write.csv(ego_up_result, "Human_mito_neuro_d7d14_deg_UPGO.csv")

# Draw GO enrichment bar plot for upregulated genes
ego_up_draw <- ego_up_result[1:10, ]
p <- ggplot(ego_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(ego_up_draw$Count) + 1))) + coord_flip() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria Up GO Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d14_deg_UPGO.pdf", height = 5, width = 10)

# KEGG enrichment for upregulated genes
kegg_up <- enrichKEGG(organism = "hsa", gene = id_up, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_up_result <- kegg_up@result[kegg_up@result$Count >= 10 & kegg_up@result$pvalue < 0.05, ]
write.csv(kegg_up_result, "Human_mito_neuro_d7d14_deg_UPKEGG.csv")

# Draw KEGG enrichment bar plot for upregulated genes
kegg_up_draw <- kegg_up_result[1:10, ]
p <- ggplot(kegg_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(kegg_up_draw$Count) + 1))) + coord_flip() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria neuro_d7d14 Up KEGG Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d14_deg_UPKEGG.pdf", height = 5, width = 10)

# Downregulated genes
diff_gene_p_down <- diff_gene_p$genename[diff_gene_p$stat == "Down"]
diff_gene_p_down <- na.omit(diff_gene_p_down)
fid <- as.character(diff_gene_p_down)
geneLists <- data.frame(symbol = fid)
results_down <- merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_down <- na.omit(results_down$gene_id)
ego_down <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_down, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_down_result <- ego_down@result[ego_down@result$Count >= 10 & ego_down@result$pvalue < 0.05, ]
write.csv(ego_down_result, "Human_mito_neuro_d7d14_deg_downGO.csv")

# Draw GO enrichment bar plot for downregulated genes
ego_down_draw <- ego_down_result[1:10, ]
p <- ggplot(ego_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(ego_down_draw$Count) + 1))) + coord_flip() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria neuro_d7d14 Down GO Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d14_deg_downGO.pdf", height = 5, width = 10)

# KEGG enrichment for downregulated genes
kegg_down <- enrichKEGG(organism = "hsa", gene = id_down, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down_result <- kegg_down@result[kegg_down@result$Count >= 10 & kegg_down@result$pvalue < 0.05, ]
write.csv(kegg_down_result, "Human_mito_neuro_d7d14_deg_downKEGG.csv")

# Draw KEGG enrichment bar plot for downregulated genes
kegg_down_draw <- kegg_down_result[1:10, ]
p <- ggplot(kegg_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(kegg_down_draw$Count) + 1))) + coord_flip() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria neuro_d7d14 Down KEGG Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d14_deg_downKEGG.pdf", height = 5, width = 10)

##-------------------- merge D1 D7 D14 ----------------------------------
data0<-as.matrix(read.csv("~/Human_Mito-depletion/RNAdata/genematrix/gene_count_matrix.csv", row.names="gene_id"))
data0<-data0[,c(3,4,5,6)]
data0_common   <- data0[rownames(data0)%in%rownames(dataall),]
dataall_common <- dataall[rownames(dataall)%in%rownames(data0),]
data0_common   <- data0_common[order(rownames(data0_common)),]
dataall_common <- dataall_common[order(rownames(dataall_common)),]
if (identical(rownames(data0_common),rownames(dataall_common))) {
  dataall2<-cbind(data0_common,dataall_common)
}else{
  print("Genes are not same!")
}
write.table(file = "M_D_hiPSC_counts.txt",x = dataall2,sep = "\t")
dataall2<-dataall2[,-c(1,2)]
colnames(dataall2)[1:2]<-c("MD_D1_1","MD_D1_2")

##--------------------- compare D1 vs D7 --------------------------------

# Prepare data for D1 vs D7 comparison
dataall3 <- dataall2[, c(3, 4, 1, 2)]
setwd("./D7vsD1")
sampletitle <- "Human MT_neuro_D7 vs MT_neuro_D1"
countData <- dataall3
condition <- factor(c(rep("treat", casenum), rep("untreat", controlnum)))
colData <- data.frame(row.names = colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize = TRUE)), by = "row.names", sort = FALSE)
resdata <- resdata[!is.na(resdata$padj), ]
resdata <- resdata[order(resdata$padj), ]

# Annotate genes
resdata$genename <- sapply(resdata$Row.names, function(x) strsplit(x, "\\|")[[1]][2])
resdata$genetype <- NA
resdata$genetype[resdata$genename %in% MT_coding[, 2]] <- "MT_coding"
resdata$genetype[resdata$genename %in% MT_asso[, 2]] <- "MT_associated"
resdata$genetype[resdata$genename %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"
resdata$genetype[resdata$genename %in% neurogene[, 1]] <- "NeuroDevelopment"
resdata$genename[is.na(resdata$genename)] <- resdata$Row.names[is.na(resdata$genename)]
resdata$log2FoldChange <- -resdata$log2FoldChange

# Identify differentially expressed genes
diff_gene_p <- resdata[resdata$padj < 0.05, ]
diff_gene_p$stat[diff_gene_p$log2FoldChange > 0] <- "Up"
diff_gene_p$stat[diff_gene_p$log2FoldChange < 0] <- "Down"
diff_gene_p <- diff_gene_p[abs(diff_gene_p$log2FoldChange) > 1.5, ]
write.table(diff_gene_p, sep = "\t", row.names = FALSE, quote = FALSE, "Human_MT_neuro_D7vsD1.txt")

# Generate volcano plot
data <- resdata
data$threshold <- as.factor(ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1.5, ifelse(data$log2FoldChange >= 1.5, 'Up', 'Down'), 'Not'))
data$genetype[!data$Row.names %in% diff_gene_p$Row.names] <- NA
data$genename[is.na(data$genetype)] <- NA
data$pos <- rowMeans(data[, 8:(8 + casenum - 1)])
data$neg <- rowMeans(data[, (8 + casenum):(8 + num - 1)])
p <- ggplot(data, aes(x = log2(pos), y = log2(neg), colour = threshold)) + geom_point()
p1 <- p + geom_point(aes(color = threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha = 0.4, size = 1)
p3 <- p2 + theme_bw() + labs(x = "Case_genes", y = "Control_genes", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 0.6)
a <- cor.test(data$pos, data$neg, method = "pearson")
p6 <- p5 + annotate("text", x = 14, y = -0.5, label = paste("Pearson's correlation ", round(a$estimate, 3)), size = 3)
ggsave(p6, file = "Human_neurod7_vs_d0_RNA-seq_correlation2.pdf", width = 8, height = 8)

# Generate detailed volcano plot
vdata <- arrange(data, padj)
vdata$threshold <- as.factor(ifelse(vdata$pvalue < 0.05 & abs(vdata$log2FoldChange) >= 1.5, ifelse(vdata$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
vdata$genetype[vdata$threshold == "Not"] <- "Not"
vdata$genetype <- factor(vdata$genetype, levels = c("MT_associated", "Not", "MT_coding", "PluripotencyGenes", "NeuroDevelopment"))
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(pvalue)))
p1 <- p + geom_point(aes(shape = threshold, color = genetype), size = 1) +
  geom_point(data = subset(vdata, genetype %in% c("MT_coding", "MT_associated", "PluripotencyGenes", "NeuroDevelopment")),
             aes(shape = threshold, color = genetype), size = 1) +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#f8766d", "lightgray", "green", "#0033A0", "#7F27FF"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlim(-15, 15) + ylim(-5, 200)
ggsave(v4, file = paste(sampletitle, "_volcano2.pdf", sep = ""), width = 5, height = 4)

# Generate heatmap
datamat <- as.matrix(diff_gene_p[, c(8:11)])
rownames(datamat) <- diff_gene_p$genename
anno_gene2 <- diff_gene_p[!is.na(diff_gene_p$genetype), c(12, 13)]
anno_gene <- as.factor(anno_gene2$genetype)
anno_gene <- as.matrix(anno_gene)
anno_gene <- as.data.frame(anno_gene)
rownames(anno_gene) <- anno_gene2$genename
rownames(datamat)[rownames(datamat) == "NA"] <- ""

anno_gene2 <- anno_gene2[-sample(which(anno_gene2$genetype == "MT_associated"), ceiling(sum(anno_gene2$genetype == "MT_associated") / 4) * 3), ]

# Create heatmap using ComplexHeatmap
A <- as.matrix(datamat)
samples <- rep(c('MD_neuro_D7', 'MD_neuro_D1'), c(2, 2))
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))
B <- Heatmap(A,
             col = colorRampPalette(c("blue", "white", "red"))(100),
             show_row_names = FALSE,
             top_annotation = HeatmapAnnotation(Group = samples, 
                                                simple_anno_size = unit(2, 'mm'), 
                                                col = list(Group = c('MD_neuro_D7' = '#F2AFEF', 'MD_neuro_D1' = '#FF9289')),
                                                show_annotation_name = FALSE))
C <- B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "MT_associated"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "MT_associated"], labels_gp = gpar(fontsize = 4, col = "#f8766d")))
D <- C + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% anno_gene2$genename[anno_gene2$genetype == "MT_coding"]),
                                        labels = anno_gene2$genename[anno_gene2$genetype == "MT_coding"], labels_gp = gpar(fontsize = 4, col = "green")))
FF <- D + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% c("SOX2", "NANOG", "POU5F1")),
                                         labels = c("SOX2", "NANOG", "POU5F1"), labels_gp = gpar(fontsize = 6, col = "#01204E")))
pdf(file = "Human_neuro7_04.pdf", width = 8, height = 10)
print(FF)
dev.off()

# Generate gene set percentage bar plot
resdata$threshold <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1.5, ifelse(resdata$log2FoldChange > 0, 'Up', 'Down'), 'Not'))
merge_bar_all <- cbind(as.matrix(resdata$threshold), as.matrix(resdata$genetype))
for (i in c("PluripotencyGenes", "MT_associated", "MT_coding", "NeuroDevelopment")) {
  sign <- table(merge_bar_all[merge_bar_all[, 2] == i, 1])
  singletype_table <- cbind(rep(i, length(sign)), names(sign), as.matrix(sign))
  rownames(singletype_table) <- NULL
  if (i == "PluripotencyGenes") {
    table_type <- singletype_table
  } else {
    table_type <- rbind(table_type, singletype_table)
  }
}

colnames(table_type) <- c("Celltype", "Cancertype", "Cellnumber")
table_sample_type <- as.data.frame(table_type)
table_sample_type$Cancertype <- as.factor(table_sample_type$Cancertype)
table_sample_type$Cellnumber <- as.numeric(as.matrix(table_sample_type$Cellnumber))

ce <- ddply(table_sample_type, "Celltype", transform, percent_weight = Cellnumber / sum(Cellnumber) * 100)
ce$Cancertype <- factor(ce$Cancertype, levels = c("Up", "Not", "Down"))
ce$Celltype <- factor(ce$Celltype, levels = c("MT_coding", "PluripotencyGenes", "NeuroDevelopment", "MT_associated"))
ce <- ddply(ce, "Celltype", transform, labely = cumsum(percent_weight))
p <- ggplot(ce, aes(x = Celltype, y = percent_weight, fill = Cancertype)) + geom_bar(stat = "identity", width = 0.7) + theme_bw()
p_bar2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#FF9289", "lightgray", "#00DAE0")) +
  scale_y_continuous(name = "Gene Percent", expand = c(0, 0)) + labs(x = "Gene Sets", colour = "") +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(y = labely, label = Cellnumber), vjust = 1.5, colour = "black", size = 6)
ggsave("Human_neuro_Geneset_percent2.pdf", p_bar2, width = 5, height = 4)

# Generate Venn diagrams
MD_D7 <- rownames(dataall3[rowSums(dataall3[, 1:2]) > 0, ])
MD_D1 <- rownames(dataall3[rowSums(dataall3[, 3:4]) > 0, ])
vd <- euler(c(MDD7 = length(setdiff(MD_D7, MD_D1)), MDD1 = length(setdiff(MD_D1, MD_D7)), "MDD7&MDD1" = length(MD_D7[MD_D7 %in% MD_D1])))
pdf(file = "Human_neuro_venn2.pdf", width = 5, height = 5)
plot(vd,
     fills = list(fill = c("#F2AFEF", "#FF9289", "lightgray"), alpha = 0.6),
     labels = list(col = "black", font = 10),
     edges = FALSE,
     quantities = TRUE)
dev.off()

# Generate DEG Venn diagram
UP <- rownames(resdata[resdata$threshold == "Up", ])
DOWN <- rownames(resdata[resdata$threshold == "Down", ])
NOT <- rownames(resdata[resdata$threshold == "Not", ])
vd2 <- euler(c(UP = length(UP), Down = length(DOWN), "UP&Down" = length(NOT)))
pdf(file = "Human_DEGs_venn2.pdf", width = 5, height = 5)
plot(vd2,
     fills = list(fill = c("#FF9289", "#00DAE0", "lightgray"), alpha = 0.6),
     labels = list(col = "black", font = 15),
     edges = FALSE,
     quantities = TRUE)
dev.off()

# Perform GO and KEGG enrichment analysis
diff_gene_p_up <- diff_gene_p$genename[diff_gene_p$stat == "Up"]
EG2symbol <- toTable(org.Hs.egSYMBOL)
fid <- as.character(diff_gene_p_up)
geneLists <- data.frame(symbol = fid)
results_up <- merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_up <- na.omit(results_up$gene_id)
ego_up <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_up, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_up_result <- ego_up@result[ego_up@result$Count >= 10 & ego_up@result$pvalue < 0.05, ]
write.csv(ego_up_result, "Human_mito_neuro_d7d0_deg_UPGO2.csv")

# Draw up GO bar plot
ego_up_draw <- ego_up_result[1:10, ]
p <- ggplot(ego_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7))
p <- p + xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(ego_up_draw$Count) + 1))) + coord_flip()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes")
p <- p + theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria Up GO Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d0_deg_UPGO2.pdf", height = 5, width = 10)

kegg_up <- enrichKEGG(organism = "hsa", gene = id_up, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_up_result <- kegg_up@result[kegg_up@result$Count >= 10 & kegg_up@result$pvalue < 0.05, ]
write.csv(kegg_up_result, "Human_mito_neuro_d7d0_deg_UPKEGG2.csv")

# Draw up KEGG bar plot
kegg_up_draw <- kegg_up_result[1:10, ]
p <- ggplot(kegg_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7))
p <- p + xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(kegg_up_draw$Count) + 1))) + coord_flip()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes")
p <- p + theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria neuro_d7d0 Up KEGG Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d0_deg_UPKEGG2.pdf", height = 5, width = 10)

# Down GO enrichment analysis
diff_gene_p_down <- diff_gene_p$genename[diff_gene_p$stat == "Down"]
diff_gene_p_down <- na.omit(diff_gene_p_down)
fid <- as.character(diff_gene_p_down)
geneLists <- data.frame(symbol = fid)
results_down <- merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_down <- na.omit(results_down$gene_id)
ego_down <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_down, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_down_result <- ego_down@result[ego_down@result$Count >= 10 & ego_down@result$pvalue < 0.05, ]
write.csv(ego_down_result, "Human_mito_neuro_d7d0_deg_downGO2.csv")

# Draw down GO bar plot
ego_down_draw <- ego_down_result[1:10, ]
p <- ggplot(ego_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7))
p <- p + xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(ego_down_draw$Count) + 1))) + coord_flip()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes")
p <- p + theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria neuro_d7d0 Down GO Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d0_deg_downGO2.pdf", height = 5, width = 10)

# Down KEGG enrichment analysis
kegg_down <- enrichKEGG(organism = "hsa", gene = id_down, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down_result <- kegg_down@result[kegg_down@result$Count >= 10 & kegg_down@result$pvalue < 0.05, ]
write.csv(kegg_down_result, "Human_mito_neuro_d7d0_deg_downKEGG2.csv")

# Draw down KEGG bar plot
kegg_down_draw <- kegg_down_result[1:10, ]
p <- ggplot(kegg_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7))
p <- p + xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (max(kegg_down_draw$Count) + 1))) + coord_flip()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes")
p <- p + theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria neuro_d7d0 Down KEGG Annotation")
ggsave(p, filename = "Human_mito_neuro_d7d0_deg_downKEGG2.pdf", height = 5, width = 10)

##-------------------Show all data: Control, D1, D7, and D14 ------------------------------
# standard for all data
countData <- dataall2
condition <- factor(c(rep("day0",2),rep("day7",2),rep("day14",2)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
resdata <- resdata[!is.na(resdata$padj),]
resdata<-resdata[order(resdata$padj),]
# gene annotation
resdata$genename=sapply(resdata$Row.names,function(x) strsplit(x, "\\|")[[1]][2])
resdata$genetype=NA
resdata$genetype[resdata$genename%in%MT_coding[,2]]<-"MT_coding"
resdata$genetype[resdata$genename%in%MT_asso[,2]]<-"MT_associated"
resdata$genetype[resdata$genename%in%PluripotencyGene[,1]]<-"PluripotencyGenes"
resdata$genetype[resdata$genename%in%neurogene[,1]]<-"NeuroDevelopment"
resdata$genename[is.na(resdata$genename)]<-resdata$Row.names[is.na(resdata$genename)]

# 1. pca
library(ggplot2)
library(ggbiplot)
# Transpose and prepare data for PCA
pcamat <- as.data.frame(t(resdata[, 8:13]))
rownames(pcamat) <- c("MD_D1_S1", "MD_D1_S2", "MD_D7_S1", "MD_D7_S2", "MD_D14_S1", "MD_D14_S2")
pcamatgroup <- c(rep("MD_D1", 2), rep("MD_D7", 2), rep("MD_D14", 2))
# Perform PCA with scaling and centering
pca_result <- prcomp(pcamat, scale = TRUE, center = TRUE)
# Plot PCA results and save to PDF
pdf(file = "Human_neuro_pcaplot.pdf", width = 6, height = 6)
ggbiplot(pca_result,
         var.axes = FALSE,            # Do not plot variable arrows
         obs.scale = 0.4,             # Scaling for the observations
         groups = pcamatgroup,        # Add grouping information
         ellipse = FALSE,
         circle = FALSE) +
  geom_text(aes(label = rownames(pcamat)), vjust = 1.5, size = 2) + # Add labels
  theme_bw() + 
  geom_point(aes(color = pcamatgroup)) + 
  scale_color_manual(values = c("MD_D1" = "#FF9289", "MD_D7" = "#F2AFEF", "MD_D14" = "#9F70FD"))
dev.off()

# 2. Gene Trend Color Merge (Upset Plot)
library(UpSetR)
# Load DEGs data
deg71 <- read.table("../D7vsD1/Human_MT_neuro_D7vsD1.txt", header = TRUE, sep = "\t")
deg147 <- read.table("../Human_MT_neuro_D7vsD14.txt", header = TRUE, sep = "\t")
# Find common DEGs with the same and opposite trends
deg71_up <- deg71$genename[deg71$stat == "Up"]
deg71_down <- deg71$genename[deg71$stat == "Down"]
deg147_up <- deg147$genename[deg147$stat == "Down"]
deg147_down <- deg147$genename[deg147$stat == "Up"]
# Find intersections
UU <- intersect(deg71_up, deg147_up)
DD <- intersect(deg71_down, deg147_down)
U71_D147 <- intersect(deg71_up, deg147_down)
D71_U147 <- intersect(deg71_down, deg147_up)
# Prepare data for UpSet plot
dataForUpSetPlot <- list(D7_D1_up = deg71_up, D7_D1_down = deg71_down,
                         D14_D7_up = deg147_up, D14_D7_down = deg147_down)
setsBarColors <- c('#00DAE0', '#FF9289', '#4285F4', '#EA4335')

# Plot UpSet plot and save to PDF
p <- upset(fromList(dataForUpSetPlot),
           nsets = length(dataForUpSetPlot),
           nintersects = 1000,
           sets = c("D14_D7_down", "D14_D7_up", 'D7_D1_down', "D7_D1_up"),
           keep.order = TRUE,
           point.size = 3,
           line.size = 1,
           number.angles = 0,
           text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1), # ytitle, ylabel, xtitle, xlabel, sets, number
           order.by = "freq",
           matrix.color = "black",
           main.bar.color = 'black',
           mainbar.y.label = 'Intersection Size',
           sets.bar.color = setsBarColors)
pdf(file = "Human_d1471_mergedegs.pdf", width = 8, height = 6)
print(p)
dev.off()

# 3. Gene Trend Analysis
library(reshape2)
# Prepare data for trend analysis
commonall <- dataall2
commonall[, 1] <- rowMeans(commonall[, 1:2])
commonall[, 3] <- rowMeans(commonall[, 3:4])
commonall[, 5] <- rowMeans(commonall[, 5:6])
commonall <- commonall[, -c(2, 4, 6)]
colnames(commonall) <- paste("MD_D", c(0, 7, 14), sep = "")
commonall <- as.data.frame(commonall)
commonall$genename <- sapply(rownames(commonall), function(x) strsplit(x, "\\|")[[1]][2])
commonall$genetype <- "NOT"
commonall$genetype[commonall$genename %in% MT_coding[, 2]] <- "MT_coding"
commonall$genetype[commonall$genename %in% MT_asso[, 2]] <- "MT_associated"
commonall$genetype[commonall$genename %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"
commonall$genetype[commonall$genename %in% neurogene[, 1]] <- "NeuroDevelopment"
commonall$genename[is.na(commonall$genename)] <- rownames(commonall)[is.na(commonall$genename)]
annoall <- commonall[!commonall$genetype == "NOT", ]
annoall2 <- annoall[annoall$genename %in% deg71$genename | annoall$genename %in% deg147$genename, ]
# Prepare colors for plotting
colorbar <- data.frame(name = c("MT_associated", "MT_coding", "NeuroDevelopment", "PluripotencyGenes"),
                       color = c("#f8766d", "green", "#7F27FF", "#0033A0"))
# Plot gene trends for each category and save to PDF
for (i in colorbar$name) {
  imat <- annoall2[annoall2$genetype == i, ]
  imat[, 1:3] <- log10(imat[, 1:3] + 1)
  imat[, 1:3] <- imat[, 1:3] - imat[, 1]
  drawmat <- melt(imat)
  p <- ggplot(drawmat, aes(x = variable, y = value, colour = genetype, group = genename)) + 
       geom_line(alpha = 0.5) + geom_point(size = 2, alpha = 0.5)
  p <- p + scale_color_manual(values = c("MT_associated" = "#f8766d", "MT_coding" = "green",
                                         "NeuroDevelopment" = "#7F27FF", "PluripotencyGenes" = "#0033A0"))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
       scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
       labs(x = "Sample Sets", y = "Standard Gene Exp", colour = "")
  ggsave(p, file = paste(i, "_line.pdf"), width = 4, height = 5)
}

#Function annotation bar plot
library("org.Hs.eg.db")
genesets<-c("UU","DD","U71_D147","D71_U147")
EG2symbol=toTable(org.Hs.egSYMBOL)
EG2symbol$symbol<-toupper(EG2symbol$symbol)

for (i in genesets[1:4]) {
  eval(parse(text = paste("iset<-",i,sep = "")))
  print(i)
  # KEGG、GO
  fid=as.character(iset)
  geneLists=data.frame(symbol=fid)
  results_up=merge(geneLists,EG2symbol,by='symbol',all.x=T)
  id_up=na.omit(results_up$gene_id)
  ego_up <- enrichGO(OrgDb="org.Hs.eg.db", gene = id_up, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 3, readable= TRUE)
  ego_up_result<-ego_up@result[ego_up@result$Count>=3&ego_up@result$pvalue<0.05,]
  write.csv(ego_up_result,paste(i,"_deg_GO.csv",sep = ""))
  ego_up_draw<-ego_up_result[1:10,]
  p<-ggplot(ego_up_draw,aes(x=reorder(Description,Count),y=Count,fill=Count))+
    geom_bar(stat = "identity",colour="NA",width = 0.9,position = position_dodge(0.7))
  p<-p+xlab("")+scale_y_continuous(name = "Gene Count",expand = c(0,0),limits=c(0,(max(ego_up_draw$Count)+1)))+coord_flip()
  p<-p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))
  p<-p+scale_fill_gradient(low="#FF9289",high="#FF9289",name="Genes")
  p<-p+theme(text=element_text(size=16, family="serif"))+labs(title = i)
  ggsave(p,filename = paste(i,"_deg_GO.pdf",sep = ""),height = 5,width = 12)
  kegg_up<-enrichKEGG(organism = "hsa",gene=id_up,keyType = "kegg",pvalueCutoff = 0.05,minGSSize = 3)
  kegg_up<-setReadable(kegg_up,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
  kegg_up_result<-kegg_up@result[kegg_up@result$Count>=3&kegg_up@result$pvalue<0.05,]
  write.csv(kegg_up_result,paste(i,"_deg_KEGG.csv",sep = ""))
  kegg_up_draw<-kegg_up_result[1:10,]
  p<-ggplot(kegg_up_draw,aes(x=reorder(Description,Count),y=Count,fill=Count))+
    geom_bar(stat = "identity",colour="NA",width = 0.9,position = position_dodge(0.7))
  p<-p+xlab("")+scale_y_continuous(name = "Gene Count",expand = c(0,0),limits=c(0,(max(kegg_up_draw$Count)+1)))+coord_flip()
  p<-p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))
  p<-p+scale_fill_gradient(low="#FF9289",high="#FF9289",name="Genes")
  p<-p+theme(text=element_text(size=16, family="serif"))+labs(title = i)
  ggsave(p,filename = paste(i,"_deg_KEGG.pdf",sep = ""),height = 5,width = 10)
}
