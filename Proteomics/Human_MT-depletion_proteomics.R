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
library(ggbiplot)
library(ComplexHeatmap)
library(eulerr)
library(plyr)

# Set working directory
setwd("~/Human proteomic")

# Load and clean data
human_all <- read.table(file = "human_filtered_norm.txt", header = TRUE, sep = "\t")
human_all <- na.omit(human_all)
# Remove rows with more than 2 NAs in columns 3 to 6
rmlist <- 0
for (i in 1:nrow(human_all)) {
  if (sum(is.na(human_all[i, 3:6])) > 2) {
    if (sum(rmlist) == 0) {
      rmlist <- i
    } else {
      rmlist <- c(rmlist, i)
    }
  }
}
if (rmlist == 0) {
  human_all2_org <- human_all
} else {
  human_all2_org <- human_all[-rmlist, ]
}

# Normalize data
human_all2 <- human_all2_org
human_all2[, 3:6] <- log2(human_all2_org[, 3:6])
human_all2[, 3:6][human_all2[, 3:6] == -Inf] <- 0

for (i in 3:6) {
  human_all2[, i] <- human_all2[, i] - mean(human_all2[, i])
}

# Standardize data
human_all2$avg <- rowMeans(human_all2[, 3:6])
human_all3 <- human_all2
for (i in 3:6) {
  imodel <- lm(human_all2$avg ~ human_all2[, i])
  islope <- coef(imodel)[2]
  human_all3[, i] <- human_all2[, i] / islope
}

# Calculate DEPs
for (i in 1:nrow(human_all3)) {
  iresult <- t.test(human_all3[i, 3:4], human_all3[i, 5:6])
  human_all3$pvalue[i] <- iresult$p.value
  if (human_all3[i, 1] == human_all2_org[i, 1]) {
    human_all3$log2fc[i] <- log2(mean(as.numeric(human_all2_org[i, 5:6])) / mean(as.numeric(human_all2_org[i, 3:4])))
  }
}

# Filter DEPs
PluripotencyGene <- read.table("Pluripotency_Genes.txt", header = FALSE)
MT_coding <- read.table("Human mitochondria genome coding genes.csv", sep = ",", header = TRUE)
MT_asso <- read.table("human_chromosome_mito_genelist.csv", sep = ",", header = TRUE)
human_dep2 <- human_all3[abs(human_all3$log2fc) > 1.5 & human_all3$pvalue < 0.05, ]
human_all3$proteintype <- NA
human_all3$proteintype[human_all3$Gene.Symbol %in% MT_coding[, 2]] <- "MT_coding"
human_all3$proteintype[human_all3$Gene.Symbol %in% MT_asso[, 2]] <- "MT_associated"
human_all3$proteintype[human_all3$Gene.Symbol %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"
write.csv(human_all3, "human_all_proteins.csv")

# Annotate DEPs
human_dep2$states <- "Up"
human_dep2$states[human_dep2$log2fc < 0] <- "Down"
a <- human_all2_org$Accession[human_all2_org$Accession %in% human_dep2$Accession]
if (identical(human_dep2$Accession, a)) {
  human_dep2[, 3:6] <- human_all2_org[human_all2_org$Accession %in% human_dep2$Accession, 3:6]
}
human_dep2$avg <- NULL
human_dep2$proteintype <- NA
human_dep2$proteintype[human_dep2$Gene.Symbol %in% MT_coding[, 2]] <- "MT_coding"
human_dep2$proteintype[human_dep2$Gene.Symbol %in% MT_asso[, 2]] <- "MT_associated"
human_dep2$proteintype[human_dep2$Gene.Symbol %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"
write.csv(human_dep2, "human_dep2.csv")

# Volcano plot
sampletitle <- "Human Mito vs H9 Proteome"
data <- human_all3
if (identical(human_all2_org$Accession, data$Accession)) {
  data[, 3:6] <- human_all2_org[, 3:6]
}
data$threshold <- as.factor(ifelse(data$pvalue < 0.05 & abs(data$log2fc) >= 1.5, ifelse(data$log2fc > 0, 'Up', 'Down'), 'Not'))
data$neg <- rowMeans(data[, 3:4])
data$pos <- rowMeans(data[, 5:6])
colnames(data)[3:6] <- c("WT_sample1", "WT_sample2", "Mito_Depletion1", "Mito_Depletion2")

# Scatter plot with correlation
p <- ggplot(data, aes(x = log2(pos), y = log2(neg), colour = threshold)) + geom_point()
p1 <- p + geom_point(aes(color = threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha = 0.4, size = 1)
p3 <- p2 + theme_bw() + labs(x = "Case_genes", y = "Control_genes", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 0.6)
a <- cor.test(data$pos, data$neg, method = "pearson")
p6 <- p5 + annotate("text", x = 25, y = 10, label = paste("Pearson's correlation ", round(a$estimate, 3)), size = 3)
ggsave(p6, file = "Human_protein_correlation.pdf", width = 9, height = 8)

# Volcano plot for DEPs
vdata <- arrange(data, pvalue)
vdata$threshold <- as.factor(ifelse(vdata$pvalue < 0.05 & abs(vdata$log2fc) >= 1.5, ifelse(vdata$log2fc > 0, 'Up', 'Down'), 'Not'))
vdata$proteintype[vdata$threshold == "Not"] <- "Not"
vdata$proteintype <- factor(vdata$proteintype, levels = c("MT_associated", "Not", "MT_coding", "PluripotencyGenes"))

p <- ggplot(vdata, aes(x = log2fc, y = -log10(pvalue)))
p1 <- p + geom_point(aes(shape = threshold, color = proteintype), size = 1) +
  geom_point(data = subset(vdata, proteintype %in% c("MT_coding", "MT_associated", "PluripotencyGenes")),
             aes(shape = threshold, color = proteintype), size = 1) +
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  scale_x_continuous(limits = c(-7, 7), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#f8766d", "lightgray", "#0033A0"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(v4, file = "Human_protein_volcano.pdf", width = 5, height = 4)

# PCA plot
pcamat <- as.data.frame(t(human_all2_org[, 3:6]))
rownames(pcamat) <- c("WT_Day0_S1", "WT_Day0_S2", "MT_Depletion_S1", "MT_Depletion_S2")
pcamatgroup <- c(rep("WT", 2), rep("MT_Depletion", 2))
pcamat <- pcamat[, colSums(pcamat) != 0]
pca_result <- prcomp(pcamat, scale = FALSE, center = TRUE)
pca_result$x <- scale(pca_result$x)

# Plot PCA results
pdf(file = "Human_protein_MT_Depletion_pcaplot3.pdf", width = 4, height = 6)
ggbiplot(pca_result,
         var.axes = FALSE,
         groups = pcamatgroup,
         ellipse = FALSE,
         circle = FALSE) +
  geom_text(aes(label = rownames(pcamat)), vjust = 1.5, size = 2) +
  theme_bw()
dev.off()

# Heatmap
datamat <- as.matrix(human_dep2[, 3:6])
rownames(datamat) <- human_dep2$Gene.Symbol
human_dep2$proteintype[is.na(human_dep2$proteintype)] <- "None"
A <- as.matrix(datamat)
samples <- rep(c('WT_protein', 'MD_D4'), c(2, 2))
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))

B <- Heatmap(A,
             col = colorRampPalette(c("blue", "white", "red"))(100),
             show_row_names = FALSE,
             top_annotation = HeatmapAnnotation(Group = samples,
                                                simple_anno_size = unit(2, 'mm'),
                                                col = list(Group = c('MD_D4' = '#FF9289', 'WT_protein' = '#00DAE0')),
                                                show_annotation_name = FALSE))

C <- B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% human_dep2$Gene.Symbol[human_dep2$proteintype %in% c("MT_associated", "PluripotencyGenes")]),
                                        labels = human_dep2$Gene.Symbol[human_dep2$proteintype %in% c("MT_associated", "PluripotencyGenes")], labels_gp = gpar(fontsize = 4, col = "#f8766d")))

pdf(file = "Human_DEP_heatmap2.pdf", width = 8, height = 15)
print(C)
dev.off()

# Gene percentage bar plot
human_all3$threshold <- as.factor(ifelse(human_all3$pvalue < 0.05 & abs(human_all3$log2fc) >= 1.5, ifelse(human_all3$log2fc > 0, 'Up', 'Down'), 'Not'))
merge_bar_all <- cbind(as.matrix(human_all3$threshold), as.matrix(human_all3$proteintype))

for (i in c("PluripotencyGenes", "MT_associated")) {
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

# Calculate percent weight for each cell type
ce <- ddply(table_sample_type, "Celltype", transform, percent_weight = Cellnumber / sum(Cellnumber) * 100)
ce$Cancertype <- factor(ce$Cancertype, levels = c("Up", "Not", "Down"))
ce <- ddply(ce, "Celltype", transform, labely = cumsum(percent_weight))
p <- ggplot(ce, aes(x = Celltype, y = percent_weight, fill = Cancertype)) + geom_bar(stat = "identity") + theme_bw()
p_bar2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#FF9289", "lightgray", "#00DAE0")) +
  scale_y_continuous(name = "Gene Percent", expand = c(0, 0)) + labs(x = "Gene Sets", colour = "") +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(aes(y = labely, label = Cellnumber), vjust = 1.5, colour = "black", size = 6)
ggsave("Human_Geneset_percent313.pdf", p_bar2, width = 4, height = 6)

# Load necessary library
library(eulerr)

# Venn diagram for WT and MD proteins
WT = rownames(human_all2_org[rowSums(human_all2_org[, 3:4]) > 0, ])
MD = rownames(human_all2_org[rowSums(human_all2_org[, 5:6]) > 0, ])
vd <- euler(c(Mito = length(setdiff(MD, WT)), WT = length(setdiff(WT, MD)), "Mito&WT" = length(MD[MD %in% WT])))
pdf(file = "Human_protein_venn.pdf", width = 5, height = 5)
plot(vd, fills = list(fill = c("#00DAE0", "#FF9289", "#7fb6b4"), alpha = 0.6), labels = list(col = "black", font = 10), edges = FALSE, quantities = TRUE)
dev.off()

# Venn diagram for DEPs
UP = rownames(human_all3[human_all3$threshold == "Up", ])
DOWN = rownames(human_all3[human_all3$threshold == "Down", ])
NOT = rownames(human_all3[human_all3$threshold == "Not", ])
vd2 <- euler(c(UP = length(UP), Down = length(DOWN), "UP&Down" = length(NOT)))
pdf(file = "Human_protein_DEP_venn.pdf", width = 5, height = 5)
plot(vd2, fills = list(fill = c("#FF9289", "#00DAE0", "lightgray"), alpha = 0.6), labels = list(col = "black", font = 15), edges = FALSE, quantities = TRUE)
dev.off()

# Function annotation: GO and KEGG
# GO enrichment for upregulated proteins
dep_up <- human_deg2$Gene.Symbol[human_deg2$states == "Up"]
EG2symbol = toTable(org.Hs.egSYMBOL)
fid = as.character(dep_up)
geneLists = data.frame(symbol = fid)
results_up = merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_up = na.omit(results_up$gene_id)
ego_up <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_up, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_up_result <- ego_up@result[ego_up@result$Count >= 5 & ego_up@result$p.adjust < 0.05, ]
write.csv(ego_up_result, "Human_mito_depletion_dep_UPGO.csv")

# Draw GO enrichment bar plot for upregulated proteins
ego_up_draw <- ego_up_result[1:min(10,nrow(ego_up_result)), ]
p <- ggplot(ego_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (ceiling(max(ego_up_draw$Count) / 10) * 10 + 2))) + coord_flip() +
  scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria Depletion Up proteins GO Annotation") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))
ggsave(p, filename = "Human_mito_depletion_dep_UPGO.pdf", height = 5, width = 10)

# KEGG pathway enrichment for upregulated proteins
kegg_up <- enrichKEGG(organism = "hsa", gene = id_up, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_up <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_up_result <- kegg_up@result[kegg_up@result$Count >= 5 & kegg_up@result$p.adjust < 0.05, ]
write.csv(kegg_up_result, "Human_mito_depletion_dep_UPKEGG.csv")

# Draw KEGG pathway bar plot for upregulated proteins
kegg_up_draw <- kegg_up_result[1:min(10,nrow(kegg_up_result)), ]
p <- ggplot(kegg_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (ceiling(max(kegg_up_draw$Count) / 10) * 10 + 2))) + coord_flip() +
  scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria Depletion Up proteins KEGG Annotation") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))
ggsave(p, filename = "Human_mito_depletion_dep_UPKEGG.pdf", height = 5, width = 10)

# GO enrichment for downregulated proteins
dep_down <- human_deg2$Gene.Symbol[human_deg2$states == "Down"]
dep_down <- na.omit(dep_down)
fid = as.character(dep_down)
geneLists = data.frame(symbol = fid)
results_down = merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_down = na.omit(results_down$gene_id)
ego_down <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_down, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_down_result <- ego_down@result[ego_down@result$Count >= 5 & ego_down@result$p.adjust < 0.05, ]
write.csv(ego_down_result, "Human_mito_depletion_dep_downGO.csv")

# Draw GO enrichment bar plot for downregulated proteins
ego_down_draw <- ego_down_result[1:min(10,nrow(ego_down_result)), ]
p <- ggplot(ego_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (ceiling(max(ego_down_draw$Count) / 10) * 10 + 2))) + coord_flip() +
  scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria Depletion Down proteins GO Annotation") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))
ggsave(p, filename = "Human_mito_depletion_dep_downGO.pdf", height = 5, width = 10)

# KEGG pathway enrichment for downregulated proteins
kegg_down <- enrichKEGG(organism = "hsa", gene = id_down, keyType = "kegg", pvalueCutoff = 0.05, minGSSize = 10)
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down_result <- kegg_down@result[kegg_down@result$Count >= 5 & kegg_down@result$p.adjust < 0.05, ]
write.csv(kegg_down_result, "Human_mito_depletion_dep_downKEGG.csv")

# Draw KEGG pathway bar plot for downregulated proteins
kegg_down_draw <- kegg_down_result[1:min(10,nrow(kegg_down_result)), ]
p <- ggplot(kegg_down_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (ceiling(max(kegg_down_draw$Count) / 10) * 10 + 2))) + coord_flip() +
  scale_fill_gradient(low = "#00DAE0", high = "#00DAE0", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria Depletion Down proteins KEGG Annotation") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))
ggsave(p, filename = "Human_mito_depletion_dep_downKEGG.pdf", height = 5, width = 10)
