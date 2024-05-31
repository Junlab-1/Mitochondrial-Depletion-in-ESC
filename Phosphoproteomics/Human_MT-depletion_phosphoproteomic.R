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
setwd("~/Human proteomic/phosphoprotein")

# Load and clean data
human_phosphoprotein <- read.csv("Phosphoprotein1.csv", header = TRUE)
human_phos <- na.omit(human_phosphoprotein)
rmlist = 0
for (i in 1:nrow(human_phos)) {
  if (sum(is.na(human_phos[i, 3:8])) > 3) {
    if (sum(rmlist) == 0) {
      rmlist = i
    } else {
      rmlist = c(rmlist, i)
    }
  }
}
if (rmlist == 0) {
  human_phos2_org <- human_phos
} else {
  human_phos2_org <- human_phos[-rmlist, ]
}
# Normalize data
human_phos2[, 3:6] <- log2(human_phos2_org[, 3:6])
human_phos2[, 3:6][human_phos2[, 3:6] == -Inf] <- 0
for (i in 3:6) {
  human_phos2[, i] <- human_phos2[, i] - mean(human_phos2[, i])
}
human_phos2$avg <- rowMeans(human_phos2[, 3:6])
human_phos3 = human_phos2
for (i in 3:6) {
  imodel <- lm(human_phos2$avg ~ human_phos2[, i])
  islope <- coef(imodel)[2]
  human_phos3[, i] <- human_phos2[, i] / islope
}

# Call DEGs
for (i in 1:nrow(human_phos3)) {
  iresult <- t.test(human_phos3[i, 3:4], human_phos3[i, 5:6])
  human_phos3$pvalue[i] <- iresult$p.value
  if (human_phos3[i, 1] == human_phos2_org[i, 1]) {
    human_phos3$log2fc[i] <- log2(mean(as.numeric(human_phos2_org[i, 5:6])) / mean(as.numeric(human_phos2_org[i, 3:4])))
  }
}

# Annotate MT genes and pluripotency genes
PluripotencyGene <- read.table("Pluripotency_Genes.txt", header = FALSE)
MT_coding <- read.table("Human mitochondria genome coding genes.csv", sep = ",", header = TRUE)
MT_asso <- read.table("human_chromosome_mito_genelist.csv", sep = ",", header = TRUE)
human_phosdeg <- human_phos3[abs(human_phos3$log2fc) > 1 & human_phos3$pvalue < 0.05, ]
human_phosdeg$proteintype <- NA
human_phosdeg$proteintype[human_phosdeg$Gene.Symbol %in% MT_coding[, 2]] <- "MT_coding"
human_phosdeg$proteintype[human_phosdeg$Gene.Symbol %in% MT_asso[, 2]] <- "MT_associated"
human_phosdeg$proteintype[human_phosdeg$Gene.Symbol %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"
write.csv(human_phosdeg, "human_phos_proteins.csv")
human_phosdeg$states <- "Up"
human_phosdeg$states[human_phosdeg$log2fc < 0] <- "Down"
a = human_phos2_org$Accession[human_phos2_org$Accession %in% human_phosdeg$Accession]
if (identical(human_phosdeg$Accession, a)) {
  human_phosdeg[, 3:6] <- human_phos2_org[human_phos2_org$Accession %in% human_phosdeg$Accession, 3:6]
}
human_phosdeg$avg <- NULL
write.csv(human_phosdeg, "human_phos_proteins2.csv")

# Volcano plot
sampletitle <- "Human Mito vs H9 phosphoproteins"
data <- human_phos3
if (identical(human_phos2_org$Accession, data$Accession)) {
  data[, 3:6] <- human_phos2_org[, 3:6]
}
data$proteintype <- NA
data$proteintype[data$Gene.Symbol %in% MT_coding[, 2]] <- "MT_coding"
data$proteintype[data$Gene.Symbol %in% MT_asso[, 2]] <- "MT_associated"
data$proteintype[data$Gene.Symbol %in% PluripotencyGene[, 1]] <- "PluripotencyGenes"
write.csv(data, "human_phos_proteins_all.csv")
data$threshold <- as.factor(ifelse(data$pvalue < 0.05 & abs(data$log2fc) >= 1, ifelse(data$log2fc > 0, 'Up', 'Down'), 'Not'))
data$neg <- rowMeans(data[, 3:4])
data$pos <- rowMeans(data[, 5:6])
p <- ggplot(data, aes(x = log2(pos), y = log2(neg), colour = threshold)) + geom_point()
p1 <- p + geom_point(aes(color = threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha = 0.4, size = 1)
p3 <- p2 + theme_bw() + labs(x = "Case_genes", y = "Control_genes", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 0.6)
a <- cor.test(data$pos, data$neg, method = "pearson")
p6 <- p5 + annotate("text", x = 20, y = 10, label = paste("Pearson's correlation ", round(a$estimate, 3)), size = 3)
ggsave(p6, file = "Human_phosphoprotein_correlation2.pdf", width = 9, height = 8)

# Volcano plot for DEGs
vdata <- arrange(data, pvalue)
vdata$threshold <- as.factor(ifelse(vdata$pvalue < 0.05 & abs(vdata$log2fc) >= 1, ifelse(vdata$log2fc > 0, 'Up', 'Down'), 'Not'))
vdata$proteintype[vdata$threshold == "Not"] <- "Not"
vdata$proteintype <- factor(vdata$proteintype, levels = c("MT_associated", "Not", "MT_coding", "PluripotencyGenes"))
p <- ggplot(vdata, aes(x = log2fc, y = -log10(pvalue)))
p1 <- p + geom_point(aes(shape = threshold, color = proteintype), size = 1) +
  geom_point(data = subset(vdata, proteintype %in% c("MT_coding", "MT_associated", "PluripotencyGenes")),
             aes(shape = threshold, color = proteintype), size = 1) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.5) + scale_x_continuous(limits = c(-6.2, 6.2), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_shape_manual(values = c(2, 16, 4)) + scale_color_manual(values = c("#f8766d", "lightgray", "#0033A0"))
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(v4, file = "Human_phosphoprotein_volcano.pdf", width = 4, height = 5)

# PCA plot
pcamat <- as.data.frame(t(human_phos2_org[, c(3, 4, 6, 7)]))
rownames(pcamat) <- c("WT_Day0_S1", "WT_Day0_S2", "MT_Depletion_S1", "MT_Depletion_S2")
pcamatgroup <- c(rep("WT", 2), rep("MT_Depletion", 2))
pca_result <- prcomp(pcamat, scale = FALSE, center = FALSE)
pca_result$x <- scale(pca_result$x)

# Plot PCA results
pdf(file = "Human_phosphoprotein_MT_Depletion_pcaplot2.pdf", width = 5, height = 5)
ggbiplot(pca_result,
         var.axes = FALSE,
         obs.scale = 1,
         groups = pcamatgroup,
         ellipse = FALSE,
         circle = FALSE) +
  geom_text(aes(label = rownames(pcamat)), vjust = 1.5, size = 2) +
  theme_bw()
dev.off()

# Heatmap
datamat <- as.matrix(human_phosdeg[, 3:6])
rownames(datamat) <- human_phosdeg$Gene.Symbol
human_phosdeg$proteintype[is.na(human_phosdeg$proteintype)] <- "None"
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
C <- B + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% human_phosdeg$Gene.Symbol[human_phosdeg$proteintype == "MT_associated"]),
                                        labels = human_phosdeg$Gene.Symbol[human_phosdeg$proteintype == "MT_associated"], labels_gp = gpar(fontsize = 4, col = "#f8766d")))
D <- C + rowAnnotation(link = anno_mark(at = which(rownames(A) %in% human_phosdeg$Gene.Symbol[human_phosdeg$proteintype == "PluripotencyGenes"]),
                                        labels = human_phosdeg$Gene.Symbol[human_phosdeg$proteintype == "PluripotencyGenes"], labels_gp = gpar(fontsize = 4, col = "#0033A0")))
pdf(file = "Human_phosphoDEP_heatmap.pdf", width = 8, height = 8)
print(D)
dev.off()

# Venn diagram for WT and MD phosphoproteins
library(eulerr)
WT = rownames(human_phos2_org[rowSums(human_phos2_org[, 3:4]) > 0, ])
MD = rownames(human_phos2_org[rowSums(human_phos2_org[, 5:6]) > 0, ])
vd <- euler(c(Mito = length(setdiff(MD, WT)), WT = length(setdiff(WT, MD)), "Mito&WT" = length(MD[MD %in% WT])))
pdf(file = "Human_phosphoprotein_venn.pdf", width = 5, height = 5)
plot(vd, fills = list(fill = c("#00DAE0", "#FF9289", "#7fb6b4"), alpha = 0.6), labels = list(col = "black", font = 10), edges = FALSE, quantities = TRUE)
dev.off()

# Venn diagram for DEGs in phosphoproteins
human_phos3$threshold <- as.factor(ifelse(human_phos3$pvalue < 0.05 & abs(human_phos3$log2fc) >= 1, ifelse(human_phos3$log2fc > 0, 'Up', 'Down'), 'Not'))
UP = rownames(human_phos3[human_phos3$threshold == "Up", ])
DOWN = rownames(human_phos3[human_phos3$threshold == "Down", ])
NOT = rownames(human_phos3[human_phos3$threshold == "Not", ])
vd2 <- euler(c(UP = length(UP), Down = length(DOWN), "UP&Down" = length(NOT)))
pdf(file = "Human_protein_DEP_venn.pdf", width = 5, height = 5)
plot(vd2, fills = list(fill = c("#FF9289", "#00DAE0", "lightgray"), alpha = 0.6), labels = list(col = "black", font = 15), edges = FALSE, quantities = TRUE)
dev.off()

# Function annotation
dep_up <- human_phosdeg$Gene.Symbol[human_phosdeg$states == "Up"]
EG2symbol = toTable(org.Hs.egSYMBOL)
fid = as.character(dep_up)
geneLists = data.frame(symbol = fid)
results_up = merge(geneLists, EG2symbol, by = 'symbol', all.x = TRUE)
id_up = na.omit(results_up$gene_id)
ego_up <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id_up, ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)
ego_up_result <- ego_up@result[ego_up@result$Count >= 3 & ego_up@result$pvalue < 0.05, ]
write.csv(ego_up_result, "Human_mito_depletion_phosphodep_UPGO.csv")

# Draw GO enrichment bar plot for upregulated phosphoproteins
ego_up_draw <- ego_up_result[1:min(10, nrow(ego_up_result)), ]
p <- ggplot(ego_up_draw, aes(x = reorder(Description, Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + scale_y_continuous(name = "Gene Count", expand = c(0, 0), limits = c(0, (ceiling(max(ego_up_draw$Count) / 2) * 2 + 0.2))) + coord_flip() +
  scale_fill_gradient(low = "#FF9289", high = "#FF9289", name = "Genes") +
  theme(text = element_text(size = 16, family = "serif")) + labs(title = "Human Mitochondria Depletion Up proteins GO Annotation") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))
ggsave(p, filename = "Human_mito_depletion_phosphodep_UPGO.pdf", height = 5, width = 10)


# Order data by p-value
human_phosdeg2 <- human_phosdeg[order(human_phosdeg$pvalue), ]

# Select top 10 upregulated and downregulated genes
phos_draw <- human_phosdeg2[human_phosdeg2$states == "Up", ][1:10, ]
phos_down_draw <- human_phosdeg2[human_phosdeg2$states == "Down", ][1:10, ]

# Convert p-values to -log10(p-value)
phos_draw$pvalue <- -log10(phos_draw$pvalue)
phos_down_draw$pvalue <- -log10(phos_down_draw$pvalue)

# Plot upregulated genes top10
p <- ggplot(phos_draw, aes(x = reorder(Gene.Symbol, pvalue), y = pvalue, fill = pvalue)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + coord_flip() + scale_y_continuous(name = "-log10 Pvalue", expand = c(0, 0), breaks = seq(0, 3, by = 0.5)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
p2 <- p + scale_fill_viridis_c()
ggsave(p2, filename = "Human_phospho_global_up.pdf", height = 5, width = 6)

# Plot downregulated genes top10
p <- ggplot(phos_down_draw, aes(x = reorder(Gene.Symbol, pvalue), y = pvalue, fill = pvalue)) +
  geom_bar(stat = "identity", colour = "NA", width = 0.9, position = position_dodge(0.7)) +
  xlab("") + coord_flip() + scale_y_continuous(name = "-log10 Pvalue", expand = c(0, 0), breaks = seq(0, 3, by = 0.5)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
p2 <- p + scale_fill_viridis_c()
ggsave(p2, filename = "Human_phospho_global_down.pdf", height = 5, width = 6)
