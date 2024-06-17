# Remove all existing objects from the workspace
rm(list=ls(all=T))
# Libraries loading
library(DESeq2)
library(ggplot2)
library(ggbiplot)
library(pheatmap)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(plyr)

# Read and sort TPM values for quality control
# Human data
oldh <- as.matrix(read.csv("human_all_result/gene_tpm_matrix.csv", row.names="gene_id"))
oldh <- oldh[order(rownames(oldh)),]
# Chimp data
oldc <- as.matrix(read.csv("chimp_all_result/gene_tpm_matrix.csv", row.names="gene_id"))
oldc <- oldc[order(rownames(oldc)),]

datahnew <- oldh[, c(1:5, 9)]
datacnew <- oldc[, c(1:5, 9)]

# Identify mitochondrial genes in human
hmt <- read.csv("Human_mt_gene_ids.txt", header = FALSE)
mtlist1 <- rownames(datahnew)[sapply(rownames(datahnew), function(x) strsplit(x, "\\|")[[1]][1]) %in% hmt$V1]
mtlist2 <- rownames(datahnew[grepl("\\|MT-", rownames(datahnew)),])
if (identical(mtlist1, mtlist2)) {
  datahnew_mt    <- colSums(datahnew[grepl("\\|MT-", rownames(datahnew)),])
  datahnew_other <- colSums(datahnew[!grepl("\\|MT-", rownames(datahnew)),])
}

# Identify mitochondrial genes in chimp
cmt <- read.table("MT_gene_ids.txt", header = FALSE)
cname <- sapply(rownames(datacnew), function(x) strsplit(x, "\\|")[[1]][1])
mtindex <- cname %in% cmt$V1
datacnew_mt    <- colSums(datacnew[mtindex,])
datacnew_other <- colSums(datacnew[!mtindex,])

# Create plots for QC analysis
# Nuclear genes
otherdraw <- rbind(cbind(datahnew_other, rownames(datahnew_other), c(rep("HCh", 3), rep("HCc", 3)), rep("Human", 6)),
                   cbind(datacnew_other, rownames(datacnew_other), c(rep("HCh", 3), rep("HCc", 3)), rep("Chimp", 6)))
otherdraw <- as.data.frame(otherdraw)
colnames(otherdraw) <- c("count", "sample", "species")
otherdraw$count <- as.numeric(otherdraw$count)
library(dplyr)
summary_other <- otherdraw %>%
  group_by(sample, species) %>%
  summarise(
    mean_count = mean(count),
    se = sd(count) / sqrt(3),   # Calculate standard error
    .groups = 'drop'
  )
summary_other$sample <- factor(summary_other$sample, levels=c("HCh", "HCc"))
summary_other$species <- factor(summary_other$species, levels=c("Human", "Chimp"))
p <- ggplot(summary_other, aes(x=sample, y=mean_count, fill=species)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 2))) +
  geom_bar(stat = "identity", width = 0.6, position=position_dodge(0.7)) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), width = 0.2, position = position_dodge(0.7)) +
  theme_bw() + scale_fill_manual(values = c("Human"="#1E88E5", "Chimp"="red"))
pdf(file = "HC_nuclear_bar_tpm.pdf", width = 4, height = 3)
print(p)
dev.off()

# Mitochondrial genes
mtdraw <- rbind(cbind(datahnew_mt, rownames(datahnew_mt), c(rep("HCh", 3), rep("HCc", 3)), rep("Human", 6)),
                cbind(datacnew_mt, rownames(datacnew_mt), c(rep("HCh", 3), rep("HCc", 3)), rep("Chimp", 6)))
mtdraw <- as.data.frame(mtdraw)
colnames(mtdraw) <- c("count", "sample", "species")
mtdraw$count <- as.numeric(mtdraw$count)
library(dplyr)
summary_mt <- mtdraw %>%
  group_by(sample, species) %>%
  summarise(
    mean_count = mean(count),
    se = sd(count) / sqrt(3),   # Calculate standard error
    .groups = 'drop'
  )
summary_mt$sample <- factor(summary_mt$sample, levels=c("HCh", "HCc"))
summary_mt$species <- factor(summary_mt$species, levels=c("Human", "Chimp"))
p <- ggplot(summary_mt, aes(x=sample, y=mean_count, fill=species)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 2))) +
  geom_bar(stat = "identity", width = 0.6, position=position_dodge(0.7)) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), width = 0.2, position = position_dodge(0.7)) +
  theme_bw() + scale_fill_manual(values = c("Human"="#1E88E5", "Chimp"="red"))
pdf(file = "HC_mt_bar_tpm.pdf", width = 4, height = 4)
print(p)
dev.off()


# ---------- RNA-seq data integrity check ------------------------------
controlnum = 3   # Number of control samples
casenum = 3      # Number of case samples
num = casenum + controlnum  # Total number of samples

# Load and order historical RNA-seq data for humans
oldh <- as.matrix(read.csv("H-C/result/human_all_result/gene_count_matrix.csv", row.names="gene_id"))
oldh <- oldh[order(rownames(oldh)),]

# Load and order historical RNA-seq data for chimpanzees
oldc <- as.matrix(read.csv("H-C/result/chimp_all_result/gene_count_matrix.csv", row.names="gene_id"))
oldc <- oldc[order(rownames(oldc)),]

datahnew <- oldh[,c(1:5,9)]
datacnew <- oldc[,c(1:5,9)]

# Read mitochondrial gene IDs for humans and check for matches
hmt <- read.csv("Human_mt_gene_ids.txt", header = FALSE)
mtlist1 <- rownames(datahnew)[sapply(rownames(datahnew), function(x) strsplit(x, "\\|")[[1]][1]) %in% hmt$V1]
mtlist2 <- rownames(datahnew[grepl("\\|MT-", rownames(datahnew)),])
if (identical(mtlist1, mtlist2)) {
  datahnew_mt    <- colSums(datahnew[grepl("\\|MT-", rownames(datahnew)),])
  datahnew_other <- colSums(datahnew[!grepl("\\|MT-", rownames(datahnew)),])
}

# Read mitochondrial gene IDs for chimpanzees and check for matches
cmt <- read.table("MT_gene_ids.txt", header = FALSE)
cname <- sapply(rownames(datacnew), function(x) strsplit(x, "\\|")[[1]][1])
mtindex <- cname %in% cmt$V1
datacnew_mt    <- colSums(datacnew[mtindex,])
datacnew_other <- colSums(datacnew[!mtindex,])

# Combine human and chimpanzee data into a single dataframe
dataall <- as.data.frame(rbind(datahnew, datacnew))
dataall <- cbind(dataall, sapply(rownames(dataall), function(x) strsplit(x, "\\|")[[1]][2]))
colnames(dataall) <- c("HC_h1", "HC_h2", "HC_h3", "HC_c1", "HC_c2", "HC_c3", "genename")
dataall <- dataall[!is.na(dataall$genename),]

# Identify and process duplicate gene names
dname <- unique(dataall$genename[duplicated(dataall$genename)])
for (i in 1:length(dname)) {
  iindex = which(dataall$genename == dname[i])
  mmat = colMeans(apply(dataall[iindex, 1:6], 2, as.numeric))
  if (i == 1) {
    dataall2 <- rbind(dataall, c(mmat, dname[i]))
    rmindex <- iindex
  } else {
    dataall2 <- rbind(dataall2, c(mmat, dname[i]))
    rmindex <- c(rmindex, iindex)
  }
  rownames(dataall2)[nrow(dataall2)] <- paste(rownames(dataall[iindex,]), collapse="|")
}
dataall2 <- dataall2[-rmindex,]  # Remove original rows for duplicates after merging

# all sample pca

# Load and preprocess human control gene count data
hcontrol <- as.matrix(read.csv("/home/DX6/jwulab/S227384/data/Denial/orangutan-human/run/hcontrol/hcontrolall_genecount/gene_count_matrix.csv", row.names="gene_id"))
rownames(hcontrol) <- sapply(rownames(hcontrol), function(x) strsplit(x, "\\|")[[1]][2])
hcontrol <- hcontrol[!is.na(rownames(hcontrol)),]
hcontrol2 <- process_duplicates(hcontrol)  # Assuming process_duplicates is a predefined function

# Load and preprocess chimpanzee control gene count data
ccontrol <- as.matrix(read.csv("/home/DX6/jwulab/S227384/data/Denial/orangutan-human/run/ccontrol/ccontrolall_genecount/gene_count_matrix.csv", row.names="gene_id"))
rownames(ccontrol) <- sapply(rownames(ccontrol), function(x) strsplit(x, "\\|")[[1]][2])
ccontrol <- ccontrol[!is.na(rownames(ccontrol)),]
ccontrol2 <- process_duplicates(ccontrol)

# Load and preprocess orangutan control gene count data
ocontrol <- as.matrix(read.csv("/home/DX6/jwulab/S227384/data/Denial/orangutan-human/run/ocontrol/ocontrolall_genecount/gene_count_matrix.csv", row.names="gene_id"))
rownames(ocontrol) <- sapply(rownames(ocontrol), function(x) strsplit(x, "\\|")[[1]][2])
ocontrol <- ocontrol[!is.na(rownames(ocontrol)),]
ocontrol2 <- process_duplicates(ocontrol)

# Load and preprocess HO all data for PCA
load("HO_all/hoall.rdata")
# Process row names to extract the second part of the name, assumed to be the relevant gene identifier
rownames(hoall) <- sapply(rownames(hoall), function(x) strsplit(x, "\\|")[[1]][2])
# Filter out rows with NA in row names to ensure data integrity
hoall <- hoall[!is.na(rownames(hoall)),]
# Order rows by gene names to ensure consistency in analysis
hoall <- hoall[order(rownames(hoall)),]
# Convert the first 6 columns to numeric for PCA analysis and assign row names
hoall2 <- apply(as.matrix(hoall[,1:6]), 2, as.numeric)
rownames(hoall2) <- hoall$genename
# Prepare a similar matrix for control data
hcall <- apply(as.matrix(dataall2[,1:6]), 2, as.numeric)
rownames(hcall) <- dataall2$genename
hcall <- hcall[!is.na(rownames(hcall)),]
hcall <- hcall[order(rownames(hcall)),]
# Align human data to the available control data
hccom <- hcall[rownames(hcall) %in% rownames(hcontrol2),]
hcom <- hcontrol2[rownames(hcontrol2) %in% rownames(hcall),]
if (identical(rownames(hccom), rownames(hcom))) {
  pcamat <- cbind(hccom, hcom)
} else {
  print("Not same!")  # Error message if row names do not match
}
# Align chimpanzee control data to the combined human data
hccom2 <- pcamat[rownames(pcamat) %in% rownames(ccontrol2),]
ccom <- ccontrol2[rownames(ccontrol2) %in% rownames(pcamat),]
if (identical(rownames(hccom2), rownames(ccom))) {
  pcamat <- cbind(hccom2, ccom)
} else {
  print("Not same!")
}
# Align orangutan control data to the existing PCA matrix
hccom3 <- pcamat[rownames(pcamat) %in% rownames(ocontrol2),]
ocom <- ocontrol2[rownames(ocontrol2) %in% rownames(pcamat),]
if (identical(rownames(hccom3), rownames(ocom))) {
  pcamat <- cbind(hccom3, ocom)
} else {
  print("Not same!")
}
# Integrate the last set of data for PCA
hccom4 <- pcamat[rownames(pcamat) %in% rownames(hoall2),]
hocom <- hoall2[rownames(hoall2) %in% rownames(pcamat),]
if (identical(rownames(hccom4), rownames(hocom))) {
  pcamat <- cbind(hccom4, hocom)
} else {
  print("Not same!")
}
# Transpose the matrix for PCA analysis, filter columns with zero variance
pcamat <- as.data.frame(t(pcamat))
pcamat_filtered <- pcamat[, apply(pcamat, 2, function(x) var(x, na.rm = TRUE) != 0)]
# Define groups for plotting in PCA
pcamatgroup <- c(rep("HC_H", 3), rep("HC_C", 3), rep("H_control", 3), rep("C_control", 3),
                 rep("O_control", 3), rep("HO_H", 3), rep("HO_O", 3))
# Perform PCA
pca_result <- prcomp(pcamat_filtered, scale = TRUE, center = FALSE)
# Plot PCA
pdf(file = "Allsample_pca2.pdf", width = 4, height = 4)
ggbiplot(pca_result,
         var.axes = FALSE,            # Do not show arrows for variables
         obs.scale = 1,               # Scale for observations
         groups = factor(group_labels, levels = c("C_control", "HC_C", "HC_H", "H_control", "HO_H", "HO_O", "O_control")),
         ellipse = TRUE,
         circle = FALSE) +
  geom_text(aes(label = rownames(pcamat_filtered)), vjust = 1.5, size = 2) +
  theme_bw() +
  scale_color_manual(values = c("HC_C" = "#F8766C", "HC_H" = "#1E88E5", "H_control" = "blue", "C_control" = "red",
                                "O_control" = "#F16821", "HO_H" = "#008080", "HO_O" = "orange"))
dev.off()
# ---------- Compare treatment (HC_H) vs control (HC_C) groups ----------
sampletitle = "HC_H vs HC_C"
# Prepare the count data matrix and round to integer
countData <- apply(as.matrix(dataall2[, 1:6]), 2, as.numeric)
rownames(countData) <- rownames(dataall2)
countData <- round(countData, digits = 0)
# Define conditions for differential expression analysis
condition <- factor(c(rep("treat", casenum), rep("untreat", controlnum)))
colData <- data.frame(row.names = colnames(countData), condition)
# Create a DESeq2 dataset and run the analysis
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize = TRUE)), by = "row.names", sort = FALSE)
resdata <- resdata[!is.na(resdata$padj),]
resdata <- resdata[order(resdata$padj),]
# Gene annotation and adjustment
# Add column for gene names extracted from row names
resdata$genename = sapply(resdata$Row.names, function(x) strsplit(x, "\\|")[[1]][2])
# Add gene type information for mitochondrial coding and associated genes
MT_coding <- read.table("/home/DX6/jwulab/S227384/PosdocProject/Denial/genedatabase/Human mitochondria genome coding genes.csv", sep = ",", header = FALSE)
MT_asso <- read.table("/home/DX6/jwulab/S227384/PosdocProject/Denial/genedatabase/human_chromosome_mito_genelist.csv", sep = ",", header = FALSE)
resdata$genetype = NA
resdata$genetype[resdata$genename %in% MT_coding[, 2]] <- "MT_coding"
resdata$genetype[resdata$genename %in% MT_asso[, 2]] <- "MT_associated"
# Adjust the log2 fold changes for direction
resdata$log2FoldChange = -resdata$log2FoldChange
# Identify significant DEGs
diff_gene_p <- resdata[resdata$pvalue < 0.05,]
diff_gene_p$stat <- ifelse(diff_gene_p$log2FoldChange > 0, "Up", "Down")
diff_gene_p <- diff_gene_p[abs(diff_gene_p$log2FoldChange) > 1.5,]
write.table(diff_gene_p, sep = "\t", row.names = FALSE, quote = FALSE, "HCh_vs_HCc_all_DEGs.txt")
write.table(resdata, sep = "\t", row.names = FALSE, quote = FALSE, "HCh_vs_HC_allgenes.txt")

# Load the dataset into 'data' and create a binary 'threshold' variable to indicate significant genes
data <- resdata
data$threshold <- as.factor(ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1.5,
                                   ifelse(data$log2FoldChange >= 1.5, 'Up', 'Down'), 'Not'))
# Clear 'genetype' and 'genename' for rows not in 'diff_gene_p'
data$genetype[!data$Row.names %in% diff_gene_p$Row.names] <- NA
data$genename[is.na(data$genetype)] <- NA

# Calculate average expression values for cases and controls
data$pos <- rowMeans(data[,8:(8 + casenum - 1)])
data$neg <- rowMeans(data[,(8 + casenum):(8 + num - 1)])

# Scatter plot visualizing positive vs negative expression values color-coded by gene significance
p <- ggplot(data, aes(x=log2(pos), y=log2(neg), colour=threshold)) + geom_point()
p1 <- p + geom_point(aes(color=threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha=0.4, size=1)
p3 <- p2 + theme_bw() + labs(x="Case_genes", y="Control_genes", title=sampletitle) +
  theme(plot.title = element_text(hjust=0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept=0, slope=1, lty=4, col="black", lwd=0.6)
a <- cor.test(data$pos, data$neg, method="pearson")
p6 <- p5 + annotate("text", x=14, y=-0.5, label=paste("Pearson's correlation ", round(a$estimate, 3)), size=3)
ggsave(p6, file="HC_All_correlation.pdf", width=8, height=8)

# Volcano plot for gene expression data, adjusting color based on gene significance
vdata <- arrange(data, padj)
vdata$threshold <- as.factor(ifelse(vdata$pvalue < 0.05 & abs(vdata$log2FoldChange) >= 1.5,
                                    ifelse(vdata$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
p <- ggplot(vdata, aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold))
p1 <- p + geom_point(aes(color=vdata$threshold)) +
  geom_vline(xintercept=c(-1.5, 1.5), lty=4, col="black", lwd=0.5) +
  geom_hline(yintercept=-log10(0.05), lty=4, col="black", lwd=0.5)
p2 <- p1 + scale_color_manual(values=c("blue", "grey", "red")) + geom_point(alpha=0.4, size=1.2)
p3 <- p2 + theme_bw() + labs(x="log2(FoldChange)", y="-log10(p-value)", title=sampletitle) +
  theme(plot.title = element_text(hjust=0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(-2, 210) + xlim(-10, 10)
ggsave(v4, file=paste(sampletitle, "_volcano.pdf", sep=""), width=4, height=3.5)

# Prepare the data matrix from the differential gene expression results
datamat <- as.matrix(diff_gene_p[, c(8:13)])  # Assuming columns 8 to 13 contain expression data
rownames(datamat) <- diff_gene_p$genename    # Set row names to gene names
datamat = datamat[!is.na(diff_gene_p$genename),]  # Remove rows where gene names are NA

# Prepare an annotation dataframe for gene names and types
anno_gene2 <- cbind(diff_gene_p$genename, diff_gene_p$genetype)  # Combine gene names and types
colnames(anno_gene2) <- c("genename", "genetype")  # Name the columns
anno_gene2 <- as.data.frame(anno_gene2)  # Convert to dataframe
anno_gene2 <- anno_gene2[!is.na(anno_gene2$genetype),]  # Remove rows with NA gene types

# Load the ComplexHeatmap library
library(ComplexHeatmap)

# Convert expression data matrix to matrix form and remove mitochondrial genes (prefixed with "MT-")
A <- as.matrix(datamat)
A <- A[!grepl("MT-", rownames(A)),]

# Define sample groups for heatmap annotation
samples <- rep(c('HC_H', 'HC_C'), c(3, 3))  # Assuming 3 samples each for HC_H and HC_C groups

# Standardize the data row-wise after log transformation
for (i in 1:nrow(A)) {
  A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))
}

# Create a heatmap object with specified color palette and annotations
B <- Heatmap(
  A,  # Matrix of expression data
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Define color gradient
  show_row_names = TRUE,  # Show row names (gene names)
  row_names_gp = gpar(fontsize = 6),  # Font size for row names
  top_annotation = HeatmapAnnotation(  # Define top annotation for sample groups
    Group = samples, 
    simple_anno_size = unit(2, 'mm'), 
    col = list(Group = c('HC_H' = '#F8766D', 'HC_C' = '#619CFF')),
    show_annotation_name = FALSE
  )
)

# Output the heatmap to a PDF file
pdf(file = "HC_all_heatmap_rmMT.pdf", width = 6, height = 10)
print(B)  # Print the heatmap
dev.off()  # Close the PDF device

