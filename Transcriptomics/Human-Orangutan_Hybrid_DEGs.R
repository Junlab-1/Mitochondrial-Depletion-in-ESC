# Clear all existing objects in the workspace
rm(list=ls(all=T))

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(plyr)
library(dplyr)  # Load dplyr for data manipulation

# Set working directory for input/output
setwd("Rresult/HO_all")

# Load and sort human TPM data
datah <- as.matrix(read.csv("ho/all_genecount/human_ho_tpm.csv", row.names="gene_id"))
datah <- datah[order(rownames(datah)),]

# Load and sort orangutan TPM data
datao <- as.matrix(read.csv("ho/all_genecount/orangutan_ho_tpm.csv", row.names="gene_id"))
datao <- datao[order(rownames(datao)),]

# Load orangutan mitochondrial gene IDs
omt <- read.table("Orangutan_MT_gene_ids.txt", header = FALSE)
oname <- sapply(rownames(datao), function(x) strsplit(x, "\\|")[[1]][1])
mtindex <- oname %in% omt$V1

# Sum nuclear and mitochondrial gene expression separately
datahnew_mt <- colSums(datah[grepl("\\|MT-", rownames(datah)),])
datahnew_other <- colSums(datah[!grepl("\\|MT-", rownames(datah)),])
datacnew_mt <- colSums(datao[mtindex,])
datacnew_other <- colSums(datao[!mtindex,])

# Update row names for mitochondrial genes with more specific identifiers
for (i in 1:nrow(datao[mtindex,])) {
  mn <- rownames(datao)[mtindex]
  i_name <- strsplit(mn[i], "\\|")[[1]][2]
  iname <- if (is.na(i_name)) {
    paste(mn[i], "|MT-", i, sep = "")
  } else {
    paste(strsplit(mn[i], "\\|")[[1]][1], "|MT-", i_name, sep = "")
  }
  rownames(datao)[mtindex][i] <- iname
}

# Create combined data for plotting nuclear gene expression
otherdraw <- rbind(cbind(datahnew_other, rownames(datahnew_other), c(rep("HOh", 3), rep("HOo", 3)), rep("Human", 6)),
                   cbind(datacnew_other, rownames(datacnew_other), c(rep("HOh", 3), rep("HOo", 3)), rep("Orangutan", 6)))
otherdraw <- as.data.frame(otherdraw)
colnames(otherdraw) <- c("count", "sample", "specie")
otherdraw$count <- as.numeric(otherdraw$count)

# Summarize and plot nuclear gene expression data
summary_other <- otherdraw %>%
  group_by(sample, specie) %>%
  summarise(mean_count = mean(count), se = sd(count) / sqrt(3), .groups = 'drop')
summary_other$sample <- factor(summary_other$sample, levels=c("HOh", "HOo"))
summary_other$specie <- factor(summary_other$specie, levels=c("Human", "Orangutan"))
p <- ggplot(summary_other, aes(x=sample, y=mean_count, fill=specie)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 2))) +
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), width = 0.2, position = position_dodge(0.7)) +
  theme_bw() + scale_fill_manual(values = c("Human" = "#1E88E5", "Orangutan" = "#FF8000"))
pdf(file = "HO_nuclear_bar.pdf", width = 4, height = 3)
print(p)
dev.off()

# Create combined data for plotting mitochondrial gene expression
mtdraw <- rbind(cbind(datahnew_mt, rownames(datahnew_mt), c(rep("HOh", 3), rep("HOo", 3)), rep("Human", 6)),
                cbind(datacnew_mt, rownames(datacnew_mt), c(rep("HOh", 3), rep("HOo", 3)), rep("Orangutan", 6)))
mtdraw <- as.data.frame(mtdraw)
colnames(mtdraw) <- c("count", "sample", "specie")
mtdraw$count <- as.numeric(mtdraw$count)

# Summarize and plot mitochondrial gene expression data
summary_mt <- mtdraw %>%
  group_by(sample, specie) %>%
  summarise(mean_count = mean(count), se = sd(count) / sqrt(3), .groups = 'drop')
summary_mt$sample <- factor(summary_mt$sample, levels=c("HOh", "HOo"))
summary_mt$specie <- factor(summary_mt$specie, levels=c("Human", "Orangutan"))
p <- ggplot(summary_mt, aes(x=sample, y=mean_count, fill=specie)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 2))) +
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), width = 0.2, position = position_dodge(0.7)) +
  theme_bw() + scale_fill_manual(values = c("Human" = "#1E88E5", "Orangutan" = "#FF8000"))
pdf(file = "HO_mt_bar.pdf", width = 4, height=4)
print(p)
dev.off()

# ---------------------------- calculated DEGs (Differentially Expressed Genes) ----------------------------
# Set the number of control and case samples
controlnum <- 3
casenum <- 3
num <- casenum + controlnum

# Read and preprocess human gene count data
datah <- as.matrix(read.csv(
  "ho/all_genecount/human_ho_gene.csv",
  row.names="gene_id"))
datah <- datah[order(rownames(datah)),]

# Read and preprocess orangutan gene count data
datao <- as.matrix(read.csv(
  "ho/all_genecount/orangutan_ho_gene.csv",
  row.names="gene_id"))
datao <- datao[order(rownames(datao)),]

# Load orangutan mitochondrial gene IDs
omt <- read.table("MT_gene_ids.txt", header = FALSE)
oname <- sapply(rownames(datao), function(x) strsplit(x, "\\|")[[1]][1])
mtindex <- oname %in% omt$V1
datahnew_mt <- colSums(datah[grepl("\\|MT-", rownames(datah)),])
datahnew_other <- colSums(datah[!grepl("\\|MT-", rownames(datah)),])
datacnew_mt <- colSums(datao[mtindex,])
datacnew_other <- colSums(datao[!mtindex,])

# Update mitochondrial gene names in orangutan data
for (i in 1:nrow(datao[mtindex,])) {
  mn <- rownames(datao)[mtindex]
  i_name <- strsplit(mn[i], "\\|")[[1]][2]
  if (is.na(i_name)) {
    iname <- paste(mn[i], "|MT-", i, sep = "")
  } else {
    iname <- paste(strsplit(mn[i], "\\|")[[1]][1], "|MT-", i_name, sep = "")
  }
  rownames(datao)[mtindex][i] <- iname
}

# Combine human and orangutan data and extract gene names
dataall <- as.data.frame(rbind(datah, datao))
dataall <- cbind(dataall, sapply(rownames(dataall), function(x) strsplit(x, "\\|")[[1]][2]))
colnames(dataall) <- c("HO_h1", "HO_h2", "HO_h3", "HO_o1", "HO_o2", "HO_o3", "genename")
dataall <- dataall[!is.na(dataall$genename),]

# Identify and average duplicate gene entries
dname <- unique(dataall$genename[duplicated(dataall$genename)])
for (i in 1:length(dname)) {
  iindex <- which(dataall$genename == dname[i])
  mmat <- colMeans(apply(dataall[iindex, 1:6], 2, as.numeric))
  if (i == 1) {
    dataall2 <- rbind(dataall, c(mmat, dname[i]))
    rmindex <- iindex
  } else {
    dataall2 <- rbind(dataall2, c(mmat, dname[i]))
    rmindex <- c(rmindex, iindex)
  }
  rownames(dataall2)[nrow(dataall2)] <- paste(rownames(dataall[iindex,]), collapse="|")
}
dataall2 <- dataall2[-rmindex,]

# Compare human controls (HC) vs. orangutan controls (HO)
sampletitle <- "HO_H vs HO_O"
countData <- apply(as.matrix(dataall2[, 1:6]), 2, as.numeric)
rownames(countData) <- rownames(dataall2)
countData <- round(countData, digits = 0)
condition <- factor(c(rep("treat", casenum), rep("untreat", controlnum)))
colData <- data.frame(row.names = colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize = TRUE)), by = "row.names", sort = FALSE)
resdata <- resdata[!is.na(resdata$padj),]
resdata <- resdata[order(resdata$padj),]

# Gene annotation
# PluripotencyGene <- read.table("Pluripotency_Genes.txt", header = FALSE)
MT_coding <- read.table("Human mitochondria genome coding genes.csv", sep = ",", header = TRUE)
MT_asso <- read.table("human_chromosome_mito_genelist.csv", sep = ",", header = TRUE)
resdata$genename <- sapply(resdata$Row.names, function(x) strsplit(x, "\\|")[[1]][2])
resdata$genetype <- NA
resdata$genetype[resdata$genename %in% MT_coding[, 2]] <- "MT_coding"
resdata$genetype[resdata$genename %in% MT_asso[, 2]] <- "MT_associated"
resdata$log2FoldChange <- -resdata$log2FoldChange

#------2 DEGs (Differentially Expressed Genes)
# Filter DEGs based on significance and fold change direction, then export results
diff_gene_p <- resdata[resdata$padj < 0.05,]
diff_gene_p$stat[diff_gene_p$log2FoldChange > 0] <- "Up"   # Mark genes upregulated
diff_gene_p$stat[diff_gene_p$log2FoldChange < 0] <- "Down" # Mark genes downregulated
diff_gene_p <- diff_gene_p[abs(diff_gene_p$log2FoldChange) > 1.5,] # Apply fold change threshold
write.table(diff_gene_p, sep = "\t", row.names = FALSE, quote = FALSE, "HOh_vs_HOo_all_DEGs.txt")
write.table(resdata, sep = "\t", row.names = FALSE, quote = FALSE, "HOh_vs_HOo_allgenes.txt")

#-------3 Volcano plot and heatmap
# Prepare data for plotting, highlighting significant changes
data <- resdata
data$threshold <- as.factor(ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1.5,
                                   ifelse(data$log2FoldChange >= 1.5, 'Up', 'Down'), 'Not'))
data$genetype[!data$Row.names %in% diff_gene_p$Row.names] <- NA  # Annotate significant gene types
data$genename[is.na(data$genetype)] <- NA  # Remove gene names for non-significant genes
data$pos <- rowMeans(data[, 8:(8 + casenum - 1)])  # Calculate means for positive samples
data$neg <- rowMeans(data[, (8 + casenum):(8 + num - 1)])  # Calculate means for negative samples

# Generate a scatter plot with annotated regression line and Pearson correlation
p <- ggplot(data, aes(x = log2(pos), y = log2(neg), colour = threshold)) + geom_point()
p1 <- p + geom_point(aes(color = threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha = 0.4, size = 1)
p3 <- p2 + theme_bw() + labs(x = "Case_genes", y = "Control_genes", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 0.6)
a <- cor.test(data$pos, data$neg, method = "pearson")
p6 <- p5 + annotate("text", x = 14, y = -0.5, label = paste("Pearson's correlation ", round(a$estimate, 3)), size = 3)
ggsave(p6, file = "HO_All_correlation.pdf", width = 8, height = 8)

# Prepare data and generate a volcano plot
vdata <- arrange(data, padj)
vdata$threshold <- as.factor(ifelse(vdata$padj < 0.05 & abs(vdata$log2FoldChange) >= 1.5,
                                    ifelse(vdata$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
p <- ggplot(vdata, aes(x = log2FoldChange, y = -log10(padj), colour = threshold))
p1 <- p + geom_point(aes(color = vdata$threshold)) + geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5)
p2 <- p1 + scale_color_manual(values = c("blue", "grey", "red")) + geom_point(alpha = 0.4, size = 1.2)
p3 <- p2 + theme_bw() + labs(x = "log2(FoldChange)", y = "-log10(p-value)", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
v4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(-2, 110) + xlim(-16, 16)
ggsave(v4, file = paste(sampletitle, "_volcano.pdf", sep = ""), width = 4, height = 3.5)

#----- 4 Heatmap
# Prepare the data matrix with selected columns for differential expression
datamat <- as.matrix(diff_gene_p[, c(8:13)])
rownames(datamat) <- diff_gene_p$genename
datamat = datamat[!is.na(diff_gene_p$genename),]

# Prepare gene annotations for the heatmap
anno_gene2 <- cbind(diff_gene_p$genename, diff_gene_p$genetype)
colnames(anno_gene2) <- c("genename", "genetype")
anno_gene2 <- as.data.frame(anno_gene2)
anno_gene2 <- anno_gene2[!is.na(anno_gene2$genetype),]

# Utilize the ComplexHeatmap package to generate a heatmap
library(ComplexHeatmap)
A <- as.matrix(datamat)  # Convert expression matrix to matrix format
A <- A[!grepl("MT-", rownames(A)),]  # Exclude mitochondrial genes
samples <- rep(c('HO_H', 'HO_O'), c(3, 3))  # Define sample group information

# Normalize data by scaling log-transformed counts
for (i in 1:nrow(A)) {
  A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))
}

# Create the heatmap with specified color gradient and annotations
B <- Heatmap(A,  # Expression matrix
             col = colorRampPalette(c("blue", "white", "red"))(100),  # Define color gradient
             show_row_names = TRUE,  # Display row names
             row_names_gp = gpar(fontsize = 6),  # Set font size for row names
             top_annotation = HeatmapAnnotation(Group = samples,  # Top annotations for sample groups
                                                simple_anno_size = unit(2, 'mm'), 
                                                col = list(Group = c('HO_H' = '#F8766D', 'HO_O' = '#619CFF')),
                                                show_annotation_name = FALSE))  # Group color settings and hide annotation name

# Generate a PDF file for the heatmap
pdf(file = "HO_all_heatmap.pdf", width = 6, height = 8)
print(B)  # Print the heatmap object
dev.off()  # Close the PDF device

#-------5 remove pseudogenes and redraw heatmap
psedogenes<-read.table("psedogenelist4.txt")
diff_gene_p<-diff_gene_p[!diff_gene_p$genename %in% psedogenes$V1,]
datamat<-as.matrix(diff_gene_p[,c(8:13)])
rownames(datamat)<-diff_gene_p$genename
datamat=datamat[!is.na(diff_gene_p$genename),]                    
A <- as.matrix(datamat) 
glist<-read.table("H-O_coding_DEGs.txt",header = F)
A<-A[rownames(A)%in%glist$V1,]
for (i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 2))
B <- Heatmap(A,
             cluster_columns = F,
             col = colorRampPalette(c("blue","white","red"))(100),
             show_row_names = T,
             row_names_gp = gpar(fontsize = 6))
pdf(file = "HO_coding_genes_heatmap.pdf",width = 5,height = 6)
print(B)
dev.off()

#-------6 mitocondra gene expression level compare
omt<-read.table("MT_gene_ids.txt",header = F)
oname<-sapply(rownames(datao),function(x) strsplit(x, "\\|")[[1]][1])
mtindex<-oname%in%omt$V1
for (i in 1:nrow(datao[mtindex,])) {
  mn<-rownames(datao)[mtindex]
  i_name<-strsplit(mn[i],"\\|")[[1]][2]
  if (is.na(i_name)) {
    iname<-paste(mn[i],"|MT-",i,sep = "")
  }else{
    iname<-paste(strsplit(mn[i],"\\|")[[1]][1],"|MT-",i_name,sep = "")
  }
  rownames(datao)[mtindex][i]<-iname
}
dhname=sapply(rownames(datah),function(x) strsplit(x, "\\|")[[1]][2])
doname=sapply(rownames(datao),function(x) strsplit(x, "\\|")[[1]][2])
commonname<-dhname[dhname%in%doname]
commonname<-commonname[!is.na(commonname)]
datah<-as.data.frame(cbind(datah,dhname))
datao<-as.data.frame(cbind(datao,doname))
datah<-datah[dhname%in%commonname,]
datao<-datao[doname%in%commonname,]
colnames(datah)<-c("hHO_h1","hHO_h2","hHO_h3","hHO_o1","hHO_o2","hHO_o3","genename")
colnames(datao)<-c("oHO_h1","oHO_h2","oHO_h3","oHO_o1","oHO_o2","oHO_o3","genename")
datah<-datah[!is.na(datah$genename),]
dname<-unique(datah$genename[duplicated(datah$genename)])
for (i in 1:length(dname)) {
  iindex=which(datah$genename==dname[i])
  mmat=colMeans(apply(datah[iindex,1:6], 2, as.numeric))
  if (i==1) {
    datah2<-rbind(datah,c(mmat,dname[i]))
    rmindex<-iindex
  }else{
    datah2<-rbind(datah2,c(mmat,dname[i]))
    rmindex<-c(rmindex,iindex)
  }
  rownames(datah2)[nrow(datah2)] <- paste(rownames(datah[iindex,]), collapse="|")
}
datah2<-datah2[-rmindex,]
datao<-datao[!is.na(datao$genename),]
dname<-unique(datao$genename[duplicated(datao$genename)])
for (i in 1:length(dname)) {
  iindex=which(datao$genename==dname[i])
  mmat=colMeans(apply(datao[iindex,1:6], 2, as.numeric))
  if (i==1) {
    datao2<-rbind(datao,c(mmat,dname[i]))
    rmindex<-iindex
  }else{
    datao2<-rbind(datao2,c(mmat,dname[i]))
    rmindex<-c(rmindex,iindex)
  }
  rownames(datao2)[nrow(datao2)] <- paste(rownames(datao[iindex,]), collapse="|")
}
datao2<-datao2[-rmindex,]
datao2<-datao2[order(datao2$genename),]
datah2<-datah2[order(datah2$genename),]
# merge all
if (identical(datah2$genename,datao2$genename)) {
  data_deseq<-cbind(datah2,datao2)
}else{
  print("Not same!")
}
rownames(data_deseq)<-NULL
countData <- apply(as.matrix(data_deseq[,c(1:6,8:13)]), 2, as.numeric)
rownames(countData)<-data_deseq$genename
countData<-round(countData,digits = 0)
condition <- factor(c(rep("treat",6),rep("untreat",6)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
dataall_deseq2 <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
mtindex=which(grepl("MT-",dataall_deseq2$Row.names))
# calculate MT gene counts
datahnew_mt    <- colSums(dataall_deseq2[mtindex,8:13])
dataonew_mt    <- colSums(dataall_deseq2[mtindex,14:19])
mtdraw<-rbind(cbind(datahnew_mt,rownames(datahnew_mt),c(rep("HOh",3),rep("HOo",3)),rep("Human",6)),
              cbind(datacnew_mt,rownames(dataonew_mt),c(rep("HOh",3),rep("HOo",3)),rep("Orangutan",6)))
mtdraw<-as.data.frame(mtdraw)
colnames(mtdraw)<-c("count","sample","specie")
mtdraw$count<-as.numeric(mtdraw$count)
library(dplyr)
summary_mt <- mtdraw %>%
  group_by(sample, specie) %>%
  summarise(
    mean_count = mean(count),
    se = sd(count) / sqrt(3),  
    .groups = 'drop'
  )
summary_mt$sample<-factor(summary_mt$sample,levels=c("HOh","HOo"))
summary_mt$specie<-factor(summary_mt$specie,levels=c("Human","Orangutan"))
ce<-ddply(summary_mt,"sample",transform,percent_weight=mean_count/sum(mean_count)*100)
ce<-as.data.frame(ce)
p_bar2<-ggplot(ce,aes(x=percent_weight,y=sample,fill=specie))+geom_bar(stat = "identity")+
  theme_bw()+coord_flip()+ scale_x_continuous(limits = c(0, NA))+
  scale_fill_manual(values = c("Human"="#1E88E5","Orangutan"="#FF8000"))
p_bar3<-p_bar2+ geom_text(aes(label = round(percent_weight,2), x = percent_weight),size = 5,color = "black")
ggsave("HO_mt_bar_percent.pdf", p_bar3, width=6,height=6)
