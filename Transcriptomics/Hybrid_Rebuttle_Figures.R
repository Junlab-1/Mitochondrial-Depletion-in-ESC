rm(list=ls(all=T))
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(plyr)
library(ComplexHeatmap)
#================================== gene expressions after rmNUMTs ==========================
#============================================================================================
# hybrid samples gene expressions
datahrmnumt<-as.matrix(read.csv(
  "~/orangutan-human/run/rmNUMT/all_genecount/human_ho_gene_nonumt.csv"
  ,row.names="gene_id"))
datahrmnumt<-datahrmnumt[order(rownames(datahrmnumt)),]

dataormnumt<-as.matrix(read.csv(
  "~/orangutan-human/run/rmNUMT/all_genecount/orangutan_ho_gene_nonumt.csv"
  ,row.names="gene_id"))
dataormnumt<-dataormnumt[order(rownames(dataormnumt)),]
omt<-read.table("~/reference/FastaGTF/orangutan/MT_gene_ids.txt",header = F)
oname<-sapply(rownames(dataormnumt),function(x) strsplit(x, "\\|")[[1]][1])
mtindex<-oname%in%omt$V1
mtindex<-oname%in%omt$V1
datahrmnumt_mt    <- colSums(datahrmnumt[grepl("\\|MT-",rownames(datahrmnumt)),])
datahrmnumt_other <- colSums(datahrmnumt[!grepl("\\|MT-",rownames(datahrmnumt)),])
dataormnumt_mt    <- colSums(dataormnumt[mtindex,])
dataormnumt_other <- colSums(dataormnumt[!mtindex,])
for (i in 1:nrow(dataormnumt[mtindex,])) {
  mn<-rownames(dataormnumt)[mtindex]
  i_name<-strsplit(mn[i],"\\|")[[1]][2]
  if (is.na(i_name)) {
    iname<-paste(mn[i],"|MT-",i,sep = "")
  }else{
    iname<-paste(strsplit(mn[i],"\\|")[[1]][1],"|MT-",i_name,sep = "")
  }
  rownames(dataormnumt)[mtindex][i]<-iname
}
dataallrmnumt<-as.data.frame(rbind(datahrmnumt,dataormnumt))
dataallrmnumt=cbind(dataallrmnumt,sapply(rownames(dataallrmnumt),function(x) strsplit(x, "\\|")[[1]][2]))
colnames(dataallrmnumt)<-c("HO_h1","HO_h2","HO_h3","HO_o1","HO_o2","HO_o3","genename")
dataallrmnumt<-dataallrmnumt[!is.na(dataallrmnumt$genename),]
dname<-unique(dataallrmnumt$genename[duplicated(dataall$genename)])
for (i in 1:length(dname)) {
  iindex=which(dataallrmnumt$genename==dname[i])
  mmat=colMeans(apply(dataallrmnumt[iindex,1:6], 2, as.numeric))
  if (i==1) {
    dataallrmnumt2<-rbind(dataallrmnumt,c(mmat,dname[i]))
    rmindex<-iindex
  }else{
    dataallrmnumt2<-rbind(dataallrmnumt2,c(mmat,dname[i]))
    rmindex<-c(rmindex,iindex)
  }
  rownames(dataallrmnumt2)[nrow(dataallrmnumt2)] <- paste(rownames(dataall[iindex,]), collapse="|")
}
dataallrmnumt2<-dataallrmnumt2[-rmindex,]
rownames(dataallrmnumt2)<-dataallrmnumt2$genename
hoallrmnumt <- apply(as.matrix(dataallrmnumt2[,1:6]), 2, as.numeric)
rownames(hoallrmnumt)<-dataallrmnumt2$genename
hoallrmnumt<-hoallrmnumt[!is.na(rownames(hoallrmnumt)),]
hoallrmnumt<-hoallrmnumt[order(rownames(hoallrmnumt)),]

# control samples gene expressions
hcontrol<-as.matrix(read.csv(
  "~/orangutan-human/run/hcontrol/rmNUMT/all_genecount/human_control_gene_nonumt.csv"
  ,row.names="gene_id"))
rownames(hcontrol)<-sapply(rownames(hcontrol),function(x) strsplit(x, "\\|")[[1]][2])
hcontrol<-hcontrol[!is.na(rownames(hcontrol)),]
dname<-unique(rownames(hcontrol)[duplicated(rownames(hcontrol))])
for (i in 1:length(dname)) {
  iindex=which(rownames(hcontrol)==dname[i])
  mmat=colMeans(hcontrol[iindex,])
  if (i==1) {
    hcontrol2<-rbind(hcontrol,mmat)
    rmindex<-iindex
  }else{
    hcontrol2<-rbind(hcontrol,mmat)
    rmindex<-c(rmindex,iindex)
  }
  rownames(hcontrol2)[nrow(hcontrol2)] <- dname[i]
}
hcontrol2<-hcontrol2[-rmindex,]
hcontrol2_rmnumt<-hcontrol2[order(rownames(hcontrol2)),]

ocontrol<-as.matrix(read.csv(
  "~/orangutan-human/run/ocontrol/rmNUMT/all_genecount/orangutan_control_gene_nonumt.csv"
  ,row.names="gene_id"))
oname<-sapply(rownames(ocontrol),function(x) strsplit(x, "\\|")[[1]][1])
mtindex<-oname%in%omt$V1
for (i in 1:nrow(ocontrol[mtindex,])) {
  mn<-rownames(ocontrol)[mtindex]
  i_name<-strsplit(mn[i],"\\|")[[1]][2]
  if (is.na(i_name)) {
    iname<-paste(mn[i],"|MT-",i,sep = "")
  }else{
    iname<-paste(strsplit(mn[i],"\\|")[[1]][1],"|MT-",i_name,sep = "")
  }
  rownames(ocontrol)[mtindex][i]<-iname
}
rownames(ocontrol)<-sapply(rownames(ocontrol),function(x) strsplit(x, "\\|")[[1]][2])
ocontrol<-ocontrol[!is.na(rownames(ocontrol)),]
dname<-unique(rownames(ocontrol)[duplicated(rownames(ocontrol))])
for (i in 1:length(dname)) {
  iindex=which(rownames(ocontrol)==dname[i])
  mmat=colMeans(ocontrol[iindex,])
  if (i==1) {
    ocontrol2<-rbind(ocontrol,mmat)
    rmindex<-iindex
  }else{
    ocontrol2<-rbind(ocontrol,mmat)
    rmindex<-c(rmindex,iindex)
  }
  rownames(ocontrol2)[nrow(ocontrol2)] <- dname[i]
}
ocontrol2<-ocontrol2[-rmindex,]
ocontrol2_rmnumt<-ocontrol2[order(rownames(ocontrol2)),]
hccom<-hoallrmnumt[rownames(hoallrmnumt)%in%rownames(hcontrol2_rmnumt),]
hcom<-hcontrol2_rmnumt[rownames(hcontrol2_rmnumt)%in%rownames(hoallrmnumt),]
if (identical(rownames(hccom),rownames(hcom))) {
  pcamat1<-cbind(hccom,hcom)
}else{
  print("Not same!")
}
hccom2<-pcamat1[rownames(pcamat1)%in%rownames(ocontrol2_rmnumt),]
ccom<-ocontrol2_rmnumt[rownames(ocontrol2_rmnumt)%in%rownames(pcamat1),]
if (identical(rownames(hccom2),rownames(ccom))) {
  pcamat2<-cbind(hccom2,ccom)
}else{
  print("Not same!")
}
hoallwithControlcounts<-pcamat2

# calculated relative counts by deseq2
library(DESeq2)
countData <- apply(as.matrix(hoallwithControlcounts), 2, as.numeric)
rownames(countData)<-rownames(hoallwithControlcounts)
countData<-round(countData,digits = 0)
condition <- factor(c(rep("HOh",3),rep("HOo",3),rep("Hcontrol",3),rep("Ocontrol",3)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
hoallwithControldeseq2<-resdata[,8:19]

# =========== pluripotency genes heatmap on all hybrid and control samples (Figure R1C)=================
PluripotencyGene<-read.table("~/reference/genesets/Pluripotency_genes.txt",header = F)
pluripotencydraw<-hoallwithControldeseq2[rownames(hoallwithControldeseq2)%in%PluripotencyGene$V1,]
pluripotencydraw<-pluripotencydraw[rowSums(pluripotencydraw)>0,]
library(viridis)
pdf(file = "HOwithcontrol_pluripotency_heatmap_deseq2.pdf",width = 6,height = 10)
pheatmap(mat = pluripotencydraw,scale = "row",cluster_cols = F,cluster_rows = T,
         color = viridis_pal()(100)[20:100],border_color = NA,
         display_numbers = F,number_format = "%.3f",fontsize_number = 10)
dev.off()
# =========== calculated DEGs and volcano plot(Figure R1D,E) ============================================
# HOh vs HOo(Figure R1D)
countData <- apply(as.matrix(hoallrmnumt), 2, as.numeric)
rownames(countData)<-rownames(hoall)
countData<-round(countData,digits = 0)
condition <- factor(c(rep("HOh",3),rep("HOo",3)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
hoalldeseq2<-resdata
hoalldeseq2$log2FoldChange<- -hoalldeseq2$log2FoldChange
rownames(hoalldeseq2)<-resdata[,1]
hopdraw<-hoalldeseq2[rownames(hoalldeseq2)%in%PluripotencyGene$V1,]

# volcano plot
hopdraw$threshold <- as.factor(ifelse(hopdraw$padj < 0.05 & abs(hopdraw$log2FoldChange) >=1.5,ifelse(hopdraw$log2FoldChange > 1.5 ,'Up','Down'),'Not'))
p<-ggplot(hopdraw,aes(x=log2FoldChange,y =-log10(padj),colour=threshold))
p1<-p + geom_point(aes(color=hopdraw$threshold))+geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.5) +geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)
p2<-p1 + scale_color_manual(values =c("grey"))+geom_point(alpha=0.4, size=1.2)
p3<- p2+ theme_bw()+labs(x="log2(FoldChange)",y="-log10(p-value)")+theme(plot.title = element_text(hjust = 0.5))
v4<-p3+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ylim(-2, 10) + xlim(-6, 6)
ggsave(v4, file="HOh_vs_HOo_volcano.pdf",width=4,height = 3.5)
write.csv(hopdraw,file = "HOh_vs_HOo_pluripotency.csv")


# Hcontrol vs Ocontrol (Figure R1E)
library(DESeq2)
# calculated relative counts by deseq2
countData <- apply(as.matrix(hoallwithControlcounts[,7:12]), 2, as.numeric)
rownames(countData)<-rownames(hoallwithControlcounts)
countData<-round(countData,digits = 0)
condition <- factor(c(rep("Hcontrol",3),rep("Ocontrol",3)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
hocontroldeseq2<-resdata
hocontroldeseq2$log2FoldChange<- -hocontroldeseq2$log2FoldChange
rownames(hocontroldeseq2)<-resdata[,1]
controlpdraw<-hocontroldeseq2[rownames(hocontroldeseq2)%in%PluripotencyGene$V1,]
# volcano plot
controlpdraw$threshold <- as.factor(ifelse(controlpdraw$padj < 0.05 & abs(controlpdraw$log2FoldChange) >=1.5,ifelse(controlpdraw$log2FoldChange > 1.5 ,'Up','Down'),'Not'))
p<-ggplot(controlpdraw,aes(x=log2FoldChange,y =-log10(padj),colour=threshold))
p1<-p + geom_point(aes(color=controlpdraw$threshold))+geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.5) +geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)
p2<-p1 + scale_color_manual(values =c("blue", "grey", "red"))+geom_point(alpha=0.4, size=1.2)
p3<- p2+ theme_bw()+labs(x="log2(FoldChange)",y="-log10(p-value)")+theme(plot.title = element_text(hjust = 0.5))
v4<-p3+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ylim(-2, 100) + xlim(-10, 10)
ggsave(v4, file="Hcontrol_vs_Ocontrol_volcano.pdf",
       width=4,height = 3.5)
write.csv(controlpdraw,file = "Hcontrol_vs_Ocontrol_pluripotency.csv")

# ===================================== compare MT gene expressions ==============================================
# published data (Figure R2B)
# load PMID38623329 data
h1<-read.csv("GSE243579_rawCounts_hgnc.csv",header = T)
rownames(h1)<-h1$Gene
h1<-h1[,-1]
h1c<-h1[,c("iPSC3_1754_Naive_iPSC_1","iPSC3_1754_Naive_iPSC_2","iPSC3_1754_Naive_iPSC_3")]
source("~/Code/LJ.function.R")
h1c_tpm<-count2tpm(h1c,idType = "symbol")
h1c_mt_tpm    <- colSums(h1c_tpm[grepl("^MT-",rownames(h1c_tpm)),])
h1c_other_tpm <- colSums(h1c_tpm[!grepl("^MT-",rownames(h1c_tpm)),])
# load Chimpanzee
c1<-read.table("GSE155443_chimp_gTPM_and_stats_table.tsv",header = F,sep = "\t")
rownames(c1)<-c1[,1]
c1<-c1[,c(10:13)]
colnames(c1)<-c1[4,]
c1c<-c1[-c(1:4),]
c1c <- as.matrix(c1c)
mode(c1c) <- "numeric"
cmt<-read.table("~/reference/FastaGTF/chimpanzee/MT_gene_ids.txt",header = F)
cmt<-unique(cmt)
cmt_string <- paste(cmt$V1, collapse = "|")
matched_rows <- grep(cmt_string, rownames(c1c), value = TRUE)

c1c_mt_tpm    <- colSums(c1c[matched_rows,])
c1c_other_tpm <- colSums(c1c[!rownames(c1c)%in%matched_rows,])
# load PMID: 38250322
o1<-read.csv("GSE235790_all_feature_pab.txt",sep = "\t",header = T)
rownames(o1)<-o1[,1]
o1c<-o1[,c(2,6,10,12)]
temptpm <- o1c[,3:4] / o1c[,2]
temptpm <- 1e6 * t(t(temptpm) / colSums(temptpm))
temptpm <- temptpm[!duplicated(rownames(temptpm)),]
temptpm <- temptpm[!is.na(rownames(temptpm)),]
temptpm <- temptpm[!rownames(temptpm)==" ",]
temptpm <- temptpm[,!is.na(colnames(temptpm))]
temptpm <- temptpm[,!colnames(temptpm)==" "]
o1cTPM  <- o1c
if (identical(rownames(temptpm),rownames(o1cTPM))) {
  o1cTPM[,3:4]<-temptpm
}else{
  print("Not Same")
}
o1c_mt_tpm    <- colSums(o1cTPM[grepl("MT",o1cTPM[,1]),3:4])
o1c_other_tpm <- colSums(o1cTPM[!grepl("MT",o1cTPM[,1]),3:4])

#TPMdraw
publicdraw<-rbind(cbind(c(h1c_mt_tpm,h1c_other_tpm,h1c_mt_tpm+h1c_other_tpm),
                       c(rep("Mito_tpm",3),rep("Other_tpm",3),rep("All_tpm",3)),rep("Human",9)),
                  cbind(c(c1c_mt_tpm,c1c_other_tpm,c1c_mt_tpm+c1c_other_tpm),
                        c(rep("Mito_tpm",4),rep("Other_tpm",4),rep("All_tpm",4)),rep("Chimpanzee",12)),
                 cbind(c(o1c_mt_tpm,o1c_other_tpm,o1c_mt_tpm+o1c_other_tpm),
                       c(rep("Mito_tpm",2),rep("Other_tpm",2),rep("All_tpm",2)),rep("Orangutan",6)))
publicdraw<-as.data.frame(publicdraw)
colnames(publicdraw)<-c("tpm","sample","specie")
publicdraw$tpm<-as.numeric(publicdraw$tpm)
library(dplyr)
summary_public <- publicdraw %>%
  group_by(sample, specie) %>%
  summarise(
    mean_tpm = mean(tpm),
    n = n(),
    se = sd(tpm) / sqrt(n),
    .groups = 'drop'
  )
summary_public$sample<-factor(summary_public$sample,levels=c("All_tpm","Other_tpm","Mito_tpm"))
summary_public$specie<-factor(summary_public$specie,levels=c("Human","Chimpanzee","Orangutan"))
p<-ggplot(summary_public,aes(x=sample,y=mean_tpm,fill=specie))+
  geom_bar(stat = "identity",width = 0.6,position=position_dodge(0.7))+
  geom_errorbar(aes(ymin = mean_tpm - se, ymax = mean_tpm + se), width = 0.2, position = position_dodge(0.7))+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Human"="#1E88E5","Chimpanzee"="#22B573","Orangutan"="#FF6000"))
pdf(file = "Published_MTcontrol_tpm_bar.pdf",width = 7,height = 4)
print(p)
dev.off()

hcontrol<-as.matrix(read.csv("~/hcontrolall_genecount/gene_tpm_matrix.csv",row.names="gene_id"))
hcontrol_mt_tpm    <- colSums(hcontrol[grepl("\\|MT-",rownames(hcontrol)),])
hcontrol_other_tpm <- colSums(hcontrol[!grepl("\\|MT-",rownames(hcontrol)),])
ccontrol<-as.matrix(read.csv("~/ccontrolall_genecount/gene_tpm_matrix.csv",row.names="gene_id"))
cmt<-read.table("~/reference/FastaGTF/chimpanzee/MT_gene_ids.txt",header = F)
cname<-sapply(rownames(ccontrol),function(x) strsplit(x, "\\|")[[1]][1])
mtindex<-cname%in%cmt$V1
for (i in 1:nrow(ccontrol[mtindex,])) {
  mn<-rownames(ccontrol)[mtindex]
  i_name<-strsplit(mn[i],"\\|")[[1]][2]
  if (is.na(i_name)) {
    iname<-paste(mn[i],"|MT-",i,sep = "")
  }else{
    iname<-paste(strsplit(mn[i],"\\|")[[1]][1],"|MT-",i_name,sep = "")
  }
  rownames(ccontrol)[mtindex][i]<-iname
}
ccontrol_mt_tpm    <- colSums(ccontrol[grepl("\\|MT-",rownames(ccontrol)),])
ccontrol_other_tpm <- colSums(ccontrol[!grepl("\\|MT-",rownames(ccontrol)),])
ocontrol<-as.matrix(read.csv("~/ocontrolall_genecount/gene_tpm_matrix.csv",row.names="gene_id"))
omt<-read.table("~/reference/FastaGTF/orangutan/MT_gene_ids.txt",header = F)
oname<-sapply(rownames(ocontrol),function(x) strsplit(x, "\\|")[[1]][1])
mtindex<-oname%in%omt$V1
for (i in 1:nrow(ocontrol[mtindex,])) {
  mn<-rownames(ocontrol)[mtindex]
  i_name<-strsplit(mn[i],"\\|")[[1]][2]
  if (is.na(i_name)) {
    iname<-paste(mn[i],"|MT-",i,sep = "")
  }else{
    iname<-paste(strsplit(mn[i],"\\|")[[1]][1],"|MT-",i_name,sep = "")
  }
  rownames(ocontrol)[mtindex][i]<-iname
}
ocontrol_mt_tpm    <- colSums(ocontrol[grepl("\\|MT-",rownames(ocontrol)),])
ocontrol_other_tpm <- colSums(ocontrol[!grepl("\\|MT-",rownames(ocontrol)),])
# our old data MT-TPM draw(Figure R2C)
otherdraw<-rbind(cbind(c(hcontrol_mt_tpm,hcontrol_other_tpm,hcontrol_mt_tpm+hcontrol_other_tpm),
                       c(rep("Mito_tpm",3),rep("Other_tpm",3),rep("All_tpm",3)),rep("Human",9)),
                 cbind(c(ccontrol_mt_tpm,ccontrol_other_tpm,ccontrol_mt_tpm+ccontrol_other_tpm),
                       c(rep("Mito_tpm",3),rep("Other_tpm",3),rep("All_tpm",3)),rep("Chimpanzee",9)),
                 cbind(c(ocontrol_mt_tpm,ocontrol_other_tpm,ocontrol_mt_tpm+ocontrol_other_tpm),
                       c(rep("Mito_tpm",3),rep("Other_tpm",3),rep("All_tpm",3)),rep("Orangutan",9)))

otherdraw<-as.data.frame(otherdraw)
colnames(otherdraw)<-c("tpm","sample","specie")
otherdraw$tpm<-as.numeric(otherdraw$tpm)
library(dplyr)
summary_other <- otherdraw %>%
  group_by(sample, specie) %>%
  summarise(
    mean_tpm = mean(tpm),
    se = sd(tpm) / sqrt(3),
    .groups = 'drop'
  )
summary_other$sample<-factor(summary_other$sample,levels=c("All_tpm","Other_tpm","Mito_tpm"))
summary_other$specie<-factor(summary_other$specie,levels=c("Human","Chimpanzee","Orangutan"))
p<-ggplot(summary_other,aes(x=sample,y=mean_tpm,fill=specie))+
  geom_bar(stat = "identity",width = 0.6,position=position_dodge(0.7))+
  geom_errorbar(aes(ymin = mean_tpm - se, ymax = mean_tpm + se), width = 0.2, position = position_dodge(0.7))+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Human"="#1E88E5","Chimpanzee"="#22B573","Orangutan"="#FF6000"))
pdf(file = "MTcontrol_tpm_bar.pdf",width = 8,height = 4)
print(p)
dev.off()
# rmNUMT MT-coding gene Figure R2E,F
mtdraw<-rbind(cbind(datahrmnumt_mt,rownames(datahrmnumt_mt),c(rep("HOh",3),rep("HOo",3)),rep("Human",6)),
              cbind(dataormnumt_mt,rownames(dataormnumt_mt),c(rep("HOh",3),rep("HOo",3)),rep("Orangutan",6)))
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
summary_mt
p<-ggplot(summary_mt,aes(x=sample,y=mean_count,fill=specie))+
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 2)))+
  geom_bar(stat = "identity",width = 0.6,position=position_dodge(0.7))+
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), width = 0.2, position = position_dodge(0.7))+
  theme_bw()+scale_fill_manual(values = c("Human"="#1E88E5","Orangutan"="#FF8000"))+
  ggtitle("HO hybrid MT coding genes TPM")
pdf(file = "HO_rmNUMT_mt_bar.pdf",width = 4,height = 4)
print(p)
dev.off()
# percent figures
ce<-ddply(summary_mt,"sample",transform,percent_weight=mean_count/sum(mean_count)*100)
ce<-as.data.frame(ce)
p_bar2<-ggplot(ce,aes(x=percent_weight,y=sample,fill=specie))+geom_bar(stat = "identity")+
  theme_bw()+coord_flip()+ scale_x_continuous(limits = c(0, NA))+
  scale_fill_manual(values = c("Human"="#1E88E5","Orangutan"="#FF8000"))
p_bar3<-p_bar2+ geom_text(aes(label = round(percent_weight,2), x = percent_weight),size = 5,color = "black")
ggsave("HO_mt_bar_percent.pdf", p_bar3, width=6,height=6)

#Nuclear coding genes
otherdraw<-rbind(cbind(datahrmnumt_other,rownames(datahrmnumt_other),c(rep("HOh",3),rep("HOo",3)),rep("Human",6)),
              cbind(dataormnumt_other,rownames(dataormnumt_other),c(rep("HOh",3),rep("HOo",3)),rep("Orangutan",6)))
otherdraw<-as.data.frame(otherdraw)
colnames(otherdraw)<-c("count","sample","specie")
otherdraw$count<-as.numeric(otherdraw$count)
library(dplyr)
summary_other <- otherdraw %>%
  group_by(sample, specie) %>%
  summarise(
    mean_count = mean(count),
    se = sd(count) / sqrt(3),
    .groups = 'drop'
  )
summary_other$sample<-factor(summary_other$sample,levels=c("HOh","HOo"))
summary_other$specie<-factor(summary_other$specie,levels=c("Human","Orangutan"))
p<-ggplot(summary_other,aes(x=sample,y=mean_count,fill=specie))+
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05), add = c(0, 2)))+
  geom_bar(stat = "identity",width = 0.6,position=position_dodge(0.7))+
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), width = 0.2, position = position_dodge(0.7))+
  theme_bw()+scale_fill_manual(values = c("Human"="#1E88E5","Orangutan"="#FF8000"))+
  ggtitle("HO hybrid nuclear coding genes TPM")
pdf(file = "HO_rmNUMT_other_bar.pdf",width = 4,height = 4)
print(p)
dev.off()
write.csv(mtdraw, file = "MTTPMsumary.csv")
write.csv(otherdraw, file = "NuclearTPMsummary.csv")
# mapping rate
maprate<-read.table(file = "mappingrates.csv",header = T,sep = ",")
maprate<-as.data.frame(maprate)
maprate$MappingRates<-as.numeric(maprate$MappingRates)
write.csv(maprate, file = "Maprate.csv")
summary_mr <- maprate %>%
  group_by(reference, Sample) %>%
  summarise(
    mean_rate = mean(MappingRates),
    sem = sd(MappingRates) / sqrt(3),
    .groups = 'drop'
  )
p<-ggplot(summary_mr, aes(x = Sample, y = mean_rate, fill = reference)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_rate - sem, ymax = mean_rate + sem),
                width = 0.2, position = position_dodge(width = 0.7),color="black") +
  labs(x = "Group", y = "Reference Map Rate") + scale_y_continuous(expand = c(0, 0),limits = c(0, 1.0))+
  theme_bw()+scale_fill_manual(values =c("Human"="#1E87E4","Orangutan"="#FF8611"))
ggsave(p,filename = "HOmaprate_bar.pdf",width = 4,height = 3)

# ================ All Hybrid samples correlation analysis (Figure R3)================================= 
library(ggplot2)
library(ggbiplot)
# calculated relative counts by deseq2
countData <- apply(as.matrix(hoallwithControlcounts), 2, as.numeric)
rownames(countData)<-rownames(hoallwithControlcounts)
countData<-round(countData,digits = 0)
condition <- factor(c(rep("HOh",3),rep("HOo",3),rep("Hcontrol",3),rep("Ocontrol",3)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
hoallwithControldeseq2<-resdata[,8:19]
rownames(hoallwithControldeseq2)<-resdata[,1]
write.csv(hoallwithControldeseq2,file = "HoallwithControlDeseq2_exp.csv")
# correlaton plot
allmat<-as.data.frame(hoallwithControldeseq2)
library(GGally)
par(mar=c(4, 4, 2,2)+0.1)
pdf("HO_correlation_deseq2.pdf",width=10,height =10)
p<-ggpairs(allmat,
           axisLabels="none",
           aes(alpha = 0.2),
           columns =1:12)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(p)
dev.off()
# correlation heatmap
allcormat<-cor(allmat)
library(viridis)
pdf(file = "HO_correlation_heatmap_deseq2.pdf",width = 11,height = 10)
pheatmap(mat = allcormat,scale = "none",cluster_cols = F,cluster_rows = F,
         color = viridis_pal()(100)[40:100],border_color = "white",
         display_numbers = TRUE,number_format = "%.3f",fontsize_number = 10)
dev.off()
# ===================== PCA figures for HOh and HOo (Figure R4A)=====================================
countData <- apply(as.matrix(hoallrmnumt), 2, as.numeric)
rownames(countData)<-rownames(hoallrmnumt)
countData<-round(countData,digits = 0)
condition <- factor(c(rep("HOh",3),rep("HOo",3)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
hoallrmnumtdeseq2<-resdata[,8:13]
pcamat_hoallrmnumtdeseq2<-as.data.frame(t(hoallrmnumtdeseq2))
pcamat_hoallrmnumtdeseq2 <- pcamat_hoallrmnumtdeseq2[, apply(pcamat_hoallrmnumtdeseq2, 2, function(x) var(x, na.rm = TRUE) != 0)]
pcamatgroup<-c(rep("HO_H",3),rep("HO_O",3))
pcamat_ho <- prcomp(pcamat_hoallrmnumtdeseq2,scale=T,center=T)
library(ggbiplot)
pcamat_ho$x <- scale(pcamat_ho$x)
pdf(file = "H-O_rmnumtdeseq2_pca_k_cluster.pdf",width = 4,height = 4)
ggbiplot(pcamat_ho,
         var.axes=F,
         obs.scale = 1,
         groups = as.factor(k_clusters$cluster),
         ellipse = T,
         circle = F)+
  geom_text(
    aes(label=rownames(pcamat_ho$x)),
    vjust=1.5,
    size=2
  ) +theme_bw() + scale_color_manual(values = c("HO_O"="orange","HO_H"="#008080"))
dev.off()

metadata <- data.frame(
  sample = rownames(pcamat_ho$x),
  group = c(rep("HOh",3),rep("HOo",3))
)
# All sample PCA 
countData <- apply(as.matrix(hoallwithControlcounts), 2, as.numeric)
rownames(countData)<-rownames(hoallwithControlcounts)
countData<-round(countData,digits = 0)
condition <- factor(c(rep("HOh",3),rep("HOo",3),rep("hcontrol",3),rep("ocontrol",3)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
hoallnumtdeseq2<-resdata[,8:19]
rownames(hoallnumtdeseq2)<-resdata[,1]
pcamat_filtered<-as.data.frame(t(hoallnumtdeseq2))
pcamat_filtered <- pcamat_filtered[, apply(pcamat_filtered, 2, function(x) var(x, na.rm = TRUE) != 0)]
pca_result <- prcomp(pcamat_filtered,scale=T,center=T)
library(ggbiplot)
pca_result$x <- scale(pca_result$x)
pdf(file = "H-O_allsample_pca_rmnumt.pdf",width = 4,height = 4)
pcamatgroup<-c(rep("HO_H",3),rep("HO_O",3),rep("Human_control",3),rep("O_control",3))
ggbiplot(pca_result,
         var.axes=F,
         obs.scale = 1,
         groups = pcamatgroup,
         ellipse = T,
         circle = F)+
  geom_text(
    aes(label=rownames(pcamat_filtered)),
    vjust=1.5,
    size=2  
  ) +theme_bw() + scale_color_manual(values = c("HO_O"="orange","HO_H"="#008080",
                                                "Human_control"="blue","O_control"="#F16821"))
dev.off()
