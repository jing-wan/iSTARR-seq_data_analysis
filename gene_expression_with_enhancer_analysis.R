#!/usr/bin/Rscript

require(data.table)
require(ggplot2)
require(ggsignif)
require(clusterProfiler)
require(org.At.tair.db)

#========================enhancer proximal gene expression=============================

gene_fpkm <- fread("gene_fpkm.txt")
colnames(gene_fpkm) <- c("geneID", "fpkm")
#merge multiple transcript of same gene
gene_fpkm <- aggregate(gene_fpkm$fpkm, by = list(geneID = gene_fpkm$geneID), data = gene_fpkm, FUN = mean)
colnames(gene_fpkm) <- c("geneID", "fpkm")

gene_enhancer_fpkm <- fread("gene_enhancer_fpkm.txt")
colnames(gene_enhancer_fpkm) <- c("geneID", "enhancer_count", "fpkm")
#merge multiple transcript of same gene
gene_enhancer_fpkm <- aggregate(gene_enhancer_fpkm$fpkm, 
                                by = list(geneID = gene_enhancer_fpkm$geneID, enhancer_count = gene_enhancer_fpkm$enhancer_count), 
                                data = gene_enhancer_fpkm, FUN = mean)
colnames(gene_enhancer_fpkm) <- c("geneID", "enhancer_count", "fpkm")

fpkm_type1 <- data.frame(fpkm = gene_enhancer_fpkm[gene_enhancer_fpkm$enhancer_count >= 1, "fpkm"], 
                         type = rep("exist enhancer", length(gene_enhancer_fpkm[gene_enhancer_fpkm$enhancer_count >= 1, "fpkm"])))

fpkm_type0 <- data.frame(fpkm = gene_enhancer_fpkm[gene_enhancer_fpkm$enhancer_count == 0, "fpkm"], 
                         type = rep("no enhancer", length(gene_enhancer_fpkm[gene_enhancer_fpkm$enhancer_count == 0, "fpkm"])))

fpkm_control <- data.frame(fpkm = gene_fpkm$fpkm, type = rep("control", length(gene_fpkm$fpkm)))

fpkm_type_df <- rbind(fpkm_type1, fpkm_type0, fpkm_control)

fpkm_type_df$fpkm <- log(fpkm_type_df$fpkm+1)

#reference line
control_df <- fpkm_type_df[fpkm_type_df$type %in% "control", ]
control_mean <- mean(control_df$fpkm)
control_median <- median(control_df$fpkm)


pdf("Gene_expression_with_enhancer_exist_not.pdf", width = 8, height = 6)

ggplot(data = fpkm_type_df, mapping = aes(x = type, y = fpkm), color = type)+
  geom_violin(aes(fill=type), width = 1)+
  geom_boxplot(aes(fill=type), width = 0.05, outlier.shape = NA)+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "gray"))+
  geom_hline(yintercept = control_median, linetype = 2, color = "gray")+
  ylim(0, 7.5)+
  labs(title = "Gene expression with enhancer", 
       x = element_blank(),  
       y = "log(fpkm+1)",
       fill = "enhancer count")+
  scale_x_discrete(labels = c(paste("exist enhancer\n", nrow(fpkm_type1)), 
                              paste("no enhancer\n", nrow(fpkm_type0)), 
                              paste("gene control\n", nrow(control_df))))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size=15, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("exist enhancer", "no enhancer"), c("no enhancer", "control"), 
                                 c("exist enhancer", "control")), 
              y_position = c(6.5, 7, 7.5))

dev.off()


#===============================gene_with_enhancer_GO_analysis================================

gene_enhancer1_fpkm_0 <- fread("gene_enhancer1_fpkm_0.txt", header = F)
gene_enhancer1_fpkm_0_1 <- fread("gene_enhancer1_fpkm_0_1.txt", header = F)
gene_enhancer1_fpkm_1_10 <- fread("gene_enhancer1_fpkm_1_10.txt", header = F)
gene_enhancer1_fpkm_over10 <- fread("gene_enhancer1_fpkm_over10.txt", header = F)

gene_enhancer1_fpkm_0_GO <- enrichGO(gene_enhancer1_fpkm_0$V1, OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pvalueCutoff = 0.05)
gene_enhancer1_fpkm_0_1_GO <- enrichGO(gene_enhancer1_fpkm_0_1$V1, OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pvalueCutoff = 0.05)
gene_enhancer1_fpkm_1_10_GO <- enrichGO(gene_enhancer1_fpkm_1_10$V1, OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pvalueCutoff = 0.05)
gene_enhancer1_fpkm_over10_GO <- enrichGO(gene_enhancer1_fpkm_over10$V1, OrgDb = "org.At.tair.db", keyType="TAIR", ont="ALL", pvalueCutoff = 0.05)


p_enhancer1 <- dotplot(merge_result(list(FPKM0 = gene_enhancer1_fpkm_0_GO@result, 
                                         FPKM0_1 = gene_enhancer1_fpkm_0_1_GO@result, 
                                         FPKM1_10 = gene_enhancer1_fpkm_1_10_GO@result, 
                                         FPKM_over10 = gene_enhancer1_fpkm_over10_GO@result)), 
                       showCategory = 5, font.size=15, title = "GO results of different expression genes with enhancer")

pdf("Different_expression_gene_with_enhancer_GO_results.pdf", width = 12, height = 7)
print(p_enhancer1)
dev.off()



#========================gene_classify_TSS_TTS_enhancer=============================

gene_TSS_enhancer <- fread("gene_TSS_enhancer_coverage.bed")
gene_TTS_enhancer <- fread("gene_TTS_enhancer_coverage.bed")

gene_TSS_TTS_enhancer <- cbind(gene_TSS_enhancer[, 1:6], gene_TTS_enhancer[, 1:6])
colnames(gene_TSS_TTS_enhancer) <- c("chr1", "start1", "end1", "strand1", "geneID1", "enhancer1", 
                                     "chr2", "start2", "end2", "strand2", "geneID2", "enhancer2")

gene_fpkm <- fread("gene_fpkm.txt")
colnames(gene_fpkm) <- c("geneID", "fpkm")

#merge multiple transcript of same gene
gene_fpkm <- aggregate(gene_fpkm$fpkm, by = list(geneID = gene_fpkm$geneID), data = gene_fpkm, FUN = mean)
colnames(gene_fpkm) <- c("geneID1", "fpkm")

gene_TSS_TTS_enhancer_fpkm <- merge(gene_TSS_TTS_enhancer, gene_fpkm, by = "geneID1", all = T)

gene_TSS_TTS_enhancer_fpkm <- na.omit(gene_TSS_TTS_enhancer_fpkm)

control_gene <- data.frame(genetype = rep("control", length(gene_fpkm$fpkm)), fpkm = gene_fpkm$fpkm)

#gene classfied by enhancer in TSS-TTS region
gene_TSS1_TTS1_enhancer_fpkm <- gene_TSS_TTS_enhancer_fpkm[gene_TSS_TTS_enhancer_fpkm$enhancer1 >= 1 & gene_TSS_TTS_enhancer_fpkm$enhancer2 >= 1, ]
gene_TSS1_TTS0_enhancer_fpkm <- gene_TSS_TTS_enhancer_fpkm[gene_TSS_TTS_enhancer_fpkm$enhancer1 >= 1 & gene_TSS_TTS_enhancer_fpkm$enhancer2 == 0, ]
gene_TSS0_TTS1_enhancer_fpkm <- gene_TSS_TTS_enhancer_fpkm[gene_TSS_TTS_enhancer_fpkm$enhancer1 == 0 & gene_TSS_TTS_enhancer_fpkm$enhancer2 >= 1, ]
gene_TSS0_TTS0_enhancer_fpkm <- gene_TSS_TTS_enhancer_fpkm[gene_TSS_TTS_enhancer_fpkm$enhancer1 == 0 & gene_TSS_TTS_enhancer_fpkm$enhancer2 == 0, ]

gene_TSS1_TTS1_fpkm <- data.frame(genetype = rep("TSS1_TTS1", nrow(gene_TSS1_TTS1_enhancer_fpkm)), fpkm = gene_TSS1_TTS1_enhancer_fpkm$fpkm)
gene_TSS1_TTS0_fpkm <- data.frame(genetype = rep("TSS1_TTS0", nrow(gene_TSS1_TTS0_enhancer_fpkm)), fpkm = gene_TSS1_TTS0_enhancer_fpkm$fpkm)
gene_TSS0_TTS1_fpkm <- data.frame(genetype = rep("TSS0_TTS1", nrow(gene_TSS0_TTS1_enhancer_fpkm)), fpkm = gene_TSS0_TTS1_enhancer_fpkm$fpkm)
gene_TSS0_TTS0_fpkm <- data.frame(genetype = rep("TSS0_TTS0", nrow(gene_TSS0_TTS0_enhancer_fpkm)), fpkm = gene_TSS0_TTS0_enhancer_fpkm$fpkm)

gene_type_fpkm_control <- rbind(gene_TSS1_TTS1_fpkm, gene_TSS1_TTS0_fpkm, gene_TSS0_TTS1_fpkm, gene_TSS0_TTS0_fpkm, control_gene)

gene_type_fpkm_control$fpkm <- log(gene_type_fpkm_control$fpkm+1)

control_df <- gene_type_fpkm_control[gene_type_fpkm_control$genetype %in% "control", ]
control_mean <- mean(control_df$fpkm)
control_median <- median(control_df$fpkm)

gene_type_fpkm_control$genetype <- factor(gene_type_fpkm_control$genetype, levels = c("TSS1_TTS1", "TSS1_TTS0", "TSS0_TTS1", "TSS0_TTS0", "control"), ordered = T)

pdf("gene_expression_classify_by_TSS-TTS_enhancer-new-pvalue.pdf", width = 9, height = 6)

ggplot(data = gene_type_fpkm_control, mapping = aes(x = genetype, y = fpkm), color = genetype)+
  geom_violin(aes(fill=factor(genetype)), width = 1)+
  geom_boxplot(aes(fill=factor(genetype)), width = 0.05, outlier.shape = NA)+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "gray"))+
  geom_hline(yintercept = control_median, linetype = 2, color = "gray")+
  ylim(0, 10)+
  labs(title = "gene expression classify by TSS-TTS enhancer", 
       x = "gene type", y = "log(fpkm+1)",
       fill = "gene type")+
  scale_x_discrete(labels = c(paste("TSS1_TTS1\n", nrow(gene_TSS1_TTS1_enhancer_fpkm)), 
                              paste("TSS1_TTS0\n", nrow(gene_TSS1_TTS0_enhancer_fpkm)), 
                              paste("TSS0_TTS1\n", nrow(gene_TSS0_TTS1_enhancer_fpkm)),
                              paste("TSS0_TTS0\n", nrow(gene_TSS0_TTS0_enhancer_fpkm)),
                              paste("Gene control\n", nrow(control_df))))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size=15, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15))+
  guides(fill = FALSE)

dev.off()



#============================clustered_enhancer_gene_expression===========================

enhancer_4cluster_gene_fpkm <- data.frame(gene_ID = NA, fpkm = NA, enhancer_type = NA)

for (i in 1:4) {
  enhancer_cluster_file <- paste("enh_epi15_cluster", i, "_closest_gene.bed", sep = "")
  
  gene_enhancer <- fread(enhancer_cluster_file, select = c(9))
  colnames(gene_enhancer) <- c("geneID")
  
  gene_fpkm <- fread("gene_fpkm.txt")
  colnames(gene_fpkm) <- c("geneID", "fpkm")
  gene_fpkm <- aggregate(gene_fpkm$fpkm, by = list(geneID = gene_fpkm$geneID), data = gene_fpkm, FUN = mean)
  colnames(gene_fpkm) <- c("geneID", "fpkm")
  
  gene_enhancer_fpkm <- merge(gene_enhancer, gene_fpkm, by = "geneID", all.x = T)
  gene_enhancer_fpkm[is.na(gene_enhancer_fpkm$fpkm)==T, "fpkm"] <- 0

  enhancer_clustered_gene_fpkm <- data.frame(gene_ID = gene_enhancer_fpkm$geneID,
                                             fpkm = gene_enhancer_fpkm$fpkm,
                                             enhancer_type = rep(paste("cluster", i, sep = ""), nrow(gene_enhancer_fpkm)))
  
  enhancer_4cluster_gene_fpkm <- rbind(enhancer_4cluster_gene_fpkm, enhancer_clustered_gene_fpkm)
}
enhancer_4cluster_gene_fpkm <- enhancer_4cluster_gene_fpkm[-1, ]

enhancer_4cluster_gene_fpkm <- unique(enhancer_4cluster_gene_fpkm)

write.table(enhancer_4cluster_gene_fpkm, file = "iSTARR_enhancer_4cluster_nearest_gene_fpkm.txt", col.names = T, row.names = F, quote = F)

 

all_enhancer_clustered_gene_fpkm <- fread("iSTARR_enhancer_4cluster_nearest_gene_fpkm.txt")

gene_fpkm <- aggregate(all_enhancer_clustered_gene_fpkm$fpkm, 
                       by = list(geneID = all_enhancer_clustered_gene_fpkm$gene_ID, enhancer_type = all_enhancer_clustered_gene_fpkm$enhancer_type), 
                       data = all_enhancer_clustered_gene_fpkm, 
                       FUN = mean)
colnames(gene_fpkm) <- c("geneID", "fpkm", "enhancer_type")

gene_fpkm <- fread("gene_fpkm.txt")
colnames(gene_fpkm) <- c("geneID", "fpkm")
gene_fpkm <- aggregate(gene_fpkm$fpkm, by = list(geneID = gene_fpkm$geneID), data = gene_fpkm, FUN = mean)
colnames(gene_fpkm) <- c("geneID", "fpkm")

control_gene <- data.frame(fpkm = gene_fpkm$fpkm, enhancer_type = rep("control", length(gene_fpkm$fpkm)))

enh_cluster_gene_control <- rbind(all_enhancer_clustered_gene_fpkm[, 2:3], control_gene)

enh_cluster_gene_control$fpkm <- log(enh_cluster_gene_control$fpkm+1)

control_df <- enh_cluster_gene_control[enh_cluster_gene_control$enhancer_type %in% "control", ]
control_mean <- mean(control_df$fpkm)
control_median <- median(control_df$fpkm)


pdf("iSTARR_epi15_4cluster_enhancer_nearest_gene_expression_pvalue.pdf", width = 9, height = 6)

ggplot(data = enh_cluster_gene_control, mapping = aes(x = enhancer_type, y = fpkm), color = enhancer_type)+
  geom_violin(aes(fill=factor(enhancer_type)), scale = "width", width = 0.6)+
  geom_boxplot(aes(fill=factor(enhancer_type)), width = 0.05, outlier.shape = NA)+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "gray"))+
  geom_hline(yintercept = control_median, linetype = 2, color = "gray")+
  ylim(0, 10)+
  labs(title = "Different iSTARR-seq enhancer cluster nearest gene expression", 
       x = "Different enhancer cluster nearest gene", y = "log(fpkm+1)",
       fill = "enhancer type")+
  scale_x_discrete(labels = c(paste("cluster1\n", nrow(enh_cluster_gene_control[enh_cluster_gene_control$enhancer_type %in% "cluster1", ])), 
                              paste("cluster2\n", nrow(enh_cluster_gene_control[enh_cluster_gene_control$enhancer_type %in% "cluster2", ])), 
                              paste("cluster3\n", nrow(enh_cluster_gene_control[enh_cluster_gene_control$enhancer_type %in% "cluster3", ])),
                              paste("cluster4\n", nrow(enh_cluster_gene_control[enh_cluster_gene_control$enhancer_type %in% "cluster4", ])),
                              paste("Gene control\n", nrow(control_df))))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text=element_text(size=12, color = "black"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15))+
  guides(fill = FALSE)

dev.off()

