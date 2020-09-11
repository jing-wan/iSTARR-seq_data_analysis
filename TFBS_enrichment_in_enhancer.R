#!/usr/bin/Rscript

require(data.table)
require(ggplot2)

#====================TF enrichment in enhancer===============

#--------------------random as control--------------------------

TAIR10_chr_length <- read.table("TAIR10.chr.sizes_no_chrC_chrM", header = F)
chr_stat <- fread("STARR-seq_peaks.bed", select = 1)
chr_stat <- as.data.frame(table(chr_stat))

#randomly selected genomic sites (same number as enhancers) as control group
sample_sites_genome <- data.frame(chr = 0, sample_site = 0)
for (i in 1:5) {
  set.seed(2019)
  sample_sites <- base::sample(TAIR10_chr_length[i, 2], chr_stat[i, 2])
  sample_sites_chr <- data.frame(chr = rep(TAIR10_chr_length[i, 1], length(sample_sites)), sample_site = sample_sites)
  sample_sites_genome <- rbind(sample_sites_genome, sample_sites_chr)
}
sample_sites_genome <- sample_sites_genome[-1, ]

random_region_genome <- data.frame(chr = sample_sites_genome$chr, start = sample_sites_genome$sample_site, end = sample_sites_genome$sample_site+499)

random_region_genome <- random_region_genome[order(random_region_genome$chr, random_region_genome$start), ]

write.table(random_region_genome, file = "random_control_for_STARR-seq.bed", sep = "\t",col.names = F, row.names = F, quote = F)


#---------------------TF enrichment-----------------------

STARR_control_TFs <- merge(STARR_ol_TFs_fq, control_ol_TFs_fq, by = "TFs", all = T)
colnames(STARR_control_TFs) <- c("TFs", "STARR", "control")

STARR_control_TFs <- STARR_control_TFs[order(STARR_control_TFs$STARR, decreasing = T), ]

STARR_control_TFs_enrich <- data.frame(STARR_control_TFs, enrichment = log(STARR_control_TFs$STARR/STARR_control_TFs$control))

write.table(STARR_control_TFs_enrich, file = "STARR_control_TFs_family_enrich.txt", sep = "\t", col.names = T, row.names = F, quote = F)


#-----------------fisher.test---------------

iSTARR_oldSTARR_43TF <- read.table("iSTARR_oldSTARR_43TF_enrich_pvalue_plabel.txt", header = T, stringsAsFactors = F)

iSTARR_oldSTARR_43TF_pvalue <- data.frame(iSTARR_oldSTARR_43TF, 
                                          iSTARR_pvalue = rep(0, nrow(iSTARR_oldSTARR_43TF)),
                                          oldSTARR_pvalue = rep(0, nrow(iSTARR_oldSTARR_43TF)),
                                          stringsAsFactors = F)

for (i in 1:nrow(iSTARR_oldSTARR_43TF_pvalue)) {
  iSTARR_fisher_res_i <- fisher.test(data.frame(iSTARR = c(4956, iSTARR_oldSTARR_43TF_pvalue[i, "iSTARR"]), 
                                                control = c(4956, iSTARR_oldSTARR_43TF_pvalue[i, "iSTARR_control"])))
  oldSTARR_fisher_res_i <- fisher.test(data.frame(iSTARR = c(15862, iSTARR_oldSTARR_43TF_pvalue[i, "old_STARR"]), 
                                                  control = c(15862, iSTARR_oldSTARR_43TF_pvalue[i, "old_STARR_control"])))
  
  iSTARR_oldSTARR_43TF_pvalue[i, "iSTARR_pvalue"] <- iSTARR_fisher_res_i$p.value
  iSTARR_oldSTARR_43TF_pvalue[i, "oldSTARR_pvalue"] <- oldSTARR_fisher_res_i$p.value
  
}


iSTARR_oldSTARR_43TF_df_enrich <- iSTARR_oldSTARR_43TF_pvalue[, c("TFs", "iSTARR_enrich", "old_STARR_enrich")]
iSTARR_oldSTARR_43TF_df_enrich <- melt(iSTARR_oldSTARR_43TF_df_enrich, id = c("TFs"))
colnames(iSTARR_oldSTARR_43TF_df_enrich) <- c("TFs", "enrich_type", "enrich_value")
iSTARR_oldSTARR_43TF_df_pvalue <- iSTARR_oldSTARR_43TF_pvalue[, c("TFs", "iSTARR_pvalue", "oldSTARR_pvalue")]
iSTARR_oldSTARR_43TF_df_pvalue <- melt(iSTARR_oldSTARR_43TF_df_pvalue, id = c("TFs"))
colnames(iSTARR_oldSTARR_43TF_df_pvalue) <- c("TFs1", "pvalue_type", "pvalue")

iSTARR_oldSTARR_43TF_df_enrich_pvalue <- cbind(iSTARR_oldSTARR_43TF_df_enrich, iSTARR_oldSTARR_43TF_df_pvalue, stringsAsFactors = F)

iSTARR_oldSTARR_43TF_df_enrich_pvalue$TFs <- factor(iSTARR_oldSTARR_43TF_df_enrich_pvalue$TFs, levels = unique(iSTARR_oldSTARR_43TF_df_enrich_pvalue$TFs))


point_size <- iSTARR_oldSTARR_43TF_df_enrich_pvalue$pvalue
point_size[point_size > 0.05] <- 7
point_size[point_size < 0.05 & point_size > 0.01] <- 8
point_size[point_size < 0.01 & point_size > 1e-10] <- 9
point_size[point_size < 1e-10] <- 10


#------------plot lollipop------------

pdf("iSTARR_oldSTARR_43TF_enrich_pvalue_lollipop.pdf", width = 14, height = 7)

ggplot(data = iSTARR_oldSTARR_43TF_df_enrich_pvalue, mapping = aes(x = TFs, y = enrich_value, color = enrich_type, size = as.factor(point_size)))+
  geom_point(alpha = 0.8)+
  scale_color_manual(values = c("blue", "orange"),
                     labels = c("iSTARR-seq", "STARR-seq"))+
  scale_size_manual(values = c(3, 5, 7, 8), 
                    labels = c(">0.05", "0.01~0.05", "0.01~1e-10", "<1e-10"))+
  geom_segment(aes(x = TFs, y = enrich_value, xend = TFs, yend = 0), 
               color = "gray30", size = 0.5)+
  geom_hline(yintercept = 0, color = "black")+
  theme_classic()+
  labs(title = "43 TFBS enrichment in STARR-seq enhancer", 
       x = element_blank(), 
       y = "TFBS enrichment to control",
       color = "", size = "P-value")+
  guides(color = guide_legend(override.aes = list(size=5), order = 1),
         size = guide_legend(order = 2))+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text=element_text(size=15, color = "black"),
        axis.text.x = element_text(angle = 50, vjust = 0.95, hjust = 1),
        axis.title.y=element_text(size=15),
        legend.position = "top",
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15))

dev.off()


#------------------TF enrichment in different enhancer cluster-----------------

iSTARR_cluster_43TF <- read.table("iSTARR_4cluster_43TFs_enrichment.txt", header = T, stringsAsFactors = F)

iSTARR_cluster_43TF_pvalue <- data.frame(iSTARR_cluster_43TF, 
                                         cluster1_pvalue = rep(0, nrow(iSTARR_cluster_43TF)),
                                         cluster2_pvalue = rep(0, nrow(iSTARR_cluster_43TF)),
                                         cluster3_pvalue = rep(0, nrow(iSTARR_cluster_43TF)),
                                         cluster4_pvalue = rep(0, nrow(iSTARR_cluster_43TF)),
                                         stringsAsFactors = F)

for (i in 1:nrow(iSTARR_cluster_43TF)) {
  cluster1_fisher_res_i <- fisher.test(data.frame(cluster1 = c(1643, iSTARR_cluster_43TF_pvalue[i, "cluster1_TF_count"]),
                                                  cluster1_control = c(1643, iSTARR_cluster_43TF_pvalue[i, "cluster1_control_TF_count"])))
  cluster2_fisher_res_i <- fisher.test(data.frame(cluster2 = c(873, iSTARR_cluster_43TF_pvalue[i, "cluster2_TF_count"]),
                                                  cluster2_control = c(873, iSTARR_cluster_43TF_pvalue[i, "cluster2_control_TF_count"])))
  cluster3_fisher_res_i <- fisher.test(data.frame(cluster3 = c(525, iSTARR_cluster_43TF_pvalue[i, "cluster3_TF_count"]),
                                                  cluster3_control = c(525, iSTARR_cluster_43TF_pvalue[i, "cluster3_control_TF_count"])))
  cluster4_fisher_res_i <- fisher.test(data.frame(cluster4 = c(1915, iSTARR_cluster_43TF_pvalue[i, "cluster4_TF_count"]),
                                                  cluster4_control = c(1915, iSTARR_cluster_43TF_pvalue[i, "cluster4_control_TF_count"])))
  
  iSTARR_cluster_43TF_pvalue[i, "cluster1_pvalue"] <- cluster1_fisher_res_i$p.value
  iSTARR_cluster_43TF_pvalue[i, "cluster2_pvalue"] <- cluster2_fisher_res_i$p.value
  iSTARR_cluster_43TF_pvalue[i, "cluster3_pvalue"] <- cluster3_fisher_res_i$p.value
  iSTARR_cluster_43TF_pvalue[i, "cluster4_pvalue"] <- cluster4_fisher_res_i$p.value
}

iSTARR_cluster_43TF_pvalue$TFs <- factor(iSTARR_cluster_43TF_pvalue$TFs, levels = iSTARR_cluster_43TF_pvalue$TFs)

point_size <- iSTARR_cluster_43TF_pvalue$cluster1_pvalue
point_size[point_size > 0.05] <- 6
point_size[point_size < 0.05 & point_size > 0.01] <- 7
point_size[point_size < 0.01 & point_size > 1e-10] <- 8
point_size[point_size < 1e-10] <- 9

# cluster1, cluster2, cluster3, cluster4

pdf("iSTARR_cluster1_43TF_enrich_pvalue_lollipop.pdf", width = 14, height = 7)

ggplot(data = iSTARR_cluster_43TF_pvalue, mapping = aes(x = TFs, y = cluster1_TF_enrich, size = as.factor(point_size)))+
  # "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"
  geom_point(color = "#E41A1C")+
  scale_size_manual(values = c(3, 5, 7, 8), 
                    labels = c(">0.05", "0.01~0.05", "0.01~1e-10", "<1e-10"))+
  geom_segment(aes(x = TFs, y = cluster1_TF_enrich, xend = TFs, yend = 0), 
               color = "gray30", size = 0.5)+
  geom_hline(yintercept = 0, color = "black")+
  ylim(-2.5, 2.5)+
  theme_classic()+
  labs(title = "43 TFBS enrichment in iSTARR-seq cluster1 enhancer", 
       x = element_blank(), 
       y = "TFBS enrichment to control",
       size = "P-value")+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text=element_text(size=15, color = "black"),
        axis.text.x = element_text(angle = 50, vjust = 0.95, hjust = 1),
        axis.title.y=element_text(size=15),
        legend.position = "top",
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15))

dev.off()




