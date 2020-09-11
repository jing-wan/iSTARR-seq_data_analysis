#!/usr/bin/Rscript

require(pheatmap)
require(data.table)
require(RColorBrewer)


H3K4me3_in_enh <- fread("H3K4me3_in_enhancer.bed", header = F)
H3K9ac_in_enh <- fread("H3K9ac_in_enhancer.bed", header = F)
RNAPOII_in_enh <- fread("RNAPOII_in_enhancer.bed", header = F)
H2AZ_in_enh <- fread("H2AZ_in_enhancer.bed", header = F)
ATAC_in_enh <- fread("ATAC-seq_in_enhancer.bed", header = F)
DHS_in_enh <- fread("DNase-seq_in_enhancer.bed", header = F)
H3K4me1_in_enh <- fread("H3K4me1_in_enhancer.bed", header = F)
H3K27ac_in_enh <- fread("H3K27ac_in_enhancer.bed", header = F)
H3K27me3_in_enh <- fread("H3K27me3_in_enhancer.bed", header = F)
H3K36me3_in_enh <- fread("H3K36me3_in_enhancer.bed", header = F)
H2A_in_enh <- fread("H2A_in_enhancer.bed", header = F)
H3K36ac_in_enh <- fread("H3K36ac_in_enhancer.bed", header = F)
H3K56ac_in_enh <- fread("H3K56ac_in_enhancer.bed", header = F)
H3K9me2_in_enh <- fread("H3K9me2_in_enhancer.bed", header = F)
methy_in_enh <- fread("methy_in_enhancer.bed", header = F)


enhancer_epi_ma <- data.frame(DHS_in_enh$V11, ATAC_in_enh$V11, H3K4me3_in_enh$V11, RNAPOII_in_enh$V11, 
                              H3K27ac_in_enh$V11, H3K9ac_in_enh$V11, H3K36ac_in_enh$V6, H3K56ac_in_enh$V6, 
                              H3K36me3_in_enh$V11, H2AZ_in_enh$V11, H2A_in_enh$V6, H3K4me1_in_enh$V11, 
                              H3K27me3_in_enh$V11, H3K9me2_in_enh$V6, methy_in_enh$V4)

colnames(enhancer_epi_ma) <- c("DHS", "ATAC", "H3K4me3", "RNAPII", 
                               "H3K27ac", "H3K9ac", "H3K36ac", "H3K56ac", 
                               "H3K36me3", "H2AZ", "H2A", "H3K4me1", 
                               "H3K27me3", "H3K9me2", "methy")

enhancer_epi_ma <- enhancer_epi_ma+1
enhancer_epi_ma <- log(enhancer_epi_ma)

#extreme value processing
for (i in 1:15) {
  enhancer_epi_ma[enhancer_epi_ma[, i] > quantile(enhancer_epi_ma[, i], 0.99), i] <- quantile(enhancer_epi_ma[, i], 0.99)
}

#scale
enhancer_epi_ma <- scale(enhancer_epi_ma, center = T, scale = T)
enhancer_epi_ma <- as.data.frame(enhancer_epi_ma)

enhancer_epi_scale <- cbind(enhancer_epi_ma, chr = H3K4me3_in_enh$V1, start = H3K4me3_in_enh$V2, end = H3K4me3_in_enh$V3)

#kmeans cluster
#set initial cluster center
enhancer_kmeans <- kmeans(enhancer_epi_scale[, 1:15], centers = enhancer_epi_scale[1:4, 1:15], iter.max = 50)

#get cluster results
enhancer_epi_scale_cluster <- data.frame(enhancer_epi_scale, enhancer_kmeans$cluster)
#reorder
enhancer_epi_scale_cluster <- enhancer_epi_scale_cluster[order(enhancer_epi_scale_cluster$enhancer_kmeans.cluster), ]
rownames(enhancer_epi_scale_cluster) <- 1:nrow(enhancer_epi_scale_cluster)

write.table(enhancer_epi_scale_cluster[, 16:19], file = "iSTARR_enhancer_15epi_4cluster.bed", col.names = F, row.names = F, sep = "\t", quote = F)


#find optimal cluster results
wss <- (nrow(enhancer_epi_scale[, 1:15])-1)*sum(apply(enhancer_epi_scale[, 1:15], 2, var))

for (i in 2:15){
  wss[i] <- sum(kmeans(enhancer_epi_scale[, 1:15], centers = i, iter.max = 50)$withinss)
}

plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")


#plot heatmap
enh_num_per_cluster <- as.data.frame(table(enhancer_epi_scale_cluster$enhancer_kmeans.cluster))[, 2]
p_annotation_row = data.frame(enhancer = factor(rep(c("cluster1", "cluster2", "cluster3", "cluster4"), enh_num_per_cluster)))

#set color value range
bk <- c(-4, seq(-2, 2, length = 100), 4)

pheatmap(enhancer_epi_scale_cluster[, 1:15], 
         breaks = bk,
         color = c("blue", colorRampPalette(c("blue", "white", "red"))(100), "red"),
         cluster_rows = F, cluster_cols = F, show_rownames = F, 
         border_color = F, main = "Cluster enhancer by epigenetic signals", 
         annotation_row = p_annotation_row, gaps_row = cumsum(enh_num_per_cluster))


bk1 <- c(-2, seq(-1.5, 1.5, length = 100), 2.5)

pheatmap(enhancer_epi_scale_cluster[, 1:15], 
         breaks = bk1, border_color = NA,
         color = c("blue", colorRampPalette(c("blue", "white", "red"))(100), "red"),
         kmeans_k = enhancer_epi_scale[1:4, 1:15], 
         cluster_rows = F, cluster_cols = F, display_numbers = T, number_color = "black", fontsize_number = 10, 
         main = "Cluster enhancer by epigenetic signals", cex = 1.1)






