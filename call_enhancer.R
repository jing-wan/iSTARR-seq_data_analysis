#!/usr/bin/Rscript

#====================call enhancer by BasicSTARRseq==================

require(BasicSTARRseq)

starrseq_data11 <- BasicSTARRseq::STARRseqData(sample="cDNA_sort.bam", control="plasmid_sort.bam", pairedEnd=TRUE)
peaks <- BasicSTARRseq::getPeaks(starrseq_data11, minQuantile = 0.9, peakWidth = 500, maxPval = 0.001,
                                   deduplicate = T, model = 1)

peaks_df <- data.frame(seqnames=seqnames(peaks), start=start(peaks), end=end(peaks),
                   width=end(peaks)-start(peaks)+1,
                   strand=strand(peaks),
                   sampleCov=elementMetadata(peaks)$sampleCov,
                   controlCov=elementMetadata(peaks)$controlCov,
                   pVal=elementMetadata(peaks)$pVal,
                   enrichment=elementMetadata(peaks)$enrichment)

peaks_df_order_pval <- peaks_df[order(peaks_df$pVal), ]
#FDR
qval <- p.adjust(peaks_df_order_pval$pVal, method = "fdr", n = length(peaks_df_order_pval$pVal))
peaks_df_fdr <- cbind(peaks_df_order_pval, FDR = qval)

peaks_df_fdr <- peaks_df_fdr[order(peaks_df_fdr$seqnames, peaks_df_fdr$start), ]

write.table(peaks_df_fdr, file="STARR-seq_peaks_fdr.bed", quote=F, sep="\t", row.names=F, col.names=F)



