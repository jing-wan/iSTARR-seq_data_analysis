#!/usr/bin/Rscript

require(reshape2)
require(ggplot2)

#================================old STARR-seq=========================================

old_STARR_up_TSS <- read.table("old_STARR_aggregate_in_TSS_up_20win.bed", header = F)

old_STARR_in_TSS <- old_STARR_up_TSS[1:20, 2] + rev(old_STARR_up_TSS[21:40, 2])

old_STARR_genebody <- read.table("old_STARR_aggregate_in_gene_body_40win.bed", header = F)
old_STARR_in_genebody <- old_STARR_genebody[1:40, 2] + rev(old_STARR_genebody[41:80, 2])

old_STARR_down_TTS <- read.table("old_STARR_aggregate_in_TTS_down_20win.bed", header = F)
old_STARR_in_TTS <- old_STARR_down_TTS[1:20, 2] + rev(old_STARR_down_TTS[21:40, 2])

old_STARR_TSS_genebody_TTS <- data.frame(position = c(1:80), 
                                         density = c(old_STARR_in_TSS, old_STARR_in_genebody, old_STARR_in_TTS))

old_STARR_TSS_genebody_TTS <- data.frame(old_STARR_TSS_genebody_TTS, 
                                         percent = old_STARR_TSS_genebody_TTS$density/sum(old_STARR_TSS_genebody_TTS$density))


#---------------------------random control-----------------------------------

random_up_TSS <- read.table("random_old_STARR_aggregate_in_TSS_up_20win.bed", header = F)

random_in_TSS <- random_up_TSS[1:20, 2] + rev(random_up_TSS[21:40, 2])

random_genebody <- read.table("random_old_STARR_aggregate_in_gene_body_40win.bed", header = F)
random_in_genebody <- random_genebody[1:40, 2] + rev(random_genebody[41:80, 2])

random_down_TTS <- read.table("random_old_STARR_aggregate_in_TTS_down_20win.bed", header = F)
random_in_TTS <- random_down_TTS[1:20, 2] + rev(random_down_TTS[21:40, 2])

random_TSS_genebody_TTS <- data.frame(position = c(1:80), 
                                      density = c(random_in_TSS, random_in_genebody, random_in_TTS))

old_random_TSS_genebody_TTS <- data.frame(random_TSS_genebody_TTS,
                                          percent = random_TSS_genebody_TTS$density/sum(random_TSS_genebody_TTS$density))

#================================iSTARR-seq=========================================

STARR_up_TSS <- read.table("STARR-seq_enhancer_aggregate_in_TSS_up_20win.bed", header = F)
#把正负链反向叠加
STARR_in_TSS <- STARR_up_TSS[1:20, 2] + rev(STARR_up_TSS[21:40, 2])

STARR_genebody <- read.table("STARR-seq_enhancer_aggregate_in_gene_body_40win.bed", header = F)
STARR_in_genebody <- STARR_genebody[1:40, 2] + rev(STARR_genebody[41:80, 2])

STARR_down_TTS <- read.table("STARR-seq_enhancer_aggregate_in_TTS_down_20win.bed", header = F)
STARR_in_TTS <- STARR_down_TTS[1:20, 2] + rev(STARR_down_TTS[21:40, 2])

STARR_TSS_genebody_TTS <- data.frame(position = c(1:80), 
                                     density = c(STARR_in_TSS, STARR_in_genebody, STARR_in_TTS))

iSTARR_TSS_genebody_TTS <- data.frame(STARR_TSS_genebody_TTS, 
                                      percent = STARR_TSS_genebody_TTS$density/sum(STARR_TSS_genebody_TTS$density))

#---------------------------random control-----------------------------------
random_up_TSS <- read.table("random_STARR_aggregate_in_TSS_up_20win.bed", header = F)

random_in_TSS <- random_up_TSS[1:20, 2] + rev(random_up_TSS[21:40, 2])

random_genebody <- read.table("random_STARR_aggregate_in_gene_body_40win.bed", header = F)
random_in_genebody <- random_genebody[1:40, 2] + rev(random_genebody[41:80, 2])

random_down_TTS <- read.table("random_STARR_aggregate_in_TTS_down_20win.bed", header = F)
random_in_TTS <- random_down_TTS[1:20, 2] + rev(random_down_TTS[21:40, 2])

random_TSS_genebody_TTS <- data.frame(position = c(1:80), 
                                      density = c(random_in_TSS, random_in_genebody, random_in_TTS))

iSTARR_random_TSS_genebody_TTS <- data.frame(random_TSS_genebody_TTS,
                                             percent = random_TSS_genebody_TTS$density/sum(random_TSS_genebody_TTS$density))


#============================merge STARR and random=================================

old_iSTARR_random_TSS_genebody_TTS <- data.frame(position = c(1:80), 
                                                 iSTARR = iSTARR_TSS_genebody_TTS$percent,
                                                 iSTARR_random = iSTARR_random_TSS_genebody_TTS$percent,
                                                 STARR = old_STARR_TSS_genebody_TTS$percent, 
                                                 STARR_random = old_random_TSS_genebody_TTS$percent)

old_iSTARR_random_TSS_genebody_TTS <- melt(old_iSTARR_random_TSS_genebody_TTS, id = "position")
colnames(old_iSTARR_random_TSS_genebody_TTS) <- c("position", "type", "density")


pdf("old_iSTARR_distribution_in_TSS_genebody_TTS.pdf", height = 6, width = 12)

ggplot(old_iSTARR_random_TSS_genebody_TTS, aes(x = position, y = density, group = type, color = type, linetype = type))+
  geom_line(size = 1)+
  labs(title = "STARR-seq enhancer distribution in TSS-genebody-TTS", x = "position", y = "Enhancer density")+
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50, 60, 70, 80), 
                     labels = c("2k", "1k", "TSS", "25%", "50%", "75%", "TTS", "1k", "2k"))+
  scale_colour_manual(values=c("darkblue", "gray60", "darkred", "gray"))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "dashed"))+
  geom_vline(xintercept = c(20, 60), linetype = 2, color = "gray")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text=element_text(size=14, color = "black"),
        axis.title.x=element_text(size=16, color = "black"),
        axis.title.y=element_text(size=16, color = "black"), 
        legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14))

dev.off()






