####CIHAN OGUZ, NIAID Collaborative Bioinformatics Resource
####This script generates four figure panels for each of the three patient groups.


packages <- c("immunarch","ggpubr","tidyverse", "readxl", "lubridate", "corrr", "car","psych","gdata","writexl","psych","gdata","reshape2","reshape")

invisible(lapply(packages, library, character.only = TRUE))

load("Figure1_group1A_panelsAtoD.RData")
load("Figure2_group2A_panelsAtoD.RData")
load("Figure3_group1P_panelsAtoD.RData")


Figure1_panelA_group1A$Patient<-Figure1_panelA_group1A$Sample
Figure1_panelA_group1A$Patient<-str_extract(Figure1_panelA_group1A$Patient, regex("(?<=_)(.+)(?=_)"))
Figure1_panelA_group1A$Patient<-gsub("_.*","",Figure1_panelA_group1A$Patient)

Figure1_panelA_group1A<-Figure1_panelA_group1A[order(Figure1_panelA_group1A$Time, Figure1_panelA_group1A$Patient), ]


Figure1_panelB_group1A$Patient<-Figure1_panelB_group1A$Sample
Figure1_panelB_group1A$Patient<-str_extract(Figure1_panelB_group1A$Patient, regex("(?<=_)(.+)(?=_)"))
Figure1_panelB_group1A$Patient<-gsub("_.*","",Figure1_panelB_group1A$Patient)

Figure1_panelB_group1A<-Figure1_panelB_group1A[order(Figure1_panelB_group1A$Time, Figure1_panelB_group1A$Patient), ]



##########################
##########################

Figure2_panelA_group2A$Patient<-Figure2_panelA_group2A$Sample
Figure2_panelA_group2A$Patient<-str_extract(Figure2_panelA_group2A$Patient, regex("(?<=_)(.+)(?=_)"))
Figure2_panelA_group2A$Patient<-gsub("_.*","",Figure2_panelA_group2A$Patient)

Figure2_panelA_group2A<-Figure2_panelA_group2A[order(Figure2_panelA_group2A$Time, Figure2_panelA_group2A$Patient), ]


Figure2_panelB_group2A$Patient<-Figure2_panelB_group2A$Sample
Figure2_panelB_group2A$Patient<-str_extract(Figure2_panelB_group2A$Patient, regex("(?<=_)(.+)(?=_)"))
Figure2_panelB_group2A$Patient<-gsub("_.*","",Figure2_panelB_group2A$Patient)

Figure2_panelB_group2A<-Figure2_panelB_group2A[order(Figure2_panelB_group2A$Time, Figure2_panelB_group2A$Patient), ]




my_comparisons <- list( c("Wk0", "Wk24"))

set.seed(123)

Figure1_panelA<-ggviolin(Figure1_panelA_group1A[which(Figure1_panelA_group1A$Cell_type=="CD8" & Figure1_panelA_group1A$Patient_group!="1P" & Figure1_panelA_group1A$Patient_group!="2A"),], x = "Time", y = "breadth",color = "Time",shape = "Time", palette = c("#0000AC","#00AFBB","#E7B800"), ylab="",xlab="",title="",font.title = list(size = 14, color = "black"), add = "mean_sd", add.params = list(width = 0.8, size = 0.3), order=levels(Figure1_panelA_group1A$Time),ylim = c(10^-3.9, 10^-2.2),jitter = 0.3) + yscale("log10", .format = TRUE)+rotate_x_text(0)+ grids(linetype = "dashed")+stat_compare_means(comparisons = my_comparisons, size = 6,vjust = -1,tip.length=0.2,label.y = -2.5,paired = TRUE)+
  theme(axis.ticks.length.y=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

Figure1_panelB<-ggviolin(Figure1_panelA_group1A[which(Figure1_panelA_group1A$Cell_type=="CD8" & Figure1_panelA_group1A$Patient_group!="1P" & Figure1_panelA_group1A$Patient_group!="2A"),], x = "Time", y = "depth",color = "Time",shape = "Time",, palette = c("#0000AC","#00AFBB","#E7B800"), ylab="",xlab="",title="",font.title = list(size = 14, color = "black"), add = "mean_sd", add.params = list(width = 0.8, size = 0.3), order=levels(Figure1_panelA_group1A$Time),ylim = c(10^1.2, 10^3),jitter = 0.3) + yscale("log10", .format = TRUE)+rotate_x_text(0)+ grids(linetype = "dashed")+stat_compare_means(comparisons = my_comparisons, size = 6,tip.length=0.2,label.y = 2.6,vjust = -1,paired = TRUE)+
  theme(axis.ticks.length.y=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

Figure1_panelA<-ggpar(Figure1_panelA,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16,xtickslab = FALSE)+rremove("x.text")+rremove("legend")
Figure1_panelB<-ggpar(Figure1_panelB,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16,xtickslab = FALSE)+rremove("x.text")+rremove("legend")

Figure1_panelC<-Figure1_panelC_group1A+rremove("legend")
Figure1_panelD<-Figure1_panelD_group1A+rremove("legend")
Figure1_panelC<-ggpar(Figure1_panelC,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16)+rremove("legend")
Figure1_panelD<-ggpar(Figure1_panelD,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16)+rremove("legend")

pdf(file = "HIV_TCR_group_1A_figure_with_panels_oct15.pdf",width=9, height=6)
ggarrange(Figure1_panelA ,Figure1_panelB ,Figure1_panelC ,Figure1_panelD , ncol = 2, nrow = 2,font.label=list(size = 10), widths=c(1, 1, 1, 1, 1)) #, hjust=(-0.1)
dev.off()


set.seed(123)

Figure2_panelA<-ggviolin(Figure2_panelA_group2A[which(Figure2_panelA_group2A$Cell_type=="CD8" & Figure2_panelA_group2A$Patient_group!="1P" & Figure2_panelA_group2A$Patient_group!="1A"),], x = "Time", y = "breadth",color = "Time",shape = "Time", palette = c("#0000AC","#00AFBB","#E7B800"), ylab="",xlab="",title="",font.title = list(size = 14, color = "black"), add = "mean_sd", add.params = list(width = 0.8, size = 0.3), order=levels(Figure2_panelA_group2A$Time),ylim = c(10^-3.9, 10^-2.2),jitter = 0.3) + yscale("log10", .format = TRUE)+rotate_x_text(0)+ grids(linetype = "dashed")+stat_compare_means(comparisons = my_comparisons, size = 6,vjust = -1,tip.length=0.2,label.y = -2.5,paired = TRUE)+
  theme(axis.ticks.length.y=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

Figure2_panelB<-ggviolin(Figure2_panelA_group2A[which(Figure2_panelA_group2A$Cell_type=="CD8" & Figure2_panelA_group2A$Patient_group!="1P" & Figure2_panelA_group2A$Patient_group!="1A"),], x = "Time", y = "depth",color = "Time",shape = "Time",, palette = c("#0000AC","#00AFBB","#E7B800"), ylab="",xlab="",title="",font.title = list(size = 14, color = "black"), add = "mean_sd", add.params = list(width = 0.8, size = 0.3), order=levels(Figure2_panelA_group2A$Time),ylim = c(10^1.2, 10^3),jitter = 0.3) + yscale("log10", .format = TRUE)+rotate_x_text(0)+ grids(linetype = "dashed")+stat_compare_means(comparisons = my_comparisons, size = 6,tip.length=0.2,label.y = 2.6,vjust = -1,paired = TRUE)+
  theme(axis.ticks.length.y=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

Figure2_panelA<-ggpar(Figure2_panelA,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16,xtickslab = FALSE)+rremove("x.text")+rremove("legend")
Figure2_panelB<-ggpar(Figure2_panelB,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16,xtickslab = FALSE)+rremove("x.text")+rremove("legend")

Figure2_panelC<-Figure2_panelC_group2A+rremove("legend")
Figure2_panelD<-Figure2_panelD_group2A+rremove("legend")
Figure2_panelC<-ggpar(Figure2_panelC,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16)+rremove("legend")
Figure2_panelD<-ggpar(Figure2_panelD,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16)+rremove("legend")

pdf(file = "HIV_TCR_group_2A_figure_with_panels_oct15.pdf",width=9, height=6)
ggarrange(Figure2_panelA ,Figure2_panelB ,Figure2_panelC ,Figure2_panelD , ncol = 2, nrow = 2,font.label=list(size = 10), widths=c(1, 1, 1, 1, 1)) #, hjust=(-0.1)
dev.off()


set.seed(123)

Figure3_panelA<-ggviolin(Figure3_panelA_group1P[which(Figure3_panelA_group1P$Cell_type=="CD8" & Figure3_panelA_group1P$Patient_group!="2A" & Figure3_panelA_group1P$Patient_group!="1A"),], x = "Time", y = "breadth",color = "Time",shape = "Time", palette = c("#0000AC","#00AFBB","#E7B800"), ylab="",xlab="",title="",font.title = list(size = 14, color = "black"), add = "mean_sd", add.params = list(width = 0.8, size = 0.3), order=levels(Figure3_panelA_group1P$Time),ylim = c(10^-3.9, 10^-2.2),jitter = 0.3) + yscale("log10", .format = TRUE)+rotate_x_text(0)+ grids(linetype = "dashed")+stat_compare_means(comparisons = my_comparisons, size = 6,vjust = -1,tip.length=0.2,label.y = -2.5)+
  theme(axis.ticks.length.y=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

Figure3_panelB<-ggviolin(Figure3_panelA_group1P[which(Figure3_panelA_group1P$Cell_type=="CD8" & Figure3_panelA_group1P$Patient_group!="2A" & Figure3_panelA_group1P$Patient_group!="1A"),], x = "Time", y = "depth",color = "Time",shape = "Time",, palette = c("#0000AC","#00AFBB","#E7B800"), ylab="",xlab="",title="",font.title = list(size = 14, color = "black"), add = "mean_sd", add.params = list(width = 0.8, size = 0.3), order=levels(Figure3_panelA_group1P$Time),ylim = c(10^1.2, 10^3),jitter = 0.3) + yscale("log10", .format = TRUE)+rotate_x_text(0)+ grids(linetype = "dashed")+stat_compare_means(comparisons = my_comparisons, size = 6,tip.length=0.2,label.y = 2.6,vjust = -1)+
  theme(axis.ticks.length.y=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

Figure3_panelA<-ggpar(Figure3_panelA,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16,xtickslab = FALSE)+rremove("x.text")+rremove("legend")
Figure3_panelB<-ggpar(Figure3_panelB,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16,xtickslab = FALSE)+rremove("x.text")+rremove("legend")

Figure3_panelC<-Figure3_panelC_group1P+rremove("legend")
Figure3_panelD<-Figure3_panelD_group1P+rremove("legend")
Figure3_panelC<-ggpar(Figure3_panelC,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16)+rremove("legend")
Figure3_panelD<-ggpar(Figure3_panelD,font.y=16,font.x=16,font.legend = c(16),font.tickslab=c(16),font.title =16)+rremove("legend")


pdf(file = "HIV_TCR_group_1P_figure_with_panels_oct15.pdf",width=9, height=6)
ggarrange(Figure3_panelA ,Figure3_panelB ,Figure3_panelC ,Figure3_panelD , ncol = 2, nrow = 2,font.label=list(size = 10), widths=c(1, 1, 1, 1, 1))
dev.off()
