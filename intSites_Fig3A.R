library(GenomicRanges)
library(data.table)
library(rtracklayer)
library(gprofiler2)
library(ggplot2)

theme_set(theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")))

##### Plot out all interaction sites

int.data <- read.table("~/Parrott_Lab/AnoleAge/Misc/DSSdata/LogInteractSites_smoothed200.bed",
                       sep = "\t", header = FALSE)

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F")))

int.sites <- GRanges(seqnames = int.data$V1,
                     ranges = IRanges(start = int.data$V2, end = int.data$V3), strand = "*")

meth_file <- as.data.frame(fread(file = "E:/anole/CpG/processed/labeled_allSites_stranded.tab", header = TRUE, sep = "\t"))

meth.gr <- GRanges(seqnames = meth_file$chr, 
                   ranges = IRanges(start = meth_file$start, end = meth_file$end),
                   strand = "*", mcols = meth_file[,c(5:41)] )

overlapped <- subsetByOverlaps(meth.gr, int.sites, type = "end", ignore.strand = TRUE)

rm(meth.gr, meth_file)
gc()

overlap.frame <- as.data.frame(overlapped)

sites <- paste0(overlap.frame$seqnames, ":", overlap.frame$start)

formatted_meth <- as.data.frame(t(overlap.frame[,c(6:42)]))
colnames(formatted_meth) <- sites

missing_sites <- colSums(apply(formatted_meth, 2, is.na))

annot_meth <- data.frame(annotations, formatted_meth[,-which(missing_sites > 30)])

annot_long <- tidyr::pivot_longer(annot_meth, cols = c(4:length(annot_meth)), names_to = "CpG", values_to = "Methylation")


a <- ggplot(annot_long, aes(x = ages, y = Methylation, color = sex)) + 
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    geom_point() + facet_wrap(~CpG, ncol = 4) + ylim(0,100) +
    xlab("Age (months)") + labs(color = "Sex") +
    theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))
a

ggsave(a, filename = "~/Parrott_Lab/AnoleAge/Misc/chr3/IntSites_plotted.svg", device = "svg",
       width = 7, height = 4)
