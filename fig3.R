library(GenomicRanges)
library(data.table)
library(rtracklayer)
library(gprofiler2)
library(ggplot2)

theme_set(theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")))

## Import CG counts per chromosome - generated using fastaRegexFinder.py from bioinformatics-cafe on github

table.cg.sorted <- read.table("C:/Users/sheal/Desktop/Parrott_Lab/Anole_Data/AnoSag2.1/CG_perChr.txt")

# read in sex DMC data

sex.data <- read.table("C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/data/SexSites_smoothed200.bed",
                       sep = "\t", header = FALSE)

# tabulate and test enrichment of sex-related sites on chr3

sex.tab <- as.data.frame(table(sex.data$V1))

sex.cont <- data.frame("SexSites" = c(sex.tab$Freq[6], sum(sex.tab$Freq[c(1:5, 8, 9, 10, 12, 13)])),
                       "NotSexSites" = c(table.cg.sorted$Freq[7], sum(table.cg.sorted$Freq[c(1:5, 8:14)])))
rownames(sex.cont) <- c("Chr3", "OtherAutosomes")

chisq.test(sex.cont) ### Enrichment is significant, p < 2.2e-16


sex.sites <- GRanges(seqnames = sex.data$V1,
                     ranges = IRanges(start = sex.data$V2, end = sex.data$V3), strand = "*")

## Now perform the same test on age-related sites

age.data <- read.table("C:/Users/sheal/Desktop/AnoleAge/Fig1.1-DSSLogAge/data/LogAgeSites_smoothed200.bed",
                       sep = "\t", header = FALSE)

age.tab <- as.data.frame(table(age.data$V1))

age.cont <- data.frame("AgeSites" = c(age.tab$Freq[10], sum(age.tab$Freq[c(1:6,8,11:15,17)])),
                       "NotAgeSites" = c(table.cg.sorted$Freq[7], sum(table.cg.sorted$Freq[c(1:5, 8:14)])))
rownames(age.cont) <- c("Chr3", "OtherAutosomes")

chisq.test(age.cont) ### Enrichment is significant, p < 2.2e-16
chisq.test(age.cont)$expected

age.sites <- GRanges(seqnames = age.data$V1,
                     ranges = IRanges(start = age.data$V2, end = age.data$V3), strand = "*")

# get CpGs which are both age- and sex-related (overlaps cpgs)

overlaps <- subsetByOverlaps(age.sites, sex.sites)

## Create chromomap of this region

library(chromoMap)

annot.df <- data.frame("name" = paste0(seqnames(overlaps), ".", start(overlaps)),
                    "chr" = seqnames(overlaps),
                    "start" = start(overlaps),
                    "end" = end(overlaps),
                    "strand" = rep(1, length(overlaps)))

annot.df$chr <- stringr::str_replace_all(annot.df$chr, "scaffold_", "chr")

write.table(annot.df, "C:/Users/sheal/Desktop/AnoleAge/Misc/chr3/overlapSites_annot.txt", 
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

annot <- "C:/Users/sheal/Desktop/AnoleAge/Misc/chr3/overlapSites_annot.txt"
chr_file <- "C:/Users/sheal/Desktop/AnoleAge/Fig1.1-DSSLogAge/data/AnoSag_chrList.txt"

c <- chromoMap(chr_file, annot, ch_gap = 2, chr_color = rep("grey", 15), anno_col = "blue")

c

# fisher test for enrichment of overlap sites on chr3

overlaps.tab <- as.data.frame(table(as.data.frame(overlaps)$seqnames))

overlaps.cont <- data.frame("OverlapSites" = c(overlaps.tab$Freq[10], sum(overlaps.tab$Freq[c(1:6, 8, 11:15, 17)])),
                       "NotOverlapSites" = c(table.cg.sorted$Freq[7], sum(table.cg.sorted$Freq[c(1:5, 8:14)])))

rownames(overlaps.cont) <- c("Chr3", "OtherAutosomes")


fisher.test(overlaps.cont) ### Enrichment is significant, p < 2.2e-16

## generate venn diagram

library(ggVennDiagram)

x <- list(sDMCs = paste0(seqnames(sex.sites), ":", start(sex.sites)),
          aDMCs = paste0(seqnames(age.sites), ":", start(age.sites)))

b <- ggVennDiagram(x, label_color = "black", label_size = 4, set_size = 5) + scale_color_manual(values = c("black", "black")) +
    scale_fill_steps(low = "skyblue", high = "navy") + theme(legend.position = "none")
b

ggsave(b, filename = "C:/Users/sheal/Desktop/AnoleAge/Misc/chr3/Sex_Age_venn.png",
       width = 3, height = 2)

## aggregate and plot out trajectory of chr3-specific overlap sites

overlaps <- overlaps[which(seqnames(overlaps) == "scaffold_3")]
overlaps.df <- as.data.frame(overlaps)

meth_file <- as.data.frame(fread(file = "E:/anole/CpG/processed/labeled_allSites_stranded.tab", header = TRUE, sep = "\t"))

meth.gr <- GRanges(seqnames = meth_file$chr, 
                   ranges = IRanges(start = meth_file$start, end = meth_file$end),
                   strand = "*", mcols = meth_file[,c(5:41)] )


overlapped <- subsetByOverlaps(meth.gr, overlaps, type = "end", ignore.strand = TRUE)
overlapped.df <- as.data.frame(overlapped)[,-c(1:5)]



region.labs <- paste0(seqnames(overlaps), ".", start(overlaps))


annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F")))


region.frame <- data.frame(annotations, t(overlapped.df))
names(region.labs) <- colnames(region.frame)[-c(1,2,3)]

region.avg <- data.frame(annotations, "RegionAvg" = rowMeans(t(overlapped.df), na.rm = TRUE))

region.long <- na.omit(tidyr::pivot_longer(region.frame, cols = c(4:50), values_to = "Meth", names_to = "Site"))

region.long$bp <- start(overlaps)[as.numeric(stringr::str_remove_all(region.long$Site, "X"))]

overlap.cor <- psych::corr.test(as.data.frame(t(overlapped.df)), annotations$ages, use = "pairwise")

siteDirs <- data.frame("site" = paste0("X", c(1:47)), 
           "dir" = ifelse(psych::corr.test(as.data.frame(t(overlapped.df)), annotations$ages, use = "pairwise")$r > 0,
                                   "Increasing", "Decreasing"))

region.long.dir <- merge(region.long, siteDirs, by.x = "Site", by.y = "site", sort = FALSE, all.x = TRUE)


d <- ggplot(region.long.dir, aes(x = ages, y = Meth, color = sex)) + 
    geom_jitter(color = "black", alpha = 0.3, height = 1, width = 1) + facet_grid(cols = vars(sex), rows = vars(dir)) +
    geom_smooth(method = "lm", formula = y ~ log(x)) +
    theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14),
          legend.position = "top", strip.text = element_text(size = 12)) + labs(x = "Age (months)", color = "Sex")
d

ggsave(d, filename = "C:/Users/sheal/Desktop/AnoleAge/Misc/chr3/chr3_regionMeth_byDir.png",
       width = 5, height = 5)



