library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggpubr)

#### CpG --------------------------------------------------------------------------------------------

theme_set(theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")))

## First, read in significant sex DMR sites and overlap to get methylation data for each individual here

sex.data <- read.table("C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/data/SexSites_smoothed200.bed",
                       sep = "\t", header = FALSE)

sex.sites <- GRanges(seqnames = sex.data$V1,
                     ranges = IRanges(start = sex.data$V2, end = sex.data$V3), strand = "*")


## Subset by Sites which are hypermethylated in females
f.data <- sex.data[which(sex.data$V4 < 0),]
f.sites <- GRanges(seqnames = f.data$V1,
                   ranges = IRanges(start = f.data$V2, end = f.data$V3), strand = "*")

## And Sites which are hypermethylated in males
m.data <- sex.data[which(sex.data$V4 > 0),]
m.sites <- GRanges(seqnames = m.data$V1,
                   ranges = IRanges(start = m.data$V2, end = m.data$V3), strand = "*")


## Get methylation matrix and sample annotations
meth_file <- as.data.frame(fread(file = "E:/anole/CpG/processed/labeled_allSites_stranded.tab", header = TRUE, sep = "\t"))


meth.gr <- GRanges(seqnames = meth_file$chr, 
                   ranges = IRanges(start = meth_file$start, end = meth_file$end),
                   strand = "*", mcols = meth_file[,c(5:41)] )

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F")))

## Table sex DMCs by chromosome

table <- as.data.table(table(sex.data$V1))
table.sorted <- table[stringr::str_order(table$V1, numeric = TRUE, decreasing = FALSE),] ## Table scaffold counts and order by number
table.sorted$V1 <- factor(table.sorted$V1, levels = stringr::str_sort(table.sorted$V1, numeric = TRUE, decreasing = FALSE))

chromSizes <- c(355360412,
                314158933,
                281461234,
                253587442,
                206893075,
                135963048,
                110948685,
                47139513,
                45588774,
                39276609,
                38356760,
                36225908)

## Calculate and plot site density per Mbp
table.sorted <- table.sorted[which(table.sorted$V1 %in% paste0("scaffold_", c(1:14))),]
table.sorted$size <- chromSizes
table.sorted$sizeMb <- table.sorted$size / 1000000
table.sorted$siteDensity <- table.sorted$N / table.sorted$sizeMb
table.sorted$chrName <- paste0("chr", c(1:6, "X", 8:12))
table.sorted$chrName <- factor(table.sorted$chrName, levels = table.sorted$chrName)


a <- ggplot(table.sorted, aes(x = siteDensity, y = chrName)) + geom_col(fill = "navy") + 
    ylab("") + xlab("Sex-Related Sites per Mbp") +
    ggbreak::scale_x_break(c(0.35,7.5), expand = FALSE) +
    theme(axis.title = element_text(size = 16), strip.background = element_blank())
a

ggsave(a, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/SexSites_byMBp.png",
       width = 12, height = 4)

## Now pull out only DMCs which are found on the X chromosome

xChrom.sites <- sex.sites[seqnames(sex.sites) == "scaffold_7"]
xChrom.data <- sex.data[which(sex.data$V1 == "scaffold_7"),]

## Effect size plots
f.xChrom.sites <- xChrom.data[which(xChrom.data$V4 < 0),]
f.xChrom.sites$Sex <- rep("F", length(f.xChrom.sites$V1))

m.xChrom.sites <- xChrom.data[which(xChrom.data$V4 > 0),]

m.xChrom.sites$Sex <- rep("M", length(m.xChrom.sites$V1))

total.xChrom.sites <- rbind(f.xChrom.sites, m.xChrom.sites)

xChrom.overlapped <- subsetByOverlaps(meth.gr, xChrom.sites, type = "end", ignore.strand = TRUE)

xChrom.overlap.frame <- as.data.frame(xChrom.overlapped)

xChrom.formatted_meth <- as.data.frame(t(xChrom.overlap.frame[,c(6:42)]))
xChrom.annot_meth <- cbind(annotations, xChrom.formatted_meth)


f.Xchrome <- subsetByOverlaps(xChrom.overlapped, f.sites, type = "any", maxgap = -1L, ignore.strand = TRUE)
m.Xchrome <- subsetByOverlaps(xChrom.overlapped, m.sites, type = "any", maxgap = -1L, ignore.strand = TRUE)

f_X.overlap.frame <- as.data.frame(f.Xchrome)
f_X.formatted_meth <- as.data.frame(t(f_X.overlap.frame[,c(6:42)]))
f_X.annot_meth <- cbind(annotations, f_X.formatted_meth)

f_X.sites_females <- f_X.annot_meth[which(f_X.annot_meth$sex == "F"),]
f_X.f_avg <- data.frame("means" = colMeans(f_X.sites_females[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("F", length(f.Xchrome)))

f_X.sites_males <- f_X.annot_meth[which(f_X.annot_meth$sex == "M"),]
f_X.m_avg <- data.frame("means" = colMeans(f_X.sites_males[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("M", length(f.Xchrome)))

f_X.sites_total <- as.data.frame(rbind(f_X.f_avg, f_X.m_avg))

d <- ggplot(data = f_X.sites_total, aes(x = sex, y = means, fill = sex)) + 
    geom_violin() + geom_jitter(color = "black", size = 0.8, alpha = 0.2, height = 0) +
    labs(title = "Female Hypermethylated", y = "Average Meth %", fill = "Sex") +
    theme(title = element_text(size = 16), axis.title = element_text(size = 16), 
          axis.title.x = element_blank(), legend.key.size = unit(0.3, "inches"), 
          legend.text = element_text(size = 16))
d

ggsave(d, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/xChrom_sexDMC_FemHyperMeth.png",
       width = 4, height = 4)


m_X.overlap.frame <- as.data.frame(m.Xchrome)
m_X.formatted_meth <- as.data.frame(t(m_X.overlap.frame[,c(6:42)]))
m_X.annot_meth <- cbind(annotations, m_X.formatted_meth)

m_X.sites_females <- m_X.annot_meth[which(m_X.annot_meth$sex == "F"),]
m_X.f_avg <- data.frame("means" = colMeans(m_X.sites_females[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("F", length(m.Xchrome)))

m_X.sites_males <- m_X.annot_meth[which(m_X.annot_meth$sex == "M"),]
m_X.m_avg <- data.frame("means" = colMeans(m_X.sites_males[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("M", length(m.Xchrome)))

m_X.sites_total <- as.data.frame(rbind(m_X.f_avg, m_X.m_avg))

d2 <- ggplot(data = m_X.sites_total, aes(x = sex, y = means, fill = sex)) + 
    geom_violin() + geom_jitter(color = "black", size = 0.8, alpha = 0.3, height = 0) +
    labs(title = "Male Hypermethylated", y = "Average Meth %", fill = "Sex") +
    theme(title = element_text(size = 16), axis.title = element_text(size = 16), 
          axis.title.x = element_blank(), legend.key.size = unit(0.3, "inches"), 
          legend.text = element_text(size = 16))

d2

ggsave(d2, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/xChrom_sexDMC_MaleHyperMeth.png",
       width = 4, height = 4)

d3 <- ggarrange(d, d2, nrow = 1, common.legend = TRUE)
d3

ggsave(d3, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/xChrom_sexDMC_COMBO.png",
       width = 8, height = 4.25)

### Get pseudo-effect sizes from the average %methylation differences of our two sexes:
m_X.sites_diff <- na.omit(data.frame("maleAvg" = m_X.m_avg$means, "femAvg" = m_X.f_avg$means, "diff" = ( m_X.m_avg$means - m_X.f_avg$means )))

m_X.sites_diff.hist <- ggplot(data = m_X.sites_diff, aes(x = diff)) + 
    geom_histogram(binwidth = 5, fill = "grey40", color = "black") + labs(x = "% Meth Difference") +
    theme(axis.title.y = element_text(size = 18))
m_X.sites_diff.box <- ggplot(data = m_X.sites_diff, aes(y = diff)) + 
    geom_boxplot(fill = "blue", color = "black") + labs(y = "% Meth Difference") +
    theme(axis.ticks.x = element_blank())

c <- ggarrange(m_X.sites_diff.hist, m_X.sites_diff.box + rremove("x.text"), ncol = 2, nrow = 1, align = "v",
          widths = c(1.5,.5))
c

ggsave(c, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/xChrom_Mhypermeth_methDiff.png",
       width = 5, height = 4)

f_X.sites_diff <- na.omit(data.frame("maleAvg" = f_X.m_avg$means, "femAvg" = f_X.f_avg$means, "diff" = ( f_X.f_avg$means - f_X.m_avg$means )))

f_X.sites_diff.hist <- ggplot(data = f_X.sites_diff, aes(x = diff)) + 
    geom_histogram(binwidth = 5, fill = "grey40", color = "black") + labs(x = "% Meth Difference") +
    theme(axis.title.y = element_text(size = 18))
f_X.sites_diff.box <- ggplot(data = f_X.sites_diff, aes(y = diff)) + 
    geom_boxplot(fill = "red", color = "black") + labs(y = "% Meth Difference") +
    theme(axis.ticks.x = element_blank())

c2 <- ggarrange(f_X.sites_diff.hist, f_X.sites_diff.box + rremove("x.text"), ncol = 2, nrow = 1, align = "v",
          widths = c(1.5,.5))
c2

ggsave(c2, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/xChrom_Fhypermeth_methDiff.png",
       width = 5, height = 4)

c3 <- ggarrange(c, c2, nrow = 2)
c3

ggsave(c3, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/xChrom_COMBO_methDiff.png",
       width = 5, height = 4)

mean(f_X.sites_diff$diff)

## Now pull out only DMCs which are *NOT* found on the X chromosome

autosom.sites <- sex.sites[seqnames(sex.sites) != "scaffold_7"]
autosom.data <- sex.data[which(sex.data$V1 != "scaffold_7"),]


autosom.overlapped <- subsetByOverlaps(meth.gr, autosom.sites, type = "any", maxgap = -1L, ignore.strand = TRUE)

autosom.overlap.frame <- as.data.frame(autosom.overlapped)

autosom.formatted_meth <- as.data.frame(t(autosom.overlap.frame[,c(6:42)]))
autosom.annot_meth <- cbind(annotations, autosom.formatted_meth)

m.a.avgs <- colMeans(autosom.annot_meth[which(autosom.annot_meth$sex == "M"),c(-1,-2,-3)], na.rm = TRUE)
f.a.avgs <- colMeans(autosom.annot_meth[which(autosom.annot_meth$sex == "F"),c(-1,-2,-3)], na.rm = TRUE)

autosom.avgs <- data.frame("avgs" = c(m.a.avgs, f.a.avgs), "sex" = c(rep("M", length(m.a.avgs)), rep("F", length(m.a.avgs))))

## Now split autosomal loci into hypermethylated in males, or hypermethylated in females

f.autosome <- subsetByOverlaps(autosom.overlapped, f.sites, type = "any", maxgap = -1L, ignore.strand = TRUE)
m.autosome <- subsetByOverlaps(autosom.overlapped, m.sites, type = "any", maxgap = -1L, ignore.strand = TRUE)

fAuto.overlap.frame <- as.data.frame(f.autosome)
fAuto.formatted_meth <- as.data.frame(t(fAuto.overlap.frame[,c(6:42)]))
fAuto.annot_meth <- cbind(annotations, fAuto.formatted_meth)

fAuto.sites_females <- fAuto.annot_meth[which(fAuto.annot_meth$sex == "F"),]
fAuto.f_avg <- data.frame("means" = colMeans(fAuto.sites_females[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("F", length(f.autosome)))

fAuto.sites_males <- fAuto.annot_meth[which(fAuto.annot_meth$sex == "M"),]
fAuto.m_avg <- data.frame("means" = colMeans(fAuto.sites_males[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("M", length(f.autosome)))

fAuto.sites_total <- as.data.frame(rbind(fAuto.f_avg, fAuto.m_avg))

b <- ggplot(data = fAuto.sites_total, aes(x = sex, y = means, fill = sex)) + geom_violin() + 
    geom_jitter(color = "black", size = 0.5, alpha = 0.2) +
        labs(title = "Female Hypermethylated", y = "Average Meth %", fill = "Sex") + 
    theme(title = element_text(size = 16), axis.title = element_text(size = 16), 
          axis.title.x = element_blank(), legend.key.size = unit(0.3, "inches"), 
          legend.text = element_text(size = 16))
b

ggsave(b, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/autosom_sexDMC_FemHyperMeth.png",
       width = 5, height = 4)



mAuto.overlap.frame <- as.data.frame(m.autosome)
mAuto.formatted_meth <- as.data.frame(t(mAuto.overlap.frame[,c(6:42)]))
mAuto.annot_meth <- cbind(annotations, mAuto.formatted_meth)

mAuto.sites_females <- mAuto.annot_meth[which(mAuto.annot_meth$sex == "F"),]
mAuto.f_avg <- data.frame("means" = colMeans(mAuto.sites_females[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("F", length(m.autosome)))

mAuto.sites_males <- mAuto.annot_meth[which(mAuto.annot_meth$sex == "M"),]
mAuto.m_avg <- data.frame("means" = colMeans(mAuto.sites_males[,c(-1,-2,-3)], na.rm = TRUE), "sex" = rep("M", length(m.autosome)))

mAuto.sites_total <- as.data.frame(rbind(mAuto.f_avg, mAuto.m_avg))

b2 <- ggplot(data = mAuto.sites_total, aes(x = sex, y = means, fill = sex)) + 
    geom_violin() + geom_jitter(color = "black", size = 0.5, alpha = 0.3) +
    labs(title = "Male Hypermethylated", y = "Average Meth %", fill = "Sex") + 
    theme(title = element_text(size = 16), axis.title = element_text(size = 16), 
          axis.title.x = element_blank(), legend.key.size = unit(0.3, "inches"), 
          legend.text = element_text(size = 16))
b2

ggsave(b2, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/autosom_sexDMC_MaleHyperMeth.png",
       width = 5, height = 4)

b3 <- ggarrange(b, b2, nrow = 1, common.legend = TRUE)
b3

ggsave(b3, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/autosom_sexDMC_COMBO.png",
       width = 8, height = 4.25)


### Get pseudo-effect sizes from the average %methylation differences of our two sexes:
mAuto.sites_diff <- na.omit(data.frame("maleAvg" = mAuto.m_avg$means, "femAvg" = mAuto.f_avg$means, "diff" = ( mAuto.m_avg$means - mAuto.f_avg$means )))

mAuto.sites_diff.hist <- ggplot(data = mAuto.sites_diff, aes(x = mAuto.sites_diff$diff)) + 
    geom_histogram(binwidth = 5, fill = "grey40", color = "black") + labs(x = "% Meth Difference") +
    theme(axis.title.y = element_text(size = 18))
mAuto.sites_diff.box <- ggplot(data = mAuto.sites_diff, aes(y = mAuto.sites_diff$diff)) + 
    geom_boxplot(fill = "blue", color = "black") + labs(y = "% Meth Difference") +
    theme(axis.ticks.x = element_blank())

e <- ggarrange(mAuto.sites_diff.hist, mAuto.sites_diff.box + rremove("x.text"), ncol = 2, nrow = 1, align = "v",
          widths = c(1.5,.5))
e

ggsave(e, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/Autosom_Mhypermeth_methDiff.png",
       width = 5, height = 4)

fAuto.sites_diff <- na.omit(data.frame("maleAvg" = fAuto.m_avg$means, "femAvg" = fAuto.f_avg$means, "diff" = ( fAuto.f_avg$means - fAuto.m_avg$means )))

fAuto.sites_diff.hist <- ggplot(data = fAuto.sites_diff, aes(x = fAuto.sites_diff$diff)) + 
    geom_histogram(binwidth = 5, fill = "grey40", color = "black") + labs(x = "% Meth Difference") +
    theme(axis.title.y = element_text(size = 18))
fAuto.sites_diff.box <- ggplot(data = fAuto.sites_diff, aes(y = fAuto.sites_diff$diff)) + 
    geom_boxplot(fill = "red", color = "black") + labs(y = "% Meth Difference") +
    theme(axis.ticks.x = element_blank())

e2 <- ggarrange(fAuto.sites_diff.hist, fAuto.sites_diff.box + rremove("x.text"), ncol = 2, nrow = 1, align = "v",
          widths = c(1.5,.5))
d2

ggsave(e2, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/Autosom_Fhypermeth_methDiff.png",
       width = 5, height = 4)

e3 <- ggarrange(e, e2, nrow = 2)
e3

ggsave(d3, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/autosom_COMBO_methDiff.png",
       width = 5, height = 4)

## Now do a similar analysis, but only looking at the cluster of DMCs / DMRs located on the ~550kb region of X chromosome

#scaffold_7:62,340,362-62,891,389

x.site <- GRanges(seqnames = "scaffold_7", ranges = IRanges(start = 62400000, end = 62725000), strand = "*")

## Plot chromoplot of this site's location
library(chromoMap)
chr_file <-  read.table("C:/Users/sheal/Desktop/AnoleAge/Fig1.1-DSSLogAge/data/AnoSag_chrList.txt")
annot_file <- data.frame("V1" = "FHM", "V2" = "chrX", "V3" = 62340362, "V4" = 62891389)

mapobj <- chromoMap(list(chr_file[7,]), list(annot_file),
          n_win.factor = 1, chr_color = "grey", anno_col = "yellow",
          interactivity = F, segment_annotation = TRUE)
mapobj


x.overlapped <- subsetByOverlaps(meth.gr, x.site, type = "any", maxgap = 0L, ignore.strand = TRUE)
x.overlapped <- subsetByOverlaps(x.overlapped, sex.sites, type = "any", maxgap = 0L, ignore.strand = TRUE)
x.overlap.frame <- as.data.frame(x.overlapped)

x.formatted_meth <- as.data.frame(t(x.overlap.frame[,c(6:42)]))
x.annot_meth <- cbind(annotations, x.formatted_meth)

m.avg <- rowMeans(x.annot_meth[which(x.annot_meth$sex == "M"),c(-1,-2,-3)], na.rm = TRUE)
f.avg <- rowMeans(x.annot_meth[which(x.annot_meth$sex == "F"),c(-1,-2,-3)], na.rm = TRUE)

hist(m.avg)
hist(f.avg)

m.avgs.Xsite <- colMeans(x.annot_meth[which(x.annot_meth$sex == "M"),c(-1,-2,-3)], na.rm = TRUE)
f.avgs.Xsite <- colMeans(x.annot_meth[which(x.annot_meth$sex == "F"),c(-1,-2,-3)], na.rm = TRUE)

xSite.total <- data.frame("means" = c(m.avgs.Xsite, f.avgs.Xsite), "sex" = c(rep("M", length(m.avgs.Xsite)), rep("F", length(f.avgs.Xsite))))

f <- ggplot(data = xSite.total, aes(x = sex, y = means, fill = sex)) + 
    geom_violin() + geom_jitter(color = "black", size = 0.5, alpha = 0.3) +
    labs(title = "FHRX", y = "Average Meth %", fill = "Sex") + 
    theme(title = element_text(size = 16), axis.title = element_text(size = 16), 
          axis.title.x = element_blank(), legend.key.size = unit(0.3, "inches"), 
          legend.text = element_text(size = 16))
f

ggsave(f, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/FHRX_sexDMC.png", 
       height = 4, width = 5)

fXsite.sites_diff <- na.omit(data.frame("maleAvg" = m.avgs.Xsite, "femAvg" = f.avgs.Xsite, "diff" = ( f.avgs.Xsite - m.avgs.Xsite )))

fXsite.sites_diff.hist <- ggplot(data = fXsite.sites_diff, aes(x = fXsite.sites_diff$diff)) + 
    geom_histogram(binwidth = 5, fill = "grey40", color = "black") + labs(x = "% Meth Difference")
fXsite.sites_diff.box <- ggplot(data = fXsite.sites_diff, aes(y = fXsite.sites_diff$diff)) + 
    geom_boxplot(fill = "red", color = "black") + labs(y = "% Meth Difference") +
    theme(axis.ticks.x = element_blank())

ggarrange(fXsite.sites_diff.hist, fXsite.sites_diff.box + rremove("x.text"), ncol = 2, nrow = 1, align = "v",
          widths = c(1.5,.5))

### Plot out this region

x.site.larger <- GRanges(seqnames = "scaffold_7", ranges = IRanges(start = 62300000, end = 62800000), strand = "*")
x.larger.overlapped <- subsetByOverlaps(meth.gr, x.site.larger, type = "any", maxgap = 0L, ignore.strand = TRUE)
x.larger.overlap.frame <- as.data.frame(x.larger.overlapped)


fhrx.df <- data.frame("site" = start(x.larger.overlapped),
                      "MaleAvg" = rowMeans(x.larger.overlap.frame[,(which(annotations$sex == "M")+5)], na.rm = TRUE),
                      "FemaleAvg" = rowMeans(x.larger.overlap.frame[,(which(annotations$sex == "F")+5)], na.rm = TRUE))

fhrx.long <- tidyr::pivot_longer(fhrx.df, cols = c(2, 3), names_to = "Sex", values_to = "AvgMeth")

library(dplyr)

fhrx.long.binned <- fhrx.long %>% mutate(Bins = ntile(site, 250))
    
New <- fhrx.long.binned %>% 
    group_by(Sex, Bins, .add = TRUE) %>%
    summarize(WindowAvgs = mean(AvgMeth, na.rm = TRUE))


g <- ggplot(New, aes(x = Bins, y = WindowAvgs, fill = WindowAvgs)) + 
    facet_wrap(~Sex, nrow = 2) +
    geom_col(width = 1) + 
    #geom_smooth(method = "loess", span = 0.2, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "forestgreen", high = "red", midpoint = 50) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.position = "none") + 
    xlab("Coordinates in Mbp") +
    scale_x_continuous(breaks = c(seq(0, 250, 50)), labels = c(62.3, 62.4, 62.5, 62.6, 62.7, 62.8))
    
g

ggsave(g, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/FHRX_rawPlot_columns.png",
       width = 6, height = 2)

## Now plot out the genes in this region

genes <- import.gff("C:/Users/sheal/Desktop/Parrott_Lab/Anole_Data/AnoSag2.1/CDS.gtf")
CGI <- import.bed("C:/Users/sheal/Desktop/Parrott_Lab/Anole_Data/AnoSag2.1/CGIslands.bed")

roi.genes <- subsetByOverlaps(genes, x.site.larger, type = "within")
roi.cgi <- subsetByOverlaps(CGI, x.site.larger, type = "within")

export.gff3(roi.genes, con = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/data/FHRX.gff")
export.bed(roi.cgi, con = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/data/FHRX_CGI.bed")

library(gggenomes)

# a minimal seq track
s0 <- tibble(
  seq_id = c("FHRX"),
  length = c(500000)
)


g0 <- read_feats("C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/data/FHRX.gff")

g1 <- tibble(
  seq_id = rep("FHRX", 5),
  start = g0$start - start(x.site.larger),
  end =  g0$end - start(x.site.larger),
  label = g0$gene_name
)

f0 <- read_feats("C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/data/FHRX_CGI.bed")

f1 <- tibble(
  seq_id = rep("FHRX", length(f0$file_id)),
  start = f0$start - start(x.site.larger),
  end =  f0$end - start(x.site.larger)
)

f2 <- read_feats("C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/data/Sex_DMRs_smoothed200.bed")

f3 <- tibble(
  seq_id = rep("FHRX", length(f2$file_id)),
  start = f2$start - start(x.site.larger),
  end =  f2$end - start(x.site.larger)
)

g.sub <- gggenomes(seqs = s0, genes = g1[-1,], feats = f3) +
    geom_seq() +
    geom_gene(size = 2, fill = "blue") +
    geom_feat(color = "magenta", linewidth = 3) +
    geom_gene_label(aes(label =  g0$gene_name[-1]), size = 2, vjust = 0)

g.sub

ggsave(g.sub, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig2-DSSsex/images/FHRX_genes.png",
       width = 5, height = 1)



## Quantify relative methylation of males and females at FHRX to chromosome background

meth.x.gr <- meth.gr[which(seqnames(meth.gr) == "scaffold_7")]
avg.male <- mean(as.matrix(mcols(meth.x.gr))[,which(annotations$sex == "M")], na.rm = TRUE)
avg.female <- mean(as.matrix(mcols(meth.x.gr))[,which(annotations$sex == "F")], na.rm = TRUE)

fhmx.overlapped <- subsetByOverlaps(meth.gr, x.site, type = "within", ignore.strand = TRUE)
avg.male.fhmx <- mean(as.matrix(mcols(fhmx.overlapped))[,which(annotations$sex == "M")], na.rm = TRUE)
avg.female.fhmx <- mean(as.matrix(mcols(fhmx.overlapped))[,which(annotations$sex == "F")], na.rm = TRUE)

