library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(ggpubr)
library(dplyr)
library(psych)
library(ComplexHeatmap)

theme_set(theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")))

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F")))

age.data <- read.table("C:/Users/sheal/Desktop/AnoleAge/Fig1.1-DSSLogAge/data/LogAgeSites_smoothed200.bed",
                       sep = "\t", header = FALSE)


age.sites <- GRanges(seqnames = age.data$V1,
                     ranges = IRanges(start = age.data$V2, end = age.data$V3), strand = "*")

filt <- as.data.frame(fread("E:/anole/CpG/processed/labeled_total80perc.tab", header = TRUE, sep = "\t"))

meth.t <- t(filt[,c(5:41)])
all.cor <- corr.test(meth.t, log(annotations$ages), use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)
cor.df <- data.frame("r" = all.cor$r)


mean(cor.df$r, na.rm = TRUE); sd(cor.df$r, na.rm = TRUE) / sqrt(length(na.omit(cor.df$r))) ## Mean: 0.08528, +/- 0.00017 SE

f.cor <- corr.test(meth.t[which(annotations$sex == "F"),],
                   log(annotations$ages)[which(annotations$sex == "F")], 
                   use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)

m.cor <- corr.test(meth.t[which(annotations$sex == "M"),],
                   log(annotations$ages)[which(annotations$sex == "M")], 
                   use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)

cor.df.sex <- data.frame("Female" = f.cor$r, "Male" = m.cor$r)
cor.sex.long <- tidyr::pivot_longer(cor.df.sex, cols = c(1,2), names_to = "Sex", values_to = "r")

b <- ggplot(cor.sex.long, aes(x = r, fill = Sex)) + 
        ylab("Number of CpGs") + xlab("Age Correlation Coefficient") + 
        geom_histogram(position = "identity", alpha = 0.4, bins = 99) + 
        geom_vline(aes(xintercept = mean(f.cor$r, na.rm = TRUE)), color = "firebrick1", linetype = "dashed", linewidth = 1.1) +
        geom_vline(aes(xintercept = mean(m.cor$r, na.rm = TRUE)), color = "dodgerblue", linetype = "dashed", linewidth = 1.1)
b

ggsave(b, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/global_sexdiff_hist.png",
       width = 4.5, height = 4)

mean(f.cor$r, na.rm = TRUE); mean(m.cor$r, na.rm = TRUE)

t.test(f.cor$r, m.cor$r, paired = TRUE)

table <- data.frame("Male" = c(length(which(m.cor$r > 0.5)), length(which(m.cor$r < -0.5))),
                    "Female" = c(length(which(f.cor$r > 0.5)), length(which(f.cor$r < -0.5))))
rownames(table) <- c("Positive", "Negative")
fisher.test(table)



#### Now repeat for aDMCs

meth_file <- as.data.frame(fread(file = "E:/anole/CpG/processed/labeled_allSites_stranded.tab", header = TRUE, sep = "\t"))

meth.gr <- GRanges(seqnames = meth_file$chr, 
                   ranges = IRanges(start = meth_file$start, end = meth_file$end),
                   strand = "*", mcols = meth_file[,c(5:41)] )

age.data <- age.data[-4671,] ## For some reason this site doesn't exist in log age set??

ageSites.gr <- GRanges(age.data$V1, IRanges(age.data$V2, age.data$V3), strand = "*")


meth.ageSites <- subsetByOverlaps(meth.gr, ageSites.gr, type = "end")

ageSites.mat <- as.data.frame(t(as.data.frame(meth.ageSites)[,-c(1,2,3,4,5)]))
colnames(ageSites.mat) <- paste0(age.data$V1, ".", age.data$V2)

ageSites.df.clean <- data.frame(annotations, ageSites.mat)

all.cor <- corr.test(ageSites.mat, log(annotations$ages), use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)
cor.admc.df <- data.frame("r" = all.cor$r)

ggplot(cor.admc.df, aes(x = r)) + geom_histogram(bins = 99, fill = "grey", color = "grey20") +
        ylab("Number of CpGs") + xlab("Age Correlation Coefficient")

f.cor <- corr.test(ageSites.mat[which(annotations$sex == "F"),],
                   log(annotations$ages)[which(annotations$sex == "F")], 
                   use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)

m.cor <- corr.test(ageSites.mat[which(annotations$sex == "M"),],
                   log(annotations$ages)[which(annotations$sex == "M")], 
                   use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)

cor.admc.sex <- data.frame("Female" = f.cor$r, "Male" = m.cor$r)

ggplot(cor.admc.sex, aes(x = Female, y = Male)) + 
    geom_point() + geom_smooth(method = "lm") + 
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1.5)


cor.admc.sex.long <- tidyr::pivot_longer(cor.admc.sex, cols = c(1,2), names_to = "Sex", values_to = "r")

d <- ggplot(cor.admc.sex.long, aes(x = r, fill = Sex)) + 
        ylab("Number of CpGs") + xlab("Age Correlation Coefficient") + 
        geom_histogram(position = "identity", alpha = 0.4, bins = 99) + 
        geom_vline(aes(xintercept = mean(f.cor$r, na.rm = TRUE)), color = "firebrick1", linetype = "dashed", linewidth = 1.1) +
        geom_vline(aes(xintercept = mean(m.cor$r, na.rm = TRUE)), color = "dodgerblue", linetype = "dashed", linewidth = 1.1)
d

ggsave(d, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/aDMC_sexdiff_hist.png",
       width = 4.5, height = 4)

