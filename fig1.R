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

#### DMC overlap and quantify ---------------------------------------------------------------------------------------------

age.data <- read.table("~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/data/LogAgeSites_smoothed200.bed",
                       sep = "\t", header = FALSE) # Read in age-related DMCs derived from DSS model

# Manually generate dataframe of sample metadata 

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F")))

# Convert age sites to GRanges object

age.sites <- GRanges(seqnames = age.data$V1,
                     ranges = IRanges(start = age.data$V2, end = age.data$V3), strand = "*")

# Import methylation matrix generated and filtered in MethylKit

filt <- as.data.frame(fread("E:/anole/CpG/processed/labeled_total80perc.tab", header = TRUE, sep = "\t"))

# Quantify and plot missingness by sample

missingNo <- colSums(is.na(filt[,c(5:41)]))

hist(missingNo/length(filt$chr))

# Format methylation matrix and calculate pearson correlation coefficients between ln(age) and each site's methylation

meth.t <- t(filt[,c(5:41)])
all.cor <- corr.test(meth.t, log(annotations$ages), use = "pairwise", method = "pearson", adjust = "fdr", ci = FALSE)
cor.df <- data.frame("r" = all.cor$r)

# generate figure showing global cor coeffs

b1 <- ggplot(cor.df, aes(x = r)) + geom_histogram(bins = 99, fill = "grey", color = "grey20") +
        ylab("Number of CpGs") + xlab("Age Correlation Coefficient")
b1

ggsave(b1, filename = "~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/images/total80perc_logAgeCor.svg", device = "svg",
       width = 5, height = 5)

mean(cor.df$r, na.rm = TRUE); sd(cor.df$r, na.rm = TRUE) / sqrt(length(na.omit(cor.df$r))) ## Mean: 0.08528, +/- 0.00017 SE

## Import unfiltered version of methylation matrix for overlapping with DSS sites

meth_file <- as.data.frame(fread(file = "E:/anole/CpG/processed/labeled_allSites_stranded.tab", header = TRUE, sep = "\t"))

# calculate global methylation values for each sample, and generate / test linear model of figure 1A

global <- data.frame(annotations,
                     "Avgs" = colMeans(filt[,-c(1:4)], na.rm = TRUE))

a <- ggplot(data = global, aes(x = ages, y = Avgs, color = sex)) + 
  geom_point(size = 2) + ylab("Global Average Methylation %") + geom_smooth(aes(group = 1), method = "lm", color = "black") +
  xlab("Age in Months") + theme(legend.title = element_text(size = 16, face = "bold",),
                                legend.text = element_text(size = 14))
a

ggsave(a, filename = "~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/images/GlobalCpGMeth.svg", device = "svg",
       width = 5, height = 5)

lm.glob <- lm(formula = Avgs ~ ages, data = global[-c(4,37),])
summary(lm.glob)

# t-test for sex differences in global methylation

t.glob <- t.test(formula = Avgs ~ sex, data = global[-c(4,37),])
t.glob

# convert unfiltered meth matrix to GRanges

meth.gr <- GRanges(seqnames = meth_file$chr, 
                   ranges = IRanges(start = meth_file$start, end = meth_file$end),
                   strand = "*", mcols = meth_file[,c(5:41)] )

# overlap GRanges objects and format dataframe

overlapped <- subsetByOverlaps(meth.gr, age.sites, type = "end", ignore.strand = TRUE)

rm(meth.gr, meth_file)
gc()

overlap.frame <- as.data.frame(overlapped)

sites <- paste0(overlap.frame$seqnames, ".", overlap.frame$start, ".", overlap.frame$end)

formatted_meth <- as.data.frame(t(overlap.frame[,c(6:42)]))

annot_meth <- cbind(annotations, formatted_meth)

## plot example site

ggplot(data = annot_meth, aes(x = log(ages), y = V408, group = ages)) + geom_point() +
    ylab("Methylation Frequency") + xlab("Log(Age in Months)") + labs(title = "Chromosome 1: 191,672,618-191,672,619") + 
    geom_smooth(method = "lm", aes(group = 1), alpha = 0.2)

ggplot(data = annot_meth, aes(x = log(ages), y = V408, group = ages)) + geom_boxplot() + geom_point() +
    ylab("Methylation Frequency") + xlab("Log(Age in Months)") + labs(title = "Chromosome 1: 191,672,618-191,672,619") + 
    geom_smooth(method = "lm", aes(group = 1), alpha = 0.2)

ggplot(data = annot_meth, aes(x = log(ages), y = V408, group = ages)) + geom_boxplot() + geom_point() +
    ylab("Methylation Frequency") + xlab("Log(Age in Months)") + labs(title = "Chromosome 1: 191,672,618-191,672,619") + 
    geom_smooth(method = "lm", aes(group = 1), alpha = 0.2) + facet_wrap("sex")

## Get corr coeffs for age-related sites

cpg.ageCor <- corr.test(formatted_meth, log(annotations$ages), use = "pairwise", method = "spearman", adjust = "fdr", ci = FALSE)
cpg.ageCor <- na.omit(data.frame("chr" = overlap.frame$seqnames, "loc" = overlap.frame$start, "r" = cpg.ageCor$r, "p.adj" = cpg.ageCor$p.adj))

# Generate plot of aDMC cor coeffs

b2 <- ggplot(data = cpg.ageCor, aes(x = r)) + geom_histogram(bins = 99, color = "black", fill = "grey") + 
    ylab("Number of CpGs") + xlab("Age Correlation Coefficients")
b2

ggsave(b2, filename = "~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/images/Signif_CpG_CorPlot.svg", device = "svg",
       width = 5, height = 5)

mean(cpg.ageCor$r, na.rm = TRUE); sd(cpg.ageCor$r, na.rm = TRUE) / sqrt(length(na.omit(cpg.ageCor$r))) ## Mean: -0.1961, +/- 0.0080 SE

### Generate combined global and aDMC plot (Fig 2.B)

# Righthand axis is scaled by a factor of 343 to match lefthand scale

b <- ggplot(data = cor.df, aes(x = r)) + 
      geom_histogram(fill = "grey40", color = "black", bins = 99) +
      geom_histogram(data = cpg.ageCor, aes(x = r, y = 343*after_stat(count)), 
                     fill = "firebrick", color = "black", alpha = 0.5, bins = 99) +
      scale_y_continuous(name = "Number of Total CpGs", 
                         sec.axis = sec_axis(transform=~./343, name = "Number of Age-Related CpGs")) +
      theme(axis.title.y.left = element_text(color = "grey40"),
            axis.title.y.right = element_text(color = "firebrick"),
            axis.text.y.left = element_text(color = "grey40"), 
            axis.text.y.right = element_text(color = "firebrick")) + 
      xlab("Age Correlation Coefficients")

b

ggsave(b, filename = "~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/images/CpG_CorPlot.svg", device = "svg",
       width = 6, height = 5)

#### calculate and Plot DMCs per Mbp on each chromosome 

table <- as.data.table(table(age.data$V1))
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
                36225908,
                23386130,
                20517538)

## Calculate and plot site density per Mbp
table.sorted <- table.sorted[which(table.sorted$V1 %in% paste0("scaffold_", c(1:14))),]
table.sorted$size <- chromSizes
table.sorted$sizeMb <- table.sorted$size / 1000000
table.sorted$siteDensity <- table.sorted$N / table.sorted$sizeMb
table.sorted$chrName <- paste0("chr", c(1:6, "X", 8:14))
table.sorted$chrName <- factor(table.sorted$chrName, levels = table.sorted$chrName)

# generate and save Fig 1D

d <- ggplot(table.sorted, aes(x = siteDensity, y = chrName)) + geom_col(fill = "navy") + 
    xlab("Age-Related Sites per Mbp") + ylab("")
d

ggsave(d, filename = "~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/images/AgeSites_byMBp.svg", device = "svg",
       width = 4, height = 5)


## Calculate total CpG density per chromosome

CGsites <- fread("~/Parrott_Lab/Anole_Data/AnoSag2.1/CGsites.bed", header = F)
table.cg <- as.data.frame(table(CGsites$V1))
rm(CGsites)

## Table scaffold counts and order by number
table.cg.sorted <- table.cg[stringr::str_order(table.cg$Var1, numeric = TRUE, decreasing = FALSE),] 
table.cg.sorted <- table.cg.sorted[c(1:14),]
table.cg.sorted$Var1 <- factor(table.cg.sorted$Var1, levels = stringr::str_sort(table.cg.sorted$Var1, numeric = TRUE, decreasing = FALSE))
table.cg.sorted$size <- chromSizes
table.cg.sorted$sizeMb <- table.cg.sorted$size / 1000000
table.cg.sorted$siteDensity <- table.cg.sorted$Freq / table.cg.sorted$sizeMb

# t-test of macro vs micro-chromosomes
t.test(table.cg.sorted$siteDensity[1:6], table.cg.sorted$siteDensity[8:14])

ggplot(table.cg.sorted, aes(x = siteDensity, y = Var1)) + geom_col(fill = "navy") + 
    ylab("Scaffold") + xlab("Number of CpGs per Mbp") + theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))


## Calculate number of measured CpGs per chromosome

meth_file <- as.data.frame(fread(file = "E:/anole/CpG/processed/labeled_allSites_stranded.tab", header = TRUE, sep = "\t"))
table.mes <- as.data.frame(table(meth_file$chr))
rm(meth_file)
table.mes.sorted <- table.mes[stringr::str_order(table.mes$Var1, numeric = TRUE, decreasing = FALSE),] ## Table scaffold counts and order by number
table.mes.sorted <- table.mes.sorted[c(1:14),]
table.mes.sorted$Var1 <- factor(table.mes.sorted$Var1, levels = stringr::str_sort(table.mes.sorted$Var1, numeric = TRUE, decreasing = FALSE))
table.mes.sorted$size <- chromSizes
table.mes.sorted$sizeMb <- table.mes.sorted$size / 1000000
table.mes.sorted$siteDensity <- table.mes.sorted$Freq / table.mes.sorted$sizeMb

ggplot(table.mes.sorted, aes(x = siteDensity, y = Var1)) + geom_col(fill = "navy") + 
    ylab("Scaffold") + xlab("Measured CpGs per Mbp") + theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))


tot <- merge(table.sorted, table.cg.sorted[,c(-3,-4)], by.x = "V1", by.y = "Var1", suffixes = c(".AgeSites", ".AllCpGs"))
tot <- data.frame(tot, "N.Measured" = table.mes.sorted$Freq, "siteDensity.Measured" = table.mes.sorted$siteDensity)

s1 <- ggplot(data = tot, aes(x = N.Measured, y = N, label = V1)) + geom_smooth(method = "lm", alpha = 0.2) + geom_point(size = 4) + 
                   ggrepel::geom_text_repel() + 
                   xlab("Number of CpGs measured") + ylab("Number of Age-Related sites") + theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
s1

ggsave(s1, filename = "~/Parrott_Lab/AnoleAge/Supp/AgeSites_byMBp.svg", device = "svg",
       width = 6, height = 6)

## binomial test for enrichment on Chr3 (actual aDMCs versus prediction by regression on number of measured CpGs)
lm.binom <- lm(data = tot, formula = N ~ N.Measured)
predicted.chr3 <- predict(lm.binom)[3]
actual.chr3 <- tot$N[3]

binom.test(x = actual.chr3, n = tot$N.Measured[3], p = ( predicted.chr3 / tot$N.Measured[3] ),
           alternative = "two.sided",
           conf.level = 0.95) ## p < 2.2e-16 for enrichment on chr3


# Age DMC enrichment ------------------------------------------------------


## First positively correlated sites 

# Import CpG island gff file generated using EMBOSS's cpgplot
islands <- import("~/Parrott_Lab/Anole_Data/AnoSag2.1/CGI.gff", format = "gff")

# define shores, shelves, seas

shores.up <- trim(flank(islands, 2000))
shores.down <- trim(flank(islands, 2000, start = FALSE))
shores <- c(shores.up, shores.down)
shores.clean <- reduce(GenomicRanges::setdiff(shores, islands))

shelves.up <- trim(flank(shores.up, 2000))
shelves.down <- trim(flank(shores.down, 2000, start = FALSE))
shelves <- c(shelves.up, shelves.down)
shelves.no.islands <- GenomicRanges::setdiff(shelves, islands)
shelves.clean <- reduce(GenomicRanges::setdiff(shelves.no.islands, shores))


not.seas <- c(islands,shores,shelves)
seas <- gaps(not.seas)

## CGdensity enrichment 

## Find background likelihood that CpGs are found in each context

CGsites <- fread("E:/anole/CpG/processed/labeled_allSites_stranded.tab", sep = "\t", header = TRUE)

CGsites.clean <- GRanges(seqnames = CGsites$chr, 
                         ranges = IRanges(start = CGsites$start, end = CGsites$end),
                         strand = "*")

totalCG <- length(CGsites.clean)


sites.in.islands <- countOverlaps(CGsites.clean, islands, type = "within")
sum(sites.in.islands) # 2,121,767
cat("Proportion of CG sites that fall in islands: ", sum(sites.in.islands) / length(CGsites.clean), "\n")


sites.in.shores <- countOverlaps(CGsites.clean, shores.clean, type = "within")
sum(sites.in.shores) # 4,858,599
cat("Proportion of CG sites that fall in shores: ", sum(sites.in.shores) / length(CGsites.clean), "\n")


sites.in.shelves <- countOverlaps(CGsites.clean, shelves.clean, type = "within")
sum(sites.in.shelves) # 3,205,908
cat("Proportion of CG sites that fall in shelves: ", sum(sites.in.shelves) / length(CGsites.clean), "\n")


sites.in.seas <- countOverlaps(CGsites.clean, seas, type = "within")
sum(sites.in.seas) # 23,179,866
cat("Proportion of CG sites that fall in seas: ", sum(sites.in.seas) / length(CGsites.clean), "\n")

# Import aDMCs w/ correlation coeffs

CPG <- import("~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/data/AgeSites_cor.bed", format = "bed")
CPG <- CPG[which(CPG$name > 0),]

totalCpG <- length(CPG)

cpg.in.islands <- countOverlaps(CPG, islands, type = "within")
sum(cpg.in.islands) # 58
cat("Proportion of (+) Age-Associated CpGs that fall in islands: ", sum(cpg.in.islands) / length(CPG), "\n")


cpg.in.shores <- countOverlaps(CPG, shores.clean, type = "within")
sum(cpg.in.shores) # 327
cat("Proportion of (+) Age-Associated CpGs that fall in shores: ", sum(cpg.in.shores) / length(CPG), "\n")


cpg.in.shelves <- countOverlaps(CPG, shelves.clean, type = "within")
sum(cpg.in.shelves) # 212
cat("Proportion of (+) Age-Associated CpGs that fall in shelves: ", sum(cpg.in.shelves) / length(CPG), "\n")


cpg.in.seas <- countOverlaps(CPG, seas, type = "within")
sum(cpg.in.seas) # 1160
cat("Proportion of (+) Age-Associated CpGs that fall in seas: ", sum(cpg.in.seas) / length(CPG), "\n")


pos.density.results <- data.frame("Context" = c("islands", "shores", "shelves", "seas","islands", "shores", "shelves", "seas"), 
                      "Values" = c(sum(sites.in.islands), sum(sites.in.shores),
                                   sum(sites.in.shelves), sum(sites.in.seas),
                                   sum(cpg.in.islands), sum(cpg.in.shores), 
                                   sum(cpg.in.shelves), sum(cpg.in.seas)),
                      "Proportions" = c(sum(sites.in.islands)/totalCG, sum(sites.in.shores)/totalCG,
                                   sum(sites.in.shelves)/totalCG, sum(sites.in.seas)/totalCG,
                                   sum(cpg.in.islands)/totalCpG, sum(cpg.in.shores)/totalCpG, 
                                   sum(cpg.in.shelves)/totalCpG, sum(cpg.in.seas)/totalCpG),
                      "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

pos.density.results$Context <- factor(pos.density.results$Context, levels=c("islands","shores","shelves", "seas"))

enrichment.cgdensity <- ggplot(data = pos.density.results, mapping = aes(Context, Proportions)) + 
                               geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
                               scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
                               labs(title = "Enrichment of CpGs with (+) Age-Association") +
                               ylab("Percentage of Total CpGs") +
                               xlab("") +
                               theme(axis.line = element_line(colour = "black")) +
                               font("xlab", size = 18, color = "black") +
                               font("ylab", size = 18, color = "black") +
                               font("xy.text", size = 18, color = "black") +
                               font("title", size = 15, color = "black", face = "bold") +
                               geom_signif(stat = "identity", data = data.frame(x=c(0.65), 
                                                                                xend=c(1.3),
                                                                                y=c(0.09), 
                                                                                annotation=c("**")),
                                           aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))

enrichment.cgdensity


## Statistical test:

#Islands
binom.test(x = sum(cpg.in.islands), n = totalCpG, p = ( sum(sites.in.islands)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 5.74e-05 **

#Shores
binom.test(x = sum(cpg.in.shores), n = totalCpG, p = ( sum(sites.in.shores)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 2.951e-06 **

#Shelves
binom.test(x = sum(cpg.in.shelves), n = totalCpG, p = ( sum(sites.in.shelves)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.0006636 *

#Seas
binom.test(x = sum(cpg.in.seas), n = totalCpG, p = ( sum(sites.in.seas)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 0.002084*


## Generate genic contexts


exons <- import.gff2("~/Parrott_Lab/Anole_Data/AnoSag2.1/AnoSag2.1_exons.gtf")
exons.destranded <- reduce(exons, ignore.strand = TRUE)

CDS <- import.gff2("~/Parrott_Lab/Anole_Data/AnoSag2.1/CDS.gtf")
CDS.destranded <- reduce(CDS, ignore.strand = TRUE)
full.CDS <- GenomicRanges::union(CDS.destranded, exons.destranded)

intergenic <- gaps(CDS.destranded)
intergenic.clean <- GenomicRanges::setdiff(intergenic, exons.destranded)

exons.comp <- gaps(exons)

introns <- GenomicRanges::setdiff(CDS.destranded, exons.destranded)

promoter <- reduce(promoters(CDS), ignore.strand = TRUE)

#### Genic context Enrichment  

## Find background likelihood that CpGs are found in each context

sites.in.promo <- countOverlaps(CGsites.clean, promoter, type = "within")
sum(sites.in.promo) # 1,243,018
cat("Proportion of CG sites that fall in promoters: ", sum(sites.in.promo) / length(CGsites.clean), "\n")


sites.in.exons <- countOverlaps(CGsites.clean, exons.destranded, type = "within")
sum(sites.in.exons) # 2,213,587
cat("Proportion of CG sites that fall in exons: ", sum(sites.in.exons) / length(CGsites.clean), "\n")


sites.in.introns <- countOverlaps(CGsites.clean, introns, type = "within")
sum(sites.in.introns) # 12,543,534
cat("Proportion of CG sites that fall in introns: ", sum(sites.in.introns) / length(CGsites.clean), "\n")


sites.in.intergenic <- countOverlaps(CGsites.clean, intergenic.clean, type = "within")
sum(sites.in.intergenic) # 18,547,797
cat("Proportion of CG sites that fall in intergenic space: ", sum(sites.in.intergenic) / length(CGsites.clean), "\n")

# get likelihood of age DMCs in each context

cpg.in.promo <- countOverlaps(CPG, promoter, type = "within")
sum(cpg.in.promo) # 59
cat("Proportion of (+) Age-Associated CpGs that fall in promoters: ", sum(cpg.in.promo) / length(CPG), "\n")

cpg.in.exons <- countOverlaps(CPG, exons.destranded, type = "within")
sum(cpg.in.exons) # 45
cat("Proportion of (+) Age-Associated CpGs that fall in exons: ", sum(cpg.in.exons) / length(CPG), "\n")


cpg.in.introns <- countOverlaps(CPG, introns, type = "within")
sum(cpg.in.introns) # 669
cat("Proportion of (+) Age-Associated CpGs that fall in introns: ", sum(cpg.in.introns) / length(CPG), "\n")


cpg.in.intergenic <- countOverlaps(CPG, intergenic.clean, type = "within")
sum(cpg.in.intergenic) # 1043
cat("Proportion of (+) Age-Associated CpGs that fall in intergenic space: ", sum(cpg.in.intergenic) / length(CPG), "\n")

# combine results

pos.genic.results <- data.frame("Context" = c("promoters", "exons", "introns", "intergenic", "promoters", "exons", "introns", "intergenic"), 
                      "Values" = c(sum(sites.in.promo), sum(sites.in.exons), 
                                   sum(sites.in.introns), sum(sites.in.intergenic),
                                   sum(cpg.in.promo), sum(cpg.in.exons), 
                                   sum(cpg.in.introns), sum(cpg.in.intergenic)),
                      "Proportions" = c(sum(sites.in.promo)/totalCG, sum(sites.in.exons)/totalCG, 
                                        sum(sites.in.introns)/totalCG, sum(sites.in.intergenic)/totalCG,
                                        sum(cpg.in.promo)/totalCpG, sum(cpg.in.exons)/totalCpG, 
                                        sum(cpg.in.introns)/totalCpG, sum(cpg.in.intergenic)/totalCpG),
                      "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

pos.genic.results$Context <- factor(pos.genic.results$Context, levels=c("promoters", "exons", "introns", "intergenic"))


enrichment.cgdensity <- ggplot(data = pos.genic.results, mapping = aes(Context, Proportions)) + 
                               geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
                               scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
                               labs(title = "Enrichment of CpGs with (+) Age-Association") +
                               ylab("Percentage of Total CpGs") +
                               xlab("") +
                               theme(axis.line = element_line(colour = "black")) +
                               font("xlab", size = 18, color = "black") +
                               font("ylab", size = 18, color = "black") +
                               font("xy.text", size = 18, color = "black") +
                               font("title", size = 15, color = "black", face = "bold")

enrichment.cgdensity


## Statistical test:

#Promoters
binom.test(x = sum(cpg.in.promo), n = totalCpG, p = ( sum(sites.in.promo)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.3815

#Exons
binom.test(x = sum(cpg.in.exons), n = totalCpG, p = ( sum(sites.in.exons)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 1.445e-14

#Introns
binom.test(x = sum(cpg.in.introns), n = totalCpG, p = ( sum(sites.in.introns)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 0.6755

#Intergenic
binom.test(x = sum(cpg.in.intergenic), n = totalCpG, p = ( sum(sites.in.intergenic)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.001413


## Now negatively correlated sites


CPG <- import("~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/data/AgeSites_cor.bed", format = "bed")
CPG <- CPG[which(CPG$name < 0),]
totalCpG <- length(CPG)

cpg.in.islands <- countOverlaps(CPG, islands, type = "within")
sum(cpg.in.islands) # 34
cat("Proportion of (-) Age-Associated CpGs that fall in islands: ", sum(cpg.in.islands) / length(CPG), "\n")


cpg.in.shores <- countOverlaps(CPG, shores.clean, type = "within")
sum(cpg.in.shores) # 320
cat("Proportion of (-) Age-Associated CpGs that fall in shores: ", sum(cpg.in.shores) / length(CPG), "\n")


cpg.in.shelves <- countOverlaps(CPG, shelves.clean, type = "within")
sum(cpg.in.shelves) # 324
cat("Proportion of (-) Age-Associated CpGs that fall in shelves: ", sum(cpg.in.shelves) / length(CPG), "\n")


cpg.in.seas <- countOverlaps(CPG, seas, type = "within")
sum(cpg.in.seas) # 2546
cat("Proportion of (-) Age-Associated CpGs that fall in seas: ", sum(cpg.in.seas) / length(CPG), "\n")


neg.density.results <- data.frame("Context" = c("islands", "shores", "shelves", "seas","islands", "shores", "shelves", "seas"), 
                      "Values" = c(sum(sites.in.islands), sum(sites.in.shores),
                                   sum(sites.in.shelves), sum(sites.in.seas),
                                   sum(cpg.in.islands), sum(cpg.in.shores), 
                                   sum(cpg.in.shelves), sum(cpg.in.seas)),
                      "Proportions" = c(sum(sites.in.islands)/totalCG, sum(sites.in.shores)/totalCG,
                                       sum(sites.in.shelves)/totalCG, sum(sites.in.seas)/totalCG,
                                       sum(cpg.in.islands)/totalCpG, sum(cpg.in.shores)/totalCpG, 
                                       sum(cpg.in.shelves)/totalCpG, sum(cpg.in.seas)/totalCpG),
                      "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

neg.density.results$Context <- factor(neg.density.results$Context, levels=c("islands","shores","shelves", "seas"))


enrichment.cgdensity <- ggplot(data = neg.density.results, mapping = aes(Context, Proportions)) + 
                               geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
                               scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
                               labs(title = "Enrichment of CpGs with (-) Age-Association") +
                               ylab("Percentage of Total CpGs") +
                               xlab("") +
                               theme(axis.line = element_line(colour = "black")) +
                               font("xlab", size = 18, color = "black") +
                               font("ylab", size = 18, color = "black") +
                               font("xy.text", size = 18, color = "black") +
                               font("title", size = 15, color = "black", face = "bold") +
                               geom_signif(stat = "identity", data = data.frame(x=c(0.65, 1.65, 3.65), 
                                                                                xend=c(1.3, 2.3, 4.3),
                                                                                y=c(0.1, 0.18, 0.73), 
                                                                                annotation=c(" *** ", "*", "   **   ")),
                                           aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))

enrichment.cgdensity


## Statistical test:

#Islands
binom.test(x = sum(cpg.in.islands), n = totalCpG, p = ( sum(sites.in.islands)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 3.839e-12

#Shores
binom.test(x = sum(cpg.in.shores), n = totalCpG, p = ( sum(sites.in.shores)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 8.054e-15

#Shelves
binom.test(x = sum(cpg.in.shelves), n = totalCpG, p = ( sum(sites.in.shelves)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.5034

#Seas
binom.test(x = sum(cpg.in.seas), n = totalCpG, p = ( sum(sites.in.seas)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 2.2e-16

cpg.in.promo <- countOverlaps(CPG, promoter, type = "within")
sum(cpg.in.promo) # 51
cat("Proportion of (+) Age-Associated CpGs that fall in promoters: ", sum(cpg.in.promo) / length(CPG), "\n")

cpg.in.exons <- countOverlaps(CPG, exons.destranded, type = "within")
sum(cpg.in.exons) # 173
cat("Proportion of (-) Age-Associated CpGs that fall in exons: ", sum(cpg.in.exons) / length(CPG), "\n")


cpg.in.introns <- countOverlaps(CPG, introns, type = "within")
sum(cpg.in.introns) # 1347
cat("Proportion of (-) Age-Associated CpGs that fall in introns: ", sum(cpg.in.introns) / length(CPG), "\n")


cpg.in.intergenic <- countOverlaps(CPG, intergenic.clean, type = "within")
sum(cpg.in.intergenic) # 1702
cat("Proportion of (-) Age-Associated CpGs that fall in intergenic space: ", sum(cpg.in.intergenic) / length(CPG), "\n")


neg.genic.results <- data.frame("Context" = c("promoters", "exons", "introns", "intergenic", "promoters", "exons", "introns", "intergenic"), 
                      "Values" = c(sum(sites.in.promo), sum(sites.in.exons), 
                                   sum(sites.in.introns), sum(sites.in.intergenic),
                                   sum(cpg.in.promo), sum(cpg.in.exons), 
                                   sum(cpg.in.introns), sum(cpg.in.intergenic)),
                      "Proportions" = c(sum(sites.in.promo)/totalCG, sum(sites.in.exons)/totalCG, 
                                        sum(sites.in.introns)/totalCG, sum(sites.in.intergenic)/totalCG,
                                        sum(cpg.in.promo)/totalCpG, sum(cpg.in.exons)/totalCpG, 
                                        sum(cpg.in.introns)/totalCpG, sum(cpg.in.intergenic)/totalCpG),
                      "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

neg.genic.results$Context <- factor(neg.genic.results$Context, levels=c("promoters", "exons", "introns", "intergenic"))


enrichment.genic <- ggplot(data = neg.genic.results, mapping = aes(Context, Proportions)) + 
                               geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
                               scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
                               labs(title = "Enrichment of CpGs with (-) Age-Association") +
                               ylab("Percentage of Total CpGs") +
                               xlab("") +
                               theme(axis.line = element_line(colour = "black")) +
                               font("xlab", size = 18, color = "black") +
                               font("ylab", size = 18, color = "black") +
                               font("xy.text", size = 18, color = "black") +
                               font("title", size = 15, color = "black", face = "bold")

enrichment.genic


## Statistical test:

# Promoters
binom.test(x = sum(cpg.in.promo), n = totalCpG, p = ( sum(sites.in.promo)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 7.854e-13

#Exons
binom.test(x = sum(cpg.in.exons), n = totalCpG, p = ( sum(sites.in.exons)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.003657

#Introns
binom.test(x = sum(cpg.in.introns), n = totalCpG, p = ( sum(sites.in.introns)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 7.397e-07

#Intergenic
binom.test(x = sum(cpg.in.intergenic), n = totalCpG, p = ( sum(sites.in.intergenic)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.002169


## Now plot all of these results together


pos.density.results$sign <- paste0("+ ", pos.density.results$type)
neg.density.results$sign <- paste0("- ", pos.density.results$type)

total.density <- rbind(pos.density.results, neg.density.results)
total.density$sign[c(1,2,3,4,9,10,11,12)] <- "Background"
total.density$sign <- factor(total.density$sign, levels = c("- Age-Associated", "Background", "+ Age-Associated"))

# Plot and save fig 1F

f <- ggplot(data = total.density, aes(x = Context, y = Proportions)) + 
    geom_bar(aes(fill = sign), stat = "identity", position = "dodge", color = "black") + 
    scale_fill_manual(values=c("navy", "forestgreen", "firebrick")) +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text=element_text(size=10),
          axis.title=element_text(size=12,face="bold"))
f

total.density[,c(1,2,5)]

ggsave(f, filename = "~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/images/density_combo.svg", device = "svg",
       width = 3, height = 3)

## Get counts to manually add to plot

pos.genic.results$sign <- paste0("+ ", pos.genic.results$type)
neg.genic.results$sign <- paste0("- ", pos.genic.results$type)

total.genic <- rbind(pos.genic.results, neg.genic.results)
total.genic$sign[c(1,2,3,4,9,10,11,12)] <- "Background"
total.genic$sign <- factor(total.genic$sign, levels = c("- Age-Associated", "Background", "+ Age-Associated"))

# plot and save fig 1G

g <- ggplot(data = total.genic, aes(x = Context, y = Proportions)) + 
    geom_bar(aes(fill = sign), stat = "identity", position = "dodge",   color = "black") + 
    scale_fill_manual(values=c("navy", "forestgreen", "firebrick")) +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text=element_text(size=10),
          axis.title=element_text(size=12,face="bold"))
    
g

ggsave(g, filename = "~/Parrott_Lab/AnoleAge/Fig1.1-DSSLogAge/images/genic_combo.svg", device = "svg",
       width = 3, height = 3)

total.genic[,c(1,2,5)]

#### Repeat entire process with WGCNA blue module sites
library(GenomicRanges)

grToBed <- function(gr, filename){
    df <- as.data.frame(gr)
    df.formatted <- df[,c("seqnames", "start", "end")]
    write.table(df.formatted, file = filename, sep = "\t", row.names = F, col.names = F, quote = F)
}


islands <- import("~/Parrott_Lab/Anole_Data/AnoSag2.1/CGI.gff", format = "gff")

shores.up <- trim(flank(islands, 2000))
shores.down <- trim(flank(islands, 2000, start = FALSE))
shores <- c(shores.up, shores.down)
shores.clean <- reduce(GenomicRanges::setdiff(shores, islands))

shelves.up <- trim(flank(shores.up, 2000))
shelves.down <- trim(flank(shores.down, 2000, start = FALSE))
shelves <- c(shelves.up, shelves.down)
shelves.no.islands <- GenomicRanges::setdiff(shelves, islands)
shelves.clean <- reduce(GenomicRanges::setdiff(shelves.no.islands, shores))


not.seas <- c(islands,shores,shelves)
seas <- gaps(not.seas)

#### CGdensity enrichment ------------------------------------------------------------------------------------

## Find background likelihood that CpGs are found in each context

CGsites <- fread("~/Parrott_Lab/AnoleAge/VariableSites.meth", sep = "\t", header = TRUE)
library(stringr)


CGsites <- data.frame("chr" = str_split_fixed(colnames(CGsites)[-c(1,2,3)], ":", 2)[,1], 
                      "start" = as.numeric(str_split_fixed(colnames(CGsites)[-c(1,2,3)], ":", 2)[,2]),
                      "end" = as.numeric(str_split_fixed(colnames(CGsites)[-c(1,2,3)], ":", 2)[,2]))

CGsites.clean <- GRanges(seqnames = CGsites$chr, 
                         ranges = IRanges(start = CGsites$start, end = CGsites$end),
                         strand = "*")

totalCG <- length(CGsites.clean)


sites.in.islands <- countOverlaps(CGsites.clean, islands, type = "within")
sum(sites.in.islands) # 638
cat("Proportion of CG sites that fall in islands: ", sum(sites.in.islands) / length(CGsites.clean), "\n")


sites.in.shores <- countOverlaps(CGsites.clean, shores.clean, type = "within")
sum(sites.in.shores) # 4964
cat("Proportion of CG sites that fall in shores: ", sum(sites.in.shores) / length(CGsites.clean), "\n")


sites.in.shelves <- countOverlaps(CGsites.clean, shelves.clean, type = "within")
sum(sites.in.shelves) # 3701
cat("Proportion of CG sites that fall in shelves: ", sum(sites.in.shelves) / length(CGsites.clean), "\n")


sites.in.seas <- countOverlaps(CGsites.clean, seas, type = "within")
sum(sites.in.seas) # 33346
cat("Proportion of CG sites that fall in seas: ", sum(sites.in.seas) / length(CGsites.clean), "\n")


CPG <- import("~/Parrott_Lab/AnoleAge/WGCNA/data/blue.bedGraph", format = "bed")
CPG <- CPG[which(CPG$name > 0),]

totalCpG <- length(CPG)

cpg.in.islands <- countOverlaps(CPG, islands, type = "within")
sum(cpg.in.islands) # 14
cat("Proportion of (+) Age-Associated CpGs that fall in islands: ", sum(cpg.in.islands) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, islands, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/POSAgeAssoc_cpg_islands.bed")


cpg.in.shores <- countOverlaps(CPG, shores.clean, type = "within")
sum(cpg.in.shores) # 120
cat("Proportion of (+) Age-Associated CpGs that fall in shores: ", sum(cpg.in.shores) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, shores.clean, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/POSAgeAssoc_cpg_shores.bed")


cpg.in.shelves <- countOverlaps(CPG, shelves.clean, type = "within")
sum(cpg.in.shelves) # 78
cat("Proportion of (+) Age-Associated CpGs that fall in shelves: ", sum(cpg.in.shelves) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, shelves.clean, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/POSAgeAssoc_cpg_shelves.bed")


cpg.in.seas <- countOverlaps(CPG, seas, type = "within")
sum(cpg.in.seas) # 744
cat("Proportion of (+) Age-Associated CpGs that fall in seas: ", sum(cpg.in.seas) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, seas, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/POSAgeAssoc_cpg_seas.bed")


pos.density.results <- data.frame("Context" = c("islands", "shores", "shelves", "seas","islands", "shores", "shelves", "seas"), 
                                  "Values" = c(sum(sites.in.islands), sum(sites.in.shores),
                                               sum(sites.in.shelves), sum(sites.in.seas),
                                               sum(cpg.in.islands), sum(cpg.in.shores), 
                                               sum(cpg.in.shelves), sum(cpg.in.seas)),
                                  "Proportions" = c(sum(sites.in.islands)/totalCG, sum(sites.in.shores)/totalCG,
                                                    sum(sites.in.shelves)/totalCG, sum(sites.in.seas)/totalCG,
                                                    sum(cpg.in.islands)/totalCpG, sum(cpg.in.shores)/totalCpG, 
                                                    sum(cpg.in.shelves)/totalCpG, sum(cpg.in.seas)/totalCpG),
                                  "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

pos.density.results$Context <- factor(pos.density.results$Context, levels=c("islands","shores","shelves", "seas"))

## Statistical test:

#Islands
binom.test(x = sum(cpg.in.islands), n = totalCpG, p = ( sum(sites.in.islands)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 1

#Shores
binom.test(x = sum(cpg.in.shores), n = totalCpG, p = ( sum(sites.in.shores)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 0.36

#Shelves
binom.test(x = sum(cpg.in.shelves), n = totalCpG, p = ( sum(sites.in.shelves)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.6

#Seas
binom.test(x = sum(cpg.in.seas), n = totalCpG, p = ( sum(sites.in.seas)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 0.81


#### Generate genic contexts ----------------------------------------------------------------------------------------------

exons <- import.gff2("~/Parrott_Lab/Anole_Data/AnoSag2.1/AnoSag2.1_exons.gtf")
exons.destranded <- reduce(exons, ignore.strand = TRUE)
#export(exons.destranded, con = "~/Parrott_Lab/Parrott_Lab/Anole_Data/AnoSag2.1/AnoSag2.1_exons_destranded.gtf", format = "gff2")


CDS <- import.gff2("~/Parrott_Lab/Anole_Data/AnoSag2.1/CDS.gtf")
CDS.destranded <- reduce(CDS, ignore.strand = TRUE)
full.CDS <- GenomicRanges::union(CDS.destranded, exons.destranded)
#export(full.CDS, con = "~/Parrott_Lab/Parrott_Lab/Anole_Data/AnoSag2.1/AnoSag2.1_CDS_destranded.gtf", format = "gff2")


intergenic <- gaps(CDS.destranded)
intergenic.clean <- GenomicRanges::setdiff(intergenic, exons.destranded)
#export(intergenic.clean, con = "~/Parrott_Lab/Parrott_Lab/Anole_Data/AnoSag2.1/AnoSag2.1_intergenic.gtf", format = "gff2")

exons.comp <- gaps(exons)

introns <- GenomicRanges::setdiff(CDS.destranded, exons.destranded)
#export(introns, con = "~/Parrott_Lab/Parrott_Lab/Anole_Data/AnoSag2.1/AnoSag2.1_introns.gtf", format = "gff2")

promoter <- reduce(promoters(CDS), ignore.strand = TRUE)


#### Genic context Enrichment  

## Find background likelihood that CpGs are found in each context

sites.in.promo <- countOverlaps(CGsites.clean, promoter, type = "within")
sum(sites.in.promo) # 1,077
cat("Proportion of CG sites that fall in promoters: ", sum(sites.in.promo) / length(CGsites.clean), "\n")


sites.in.exons <- countOverlaps(CGsites.clean, exons.destranded, type = "within")
sum(sites.in.exons) # 1744
cat("Proportion of CG sites that fall in exons: ", sum(sites.in.exons) / length(CGsites.clean), "\n")


sites.in.introns <- countOverlaps(CGsites.clean, introns, type = "within")
sum(sites.in.introns) # 16514
cat("Proportion of CG sites that fall in introns: ", sum(sites.in.introns) / length(CGsites.clean), "\n")


sites.in.intergenic <- countOverlaps(CGsites.clean, intergenic.clean, type = "within")
sum(sites.in.intergenic) # 24358
cat("Proportion of CG sites that fall in intergenic space: ", sum(sites.in.intergenic) / length(CGsites.clean), "\n")



cpg.in.promo <- countOverlaps(CPG, promoter, type = "within")
sum(cpg.in.promo) # 28
cat("Proportion of (+) Age-Associated CpGs that fall in promoters: ", sum(cpg.in.promo) / length(CPG), "\n")


cpg.in.exons <- countOverlaps(CPG, exons.destranded, type = "within")
sum(cpg.in.exons) # 61
cat("Proportion of (+) Age-Associated CpGs that fall in exons: ", sum(cpg.in.exons) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, exons.destranded, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/POSAgeAssoc_cpg_exons.bed")


cpg.in.introns <- countOverlaps(CPG, introns, type = "within")
sum(cpg.in.introns) # 384
cat("Proportion of (+) Age-Associated CpGs that fall in introns: ", sum(cpg.in.introns) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, introns, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/POSAgeAssoc_cpg_introns.bed")


cpg.in.intergenic <- countOverlaps(CPG, intergenic.clean, type = "within")
sum(cpg.in.intergenic) # 510
cat("Proportion of (+) Age-Associated CpGs that fall in intergenic space: ", sum(cpg.in.intergenic) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, intergenic.clean, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/POSAgeAssoc_cpg_intergenic.bed")


pos.genic.results <- data.frame("Context" = c("promoters", "exons", "introns", "intergenic", "promoters", "exons", "introns", "intergenic"), 
                                "Values" = c(sum(sites.in.promo), sum(sites.in.exons), 
                                             sum(sites.in.introns), sum(sites.in.intergenic),
                                             sum(cpg.in.promo), sum(cpg.in.exons), 
                                             sum(cpg.in.introns), sum(cpg.in.intergenic)),
                                "Proportions" = c(sum(sites.in.promo)/totalCG, sum(sites.in.exons)/totalCG, 
                                                  sum(sites.in.introns)/totalCG, sum(sites.in.intergenic)/totalCG,
                                                  sum(cpg.in.promo)/totalCpG, sum(cpg.in.exons)/totalCpG, 
                                                  sum(cpg.in.introns)/totalCpG, sum(cpg.in.intergenic)/totalCpG),
                                "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

pos.genic.results$Context <- factor(pos.genic.results$Context, levels=c("promoters", "exons", "introns", "intergenic"))

## Statistical test:

#Promoters
binom.test(x = sum(cpg.in.promo), n = totalCpG, p = ( sum(sites.in.promo)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.4085

#Exons
binom.test(x = sum(cpg.in.exons), n = totalCpG, p = ( sum(sites.in.exons)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.00076

#Introns
binom.test(x = sum(cpg.in.introns), n = totalCpG, p = ( sum(sites.in.introns)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 0.353

#Intergenic
binom.test(x = sum(cpg.in.intergenic), n = totalCpG, p = ( sum(sites.in.intergenic)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.02

## Negative


CPG <- import("~/Parrott_Lab/AnoleAge/WGCNA/data/blue.bedGraph", format = "bed")
CPG <- CPG[which(CPG$name < 0),]
totalCpG <- length(CPG)

cpg.in.islands <- countOverlaps(CPG, islands, type = "within")
sum(cpg.in.islands) # 22
cat("Proportion of (-) Age-Associated CpGs that fall in islands: ", sum(cpg.in.islands) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, islands, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/NEGAgeAssoc_cpg_islands.bed")


cpg.in.shores <- countOverlaps(CPG, shores.clean, type = "within")
sum(cpg.in.shores) # 188
cat("Proportion of (-) Age-Associated CpGs that fall in shores: ", sum(cpg.in.shores) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, shores.clean, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/NEGAgeAssoc_cpg_shores.bed")


cpg.in.shelves <- countOverlaps(CPG, shelves.clean, type = "within")
sum(cpg.in.shelves) # 124
cat("Proportion of (-) Age-Associated CpGs that fall in shelves: ", sum(cpg.in.shelves) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, shelves.clean, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/NEGAgeAssoc_cpg_shelves.bed")


cpg.in.seas <- countOverlaps(CPG, seas, type = "within")
sum(cpg.in.seas) # 722
cat("Proportion of (-) Age-Associated CpGs that fall in seas: ", sum(cpg.in.seas) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, seas, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/NEGAgeAssoc_cpg_seas.bed")


neg.density.results <- data.frame("Context" = c("islands", "shores", "shelves", "seas","islands", "shores", "shelves", "seas"), 
                                  "Values" = c(sum(sites.in.islands), sum(sites.in.shores),
                                               sum(sites.in.shelves), sum(sites.in.seas),
                                               sum(cpg.in.islands), sum(cpg.in.shores), 
                                               sum(cpg.in.shelves), sum(cpg.in.seas)),
                                  "Proportions" = c(sum(sites.in.islands)/totalCG, sum(sites.in.shores)/totalCG,
                                                    sum(sites.in.shelves)/totalCG, sum(sites.in.seas)/totalCG,
                                                    sum(cpg.in.islands)/totalCpG, sum(cpg.in.shores)/totalCpG, 
                                                    sum(cpg.in.shelves)/totalCpG, sum(cpg.in.seas)/totalCpG),
                                  "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

neg.density.results$Context <- factor(neg.density.results$Context, levels=c("islands","shores","shelves", "seas"))

## Statistical test:

#Islands
binom.test(x = sum(cpg.in.islands), n = totalCpG, p = ( sum(sites.in.islands)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.12

#Shores
binom.test(x = sum(cpg.in.shores), n = totalCpG, p = ( sum(sites.in.shores)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 4.017e-09
#Shelves
binom.test(x = sum(cpg.in.shelves), n = totalCpG, p = ( sum(sites.in.shelves)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.0006802

#Seas
binom.test(x = sum(cpg.in.seas), n = totalCpG, p = ( sum(sites.in.seas)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 2.128e-13

#### Now genic

cpg.in.promo <- countOverlaps(CPG, promoter, type = "within")
sum(cpg.in.promo) # 31
cat("Proportion of (-) Age-Associated CpGs that fall in promoters: ", sum(cpg.in.promo) / length(CPG), "\n")


cpg.in.exons <- countOverlaps(CPG, exons.destranded, type = "within")
sum(cpg.in.exons) # 46
cat("Proportion of (-) Age-Associated CpGs that fall in exons: ", sum(cpg.in.exons) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, exons.destranded, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/NEGAgeAssoc_cpg_exons.bed")


cpg.in.introns <- countOverlaps(CPG, introns, type = "within")
sum(cpg.in.introns) # 393
cat("Proportion of (-) Age-Associated CpGs that fall in introns: ", sum(cpg.in.introns) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, introns, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/NEGAgeAssoc_cpg_introns.bed")


cpg.in.intergenic <- countOverlaps(CPG, intergenic.clean, type = "within")
sum(cpg.in.intergenic) # 615
cat("Proportion of (-) Age-Associated CpGs that fall in intergenic space: ", sum(cpg.in.intergenic) / length(CPG), "\n")
#grToBed(subsetByOverlaps(CPG, intergenic.clean, type = "within"), "~/Parrott_Lab/Anole/AnoleAge/Fig1/NEGAgeAssoc_cpg_intergenic.bed")


neg.genic.results <- data.frame("Context" = c("promoters", "exons", "introns", "intergenic", "promoters", "exons", "introns", "intergenic"), 
                                "Values" = c(sum(sites.in.promo), sum(sites.in.exons), 
                                             sum(sites.in.introns), sum(sites.in.intergenic),
                                             sum(cpg.in.promo), sum(cpg.in.exons), 
                                             sum(cpg.in.introns), sum(cpg.in.intergenic)),
                                "Proportions" = c(sum(sites.in.promo)/totalCG, sum(sites.in.exons)/totalCG, 
                                                  sum(sites.in.introns)/totalCG, sum(sites.in.intergenic)/totalCG,
                                                  sum(cpg.in.promo)/totalCpG, sum(cpg.in.exons)/totalCpG, 
                                                  sum(cpg.in.introns)/totalCpG, sum(cpg.in.intergenic)/totalCpG),
                                "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

neg.genic.results$Context <- factor(neg.genic.results$Context, levels=c("promoters", "exons", "introns", "intergenic"))

## Statistical test:

# Promoters
binom.test(x = sum(cpg.in.promo), n = totalCpG, p = ( sum(sites.in.promo)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.3766

#Exons
binom.test(x = sum(cpg.in.exons), n = totalCpG, p = ( sum(sites.in.exons)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.64

#Introns
binom.test(x = sum(cpg.in.introns), n = totalCpG, p = ( sum(sites.in.introns)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 0.32

#Intergenic
binom.test(x = sum(cpg.in.intergenic), n = totalCpG, p = ( sum(sites.in.intergenic)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.45


## Combo plots


pos.density.results$sign <- paste0("+ ", pos.density.results$type)
neg.density.results$sign <- paste0("- ", pos.density.results$type)

total.density <- rbind(pos.density.results, neg.density.results)
total.density$sign[c(1,2,3,4,9,10,11,12)] <- "Background"
total.density$sign <- factor(total.density$sign, levels = c("- Age-Associated", "Background", "+ Age-Associated"))

h <- ggplot(data = total.density, aes(x = Context, y = Proportions)) + 
    geom_bar(aes(fill = sign), stat = "identity", position = "dodge",   color = "black") + 
    scale_fill_manual(values=c("navy", "forestgreen", "firebrick")) +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text=element_text(size=10),
          axis.title=element_text(size=12,face="bold"))
h

total.density[,c(1,2,5)]

ggsave(h, filename = "~/Parrott_Lab/AnoleAge/WGCNA/images/blueMod_density_combo.svg", device = "svg",
       width = 3, height = 3)

pos.genic.results$sign <- paste0("+ ", pos.genic.results$type)
neg.genic.results$sign <- paste0("- ", pos.genic.results$type)

total.genic <- rbind(pos.genic.results, neg.genic.results)
total.genic$sign[c(1,2,3,4,9,10,11,12)] <- "Background"
total.genic$sign <- factor(total.genic$sign, levels = c("- Age-Associated", "Background", "+ Age-Associated"))

i <- ggplot(data = total.genic, aes(x = Context, y = Proportions)) + 
    geom_bar(aes(fill = sign), stat = "identity", position = "dodge",   color = "black") + 
    scale_fill_manual(values=c("navy", "forestgreen", "firebrick")) +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text=element_text(size=10),
          axis.title=element_text(size=12,face="bold"))
i

i.leg <- ggplot(data = total.genic, aes(x = Context, y = Proportions)) + 
  geom_bar(aes(fill = sign), stat = "identity", position = "dodge",   color = "black") + 
  scale_fill_manual(values=c("navy", "forestgreen", "firebrick")) +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"))

leg <- cowplot::get_legend(i.leg)

grid.newpage()

leg.plot <- grid.draw(leg)

ggsave(leg, filename = "~/Parrott_Lab/AnoleAge/WGCNA/images/enrichments_legend.svg", device = "svg",
       width = 3, height = 2)

total.genic[,c(1,2,5)]

ggsave(i, filename = "~/Parrott_Lab/AnoleAge/WGCNA/images/blueMod_genic_combo.svg", device = "svg",
       width = 2, height = 1)

