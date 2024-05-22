library(data.table)
library(GenomicRanges)
library(ggpubr)

theme_set(theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")))

## import unfiltered methylation file and manually create metadata object

meth_file <- as.data.frame(fread(file = "E:/anole/CpG/processed/labeled_allSites_stranded.tab", header = TRUE, sep = "\t"))

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F")))

meth.gr <- GRanges(seqnames = meth_file$chr, 
                   ranges = IRanges(start = meth_file$start, end = meth_file$end),
                   strand = "*", mcols = meth_file[,c(5:41)] )

## read in age DMC coordinates

ageSites <- read.table("C:/Users/sheal/Desktop/AnoleAge/Fig1.1-DSSLogAge/data/LogAgeSites_smoothed200.bed", header = FALSE)

ageSites <- ageSites[-4671,] ## For some reason this single site isn't represented in full dataset - remove it

ageSites.gr <- GRanges(ageSites$V1, IRanges(ageSites$V2, ageSites$V3), strand = "*")

## overlap age sites with total dataset to get meth values at ageDMCs

meth.ageSites <- subsetByOverlaps(meth.gr, ageSites.gr, type = "end")

ageSites.mat <- as.data.frame(t(as.data.frame(meth.ageSites)[,-c(1,2,3,4,5)]))
colnames(ageSites.mat) <- paste0(ageSites$V1, ".", ageSites$V2)

ageSites.df.clean <- data.frame(annotations, ageSites.mat)

## Split into male and female, and calculate number of sites missing in 2/3 of each group

ageSites.m <- ageSites.df.clean[which(ageSites.df.clean$sex == "M"),]
missing.m <- which(colSums(is.na(ageSites.m)) > round(length(ageSites.m$ID)*2/3, 0))

ageSites.f <- ageSites.df.clean[which(ageSites.df.clean$sex == "F"),]
missing.f <- which(colSums(is.na(ageSites.f)) > round(length(ageSites.m$ID)*2/3, 0))

## take any sites missing at >66% of either males or females and remove these

missing <- union(missing.m, missing.f)

ageSites.m <- ageSites.m[,-missing]
ageSites.f <- ageSites.f[,-missing]
ageSites.all <- ageSites.df.clean[,-missing]

# create empty vectors for rates

m.rate <- NA
f.rate <- NA
intTerm <- NA

## Loop through sites, fitting male and female age vs meth models, before extracting beta values for "rates"

for (i in 4:length(ageSites.m)){
  
    m.model <- lm(ageSites.m[,i] ~ ageSites.m$ages)
    m.rate <- tryCatch(append(m.rate, summary(m.model)$coefficients[2,1]), error = function(e) {
                       append(m.rate, NA) ## Try to get coefficients, add NA instead if model fails to fit
    })
    f.model <- lm(ageSites.f[,i] ~ ageSites.f$ages)
    f.rate <- tryCatch(append(f.rate, summary(f.model)$coefficients[2,1]), error = function(e) {
                       append(f.rate, NA)
    })
    model <- lm(ageSites.all[,i] ~ ageSites.all$sex*ageSites.all$ages)
    intTerm <- tryCatch(append(intTerm, summary(model)$coefficients[4,1]), error = function(e) {
                       append(intTerm, NA)
    })
}

# create dataframe with rates

modelSlopes <- data.frame("Site" = colnames(ageSites.m)[3:length(ageSites.m)],
                          "MaleRate" = m.rate,
                          "FemaleRate" = f.rate,
                          "SlopeDiff" = intTerm)[-1,]

## remove any sites which have different directionality in males and females

modelSlopes <- na.omit(modelSlopes)

cleanSites <- c(which(modelSlopes$MaleRate > 0 & modelSlopes$FemaleRate > 0),
                  which(modelSlopes$MaleRate < 0 & modelSlopes$FemaleRate < 0))

#### SITES WHICH HAVE OPPOSITE DIRECTIONALITY IN MALES / FEMALES ARE EXCLUDED (n = 656)

modelSlopes.clean <- modelSlopes[cleanSites,]

modelSlopes.clean$rateDiff <- modelSlopes.clean$MaleRate - modelSlopes.clean$FemaleRate

median(modelSlopes.clean$rateDiff, na.rm = TRUE) # Median rate difference = 0.11

## Split by directionality of change with age

modelSlopes.clean$Dir <- ifelse(modelSlopes.clean$MaleRate > 0 & modelSlopes.clean$FemaleRate > 0, "Increasing", "Decreasing")

modelSlopes.clean.long <- tidyr::pivot_longer(modelSlopes.clean, cols = c(2,3), values_to = "slope")

modelSlopes.clean.long <- as.data.frame(modelSlopes.clean.long)

# Plot rates comparing males and females at age DMCs

a <- ggplot(data = modelSlopes.clean.long, aes(x = name, y = slope, color = name)) + geom_boxplot() +
       xlab("Sex") + ylab("% Methylation per Month") + labs(title = "aDMCs") +
       theme(legend.position = "none", 
             strip.text = element_text(size = 12, face = "bold")) + facet_wrap("Dir")# + 
       #scale_x_discrete(labels = c("Female", "Male"))
a

ggsave(a, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/LogAgeSites_SexRates_ByDir.png",
       width = 4, height = 5)

# t-tests comparing rates seperately at sites which increase and decrease

t.test(modelSlopes.clean$FemaleRate[which(modelSlopes.clean$Dir == "Increasing")], 
       modelSlopes.clean$MaleRate[which(modelSlopes.clean$Dir == "Increasing")], paired = TRUE)
# est = 0.336 +/-0.031, p , 2.2e-16

t.test(modelSlopes.clean$FemaleRate[which(modelSlopes.clean$Dir == "Decreasing")], 
       modelSlopes.clean$MaleRate[which(modelSlopes.clean$Dir == "Decreasing")], paired = TRUE)
# est = -0.279 +/-0.020, p , 2.2e-16

#### Plot raw values

library(tidyverse)

up <- ageSites.all[,-c(1,2,3)][,cleanSites][,which(modelSlopes.clean$Dir == "Increasing")]

up <- data.frame(annotations, up)
up$ages <- factor(up$ages, levels=c(1, 7, 12, 18, 24, 30, 36, 40, 48, 54, 60))

down <- ageSites.all[,-c(1,2,3)][,cleanSites][,which(modelSlopes.clean$Dir == "Decreasing")]

down <- data.frame(annotations, down)
down$ages <- factor(down$ages, levels=c(1, 7, 12, 18, 24, 30, 36, 40, 48, 54, 60))

up.long <- pivot_longer(up, cols = c(4:length(up)), names_to = "Site", values_to = "Meth")

down.long <- pivot_longer(down, cols = c(4:length(down)), names_to = "Site", values_to = "Meth")

b1 <- ggplot(down.long, aes(x = ages, y = Meth, color = sex)) + 
    geom_jitter(height = 0, width = 0.1, color = "grey", alpha = 0.1) + 
        geom_boxplot(outlier.alpha = 0, alpha = 0.8) +
    geom_smooth(aes(group = sex), method = "lm", formula = y ~ log(x)) + xlab("Age in Months") + 
    scale_x_discrete(drop = FALSE, breaks = c(1, 7, 18, 36, 48, 60))

b2 <- ggplot(up.long, aes(x = ages, y = Meth, color = sex)) + 
    geom_jitter(height = 0, width = 0.1, color = "grey", alpha = 0.1) + 
    geom_boxplot(outlier.alpha = 0, alpha = 0.8) +
    geom_smooth(aes(group = sex), method = "lm", formula = y ~ log(x)) + xlab("Age in Months") + 
    scale_x_discrete(drop = FALSE, breaks = c(1, 7, 18, 36, 48, 60))

bFull <- ggarrange(b1, b2, ncol = 1, nrow = 2)
bFull

ggsave(bFull, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/LogAgeSites_MethvsAge.png",
       width = 6, height = 8)

## Looks like counterintuitive pattern could be driven by more "robust" initial methylation in females - look at early life differences

## limit to just earliest two age groups

young.df <- ageSites.df.clean[which(ageSites.df.clean$ages %in% c(1,7)),]

m.avgs <- colMeans(young.df[which(young.df$sex == "M"), -c(1,2,3)], na.rm = TRUE)
f.avgs <- colMeans(young.df[which(young.df$sex == "F"), -c(1,2,3)], na.rm = TRUE)

youngAvgs <- data.frame("Site" = colnames(young.df)[-c(1,2,3)],
                        "YoungMaleAvg" = m.avgs,
                        "YoungFemaleAvg" = f.avgs)

# again remove same sites which were omitted from rate comparison

youngAvgs.clean <- youngAvgs[-missing,]
youngAvgs.clean <- youngAvgs.clean[cleanSites,]

youngAvgs.long <- tidyr::pivot_longer(youngAvgs, cols = c(2,3), names_to = "Sex", values_to = "Avgs")
#youngAvgs.long2 <- tidyr::pivot_longer(youngAvgs.long, cols = )

mergedSites <- merge(modelSlopes.clean, youngAvgs.clean, by = "Site")

mergedSites$chr <- as.numeric(stringr::str_remove_all(mergedSites$chr, "scaffold_"))

mergedSites.long <- tidyr::pivot_longer(mergedSites, cols = c(9,10), names_to = "Sex", values_to = "Avgs")

# plot early life averages

b <- ggplot(mergedSites.long, aes(x = Sex, y = Avgs, color = Sex)) + geom_boxplot() + 
    facet_wrap("Dir") + theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold")) + ylab("Avg % Methylation") + 
    labs(title = "aDMCs in Early life (1mo)") + 
    scale_x_discrete(labels = c("Female", "Male"))
b

ggsave(b, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/LogAgeSites_EarlyLifeAvgs_ByDir.png",
       width = 4, height = 5)

## t-tests of early life averages

t.test(mergedSites$YoungFemaleAvg[which(mergedSites$Dir == "Increasing")], 
       mergedSites$YoungMaleAvg[which(mergedSites$Dir == "Increasing")], paired = TRUE)
#est: -12.07, 95% CI: -11.67 - -7.67; p < 2.2e-16

t.test(mergedSites$YoungFemaleAvg[which(mergedSites$Dir == "Decreasing")], 
       mergedSites$YoungMaleAvg[which(mergedSites$Dir == "Decreasing")], paired = TRUE)
#est: 9.72, 95% CI: 8.39 - 11.03; p = 4.1e-13



## Now repeat for blue module sites --------------------------------------------------------------------------------------------------------

bluemod.mat <- read.table("C:/Users/sheal/Desktop/AnoleAge/WGCNA/data/blue.mat")

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F")))

bluemod.clean <- data.frame(annotations[c(-4,-37),], bluemod.mat)

## Take only individuals above age at maturity:
#bluemod.clean <- bluemod.clean[which(bluemod.clean$ages > 12),]

bluemod.m <- bluemod.clean[which(bluemod.clean$sex == "M"),]
bluemod.f <- bluemod.clean[which(bluemod.clean$sex == "F"),]


## Declare empty variables
m.rate <- NA
f.rate <- NA

for (i in 4:length(bluemod.m)){
    
        m.model <- lm(bluemod.m[,i] ~ bluemod.m$ages)
    m.rate <- tryCatch(append(m.rate, summary(m.model)$coefficients[2,1]), error = function(e) {
                       append(m.rate, NA) ## Try to get coefficients, add NA instead if model fails to fit
    })
    f.model <- lm(bluemod.f[,i] ~ bluemod.f$ages)
    f.rate <- tryCatch(append(f.rate, summary(f.model)$coefficients[2,1]), error = function(e) {
                       append(f.rate, NA)
    })
}

modelSlopes.wgcna <- data.frame("Site" = colnames(bluemod.m)[3:length(bluemod.m)],
                                "MaleRate" = m.rate,
                                "FemaleRate" = f.rate)[-1,]

modelSlopes.wgcna2 <- data.frame("Site" = colnames(bluemod.m)[3:length(bluemod.m)],
                                "MaleRate" = m.rate,
                                "FemaleRate" = f.rate)[-1,]


cleanSites <- c(which(modelSlopes.wgcna$MaleRate > 0 & modelSlopes.wgcna$FemaleRate > 0),
                  which(modelSlopes.wgcna$MaleRate < 0 & modelSlopes.wgcna$FemaleRate < 0))

modelSlopes.wgcna <- modelSlopes.wgcna[cleanSites,]

modelSlopes.wgcna$Dir <- ifelse(modelSlopes.wgcna$MaleRate > 0, "Increasing", "Decreasing")

modelSlopes.wgcna.long <- tidyr::pivot_longer(modelSlopes.wgcna, cols = c(2,3), values_to = "slope")


d <- ggplot(data = modelSlopes.wgcna.long, aes(x = name, y = slope, color = name)) + geom_boxplot() +
       xlab("Sex") + ylab("% Methylation per Month") + labs(title = "Blue Module") +
       theme(legend.position = "none", 
             strip.text = element_text(size = 12, face = "bold")) + facet_wrap("Dir") + 
       scale_x_discrete(labels = c("Female", "Male"))
d 

ggsave(d, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/WGCNA_SexRates_ByDir.png",
       width = 4, height = 5)

## Now look at differences in methylation at young ages in these individuals

young.wgcna.df <- bluemod.clean[which(ageSites.df.clean$ages < 12),]

m.wgcna.avgs <- colMeans(young.wgcna.df[which(young.wgcna.df$sex == "M"), -c(1,2,3)], na.rm = TRUE)
f.wgcna.avgs <- colMeans(young.wgcna.df[which(young.wgcna.df$sex == "F"), -c(1,2,3)], na.rm = TRUE)

youngAvgs.wgcna <- data.frame("Site" = colnames(young.wgcna.df)[-c(1,2,3)],
                        "YoungMaleAvg" = m.wgcna.avgs,
                        "YoungFemaleAvg" = f.wgcna.avgs)


youngAvgs.wgcna.long <- tidyr::pivot_longer(youngAvgs.wgcna, cols = c(2,3), names_to = "Sex", values_to = "Avgs")


mean(m.wgcna.avgs, na.rm = TRUE); mean(f.wgcna.avgs, na.rm = TRUE)

mergedSites.wgcna <- merge(modelSlopes.wgcna, youngAvgs.wgcna, by = "Site")

mergedSites.wgcna$Dir <- ifelse(mergedSites.wgcna$MaleRate > 0, "Increasing", "Decreasing")

## Paired T-tests

t.test(mergedSites.wgcna$FemaleRate[which(mergedSites.wgcna$Dir == "Decreasing")], 
       mergedSites.wgcna$MaleRate[which(mergedSites.wgcna$Dir == "Decreasing")], paired = TRUE)
#est: -0.091, 95% CI: -0.130 - -0.0519; p = 5.994e-06

t.test(mergedSites.wgcna$FemaleRate[which(mergedSites.wgcna$Dir == "Increasing")], 
       mergedSites.wgcna$MaleRate[which(mergedSites.wgcna$Dir == "Increasing")], paired = TRUE)
#est: 0.04 , 95% CI: -0.001  0.083; p = 0.06065


mergedSites.wgcna.long <- tidyr::pivot_longer(mergedSites.wgcna, cols = c(5,6), names_to = "Sex", values_to = "Avgs")

e <- ggplot(mergedSites.wgcna.long, aes(x = Sex, y = Avgs, color = Sex)) + geom_boxplot() + 
    facet_wrap("Dir") + theme(legend.position = "none", 
             strip.text = element_text(size = 12, face = "bold")) + 
    scale_x_discrete(labels = c("Female", "Male")) +
    ylab("Avg Methylation") + labs(title = "Blue Module in Early life (1-7mo)")
e

ggsave(e, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/WGCNA_EarlyLifeAvgs_ByDir.png",
       width = 4, height = 5)

mergedSites.wgcna.long2 <- tidyr::pivot_longer(mergedSites.wgcna, cols = c(2,3), names_to = "Sex", values_to = "Rate")

t.test(mergedSites.wgcna$YoungFemaleAvg[which(mergedSites.wgcna$Dir == "Decreasing")], 
       mergedSites.wgcna$YoungMaleAvg[which(mergedSites.wgcna$Dir == "Decreasing")], paired = TRUE)
#est: -1.09, 95% CI: -11.67 - -7.67; p < 0.22

t.test(mergedSites.wgcna$YoungFemaleAvg[which(mergedSites.wgcna$Dir == "Increasing")], 
       mergedSites.wgcna$YoungMaleAvg[which(mergedSites.wgcna$Dir == "Increasing")], paired = TRUE)
#est: 6.26, 95% CI: 4.40 - 8.13; p = 9.479e-11




#### Plot raw values for blue module

up <- bluemod.clean[,-c(1,2,3)][,cleanSites][,which(modelSlopes.wgcna$Dir == "Increasing")]

up <- data.frame(annotations[c(-4,-37),], up)
up$ages <- factor(up$ages, levels=c(1, 7, 12, 18, 24, 30, 36, 40, 48, 54, 60))

down <- bluemod.clean[,-c(1,2,3)][,cleanSites][,which(modelSlopes.wgcna$Dir == "Decreasing")]

down <- data.frame(annotations[c(-4,-37),], down)
down$ages <- factor(down$ages, levels=c(1, 7, 12, 18, 24, 30, 36, 40, 48, 54, 60))

up.long <- pivot_longer(up, cols = c(4:length(up)), names_to = "Site", values_to = "Meth")

down.long <- pivot_longer(down, cols = c(4:length(down)), names_to = "Site", values_to = "Meth")

b1 <- ggplot(down.long, aes(x = ages, y = Meth, color = sex)) + 
    geom_jitter(height = 0, width = 0.1, color = "grey", alpha = 0.1) + 
        geom_boxplot(outlier.alpha = 0, alpha = 0.8) +
    geom_smooth(aes(group = sex), method = "lm", formula = y ~ log(x)) + xlab("Age in Months") + 
    scale_x_discrete(drop = FALSE, breaks = c(1, 7, 18, 36, 48, 60))

b2 <- ggplot(up.long, aes(x = ages, y = Meth, color = sex)) + 
    geom_jitter(height = 0, width = 0.1, color = "grey", alpha = 0.1) + 
    geom_boxplot(outlier.alpha = 0, alpha = 0.8) +
    geom_smooth(aes(group = sex), method = "lm", formula = y ~ log(x)) + xlab("Age in Months") + 
    scale_x_discrete(drop = FALSE, breaks = c(1, 7, 18, 36, 48, 60))

bFull <- ggarrange(b1, b2, ncol = 1, nrow = 2)
bFull

ggsave(bFull, filename = "C:/Users/sheal/Desktop/AnoleAge/Fig5-SexDiff/images/BlueModSites_MethvsAge.png",
       width = 6, height = 8)
