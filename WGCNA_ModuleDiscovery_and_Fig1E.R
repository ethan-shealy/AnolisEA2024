## WGCNA #########################################################################################################################################
library(WGCNA)
library(data.table)

options(stringsAsFactors = FALSE);

setwd("C:/Users/sheal/Desktop/AnoleAge")

## manually create annotations/metadata

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M",
                                              "F","F","F","F","M","M","M","M","M","M","F","M","F","F")),
                          "sexNum" = c(0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,1,1,0,1,0,1,1,
                                              0,0,0,0,1,1,1,1,1,1,0,1,0,0))

# remove low-coverage individuals which could not reliably be imputed

annotations <- annotations[c(-4,-37),]

## import methylation matrix which has been knn-imputed, and all sites with sd < 0.3 removed. 

meth <- fread("VariableSites.meth", sep = "\t", header = TRUE)

rownames(annotations) <- annotations$ID

mat <- as.matrix(meth[,c(-1,-2,-3)])

gsc <-  goodSamplesGenes(mat, verbose = 3)
gsc$allOK

if (!gsc$allOK) {
# Optionally, print the gene and sample names that were removed:
if (sum(!gsc$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(mat.stand)[!gsc$goodGenes], collapse = ", ")));
if (sum(!gsc$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(mat.stand)[!gsc$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
mat.stand <- mat.stand[gsc$goodSamples, gsc$goodGenes]
}


sampleTree = hclust(dist(mat), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


datExpr <- mat
rownames(datExpr) <- annotations$ID
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr);
traitRows = match(Samples, annotations$ID);
datTraits = annotations[traitRows, -1];
rownames(datTraits) = annotations[traitRows, 1];
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits[c(1,3)], signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap",
addGuide = TRUE)


save(datExpr, datTraits, file = "C:/Users/sheal/Desktop/AnoleAge/Misc/WGCNA/VarData_anoleMeth-dataInput.RData")

rm(list=ls())

#### Next Steps ##########################################

setwd("C:/Users/sheal/Desktop/AnoleAge/Misc/WGCNA")

library(WGCNA)

enableWGCNAThreads(nThreads = 4)

lnames = load(file = "VarData_anoleMeth-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2])
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## Now run the actual module construction function

net = blockwiseModules(datExpr, power = 6, maxBlockSize = 22000,
TOMType = "signed", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "anoleAgeTOM",
verbose = 3)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "AnoleMeth-networkConstruction-auto.RData")

rm(list=ls())

###### The following script was used to generate Fig 1E and it's inlaid plot -------------------------------------------------------------

library(WGCNA)
library(data.table)
library(ggpubr)

## Set fig aesthetics

theme_set(theme_bw() + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")))


setwd("~/Parrott_Lab/AnoleAge/WGCNA")
# Load the expression and trait data saved in the first part
lnames = load(file = "data/VarData_anoleMeth-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load
lnames = load(file = "data/AnoleMeth-networkConstruction-auto.RData");
#The variable lnames contains the names of loaded variables.
lnames

annotations <- data.frame("ID" = c("S10","S15","S19","S23","S28","S33","S41","S46","S49","S51","S53","S61","S62","S64","S66","S67",
                                   "S70","S76","S84","S85","S95","S96", "S107","S110","S111","S116","S117","S120","S122","S125","S126",
                                   "S127","S128","S129","S131","S132","S133"),
                          "ages" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60),
                          "sex" = as.factor(c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M","F","M","F","M","M",
                                              "F","F","F","F","M","M","M","M","M","M","F","M","F","F")),
                          "sexNum" = c(0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,1,1,0,1,0,1,1,
                                              0,0,0,0,1,1,1,1,1,1,0,1,0,0))

annotations <- annotations[c(-4,-37),]

meth <- fread("data/VariableSites.meth", sep = "\t", header = TRUE)

rownames(annotations) <- annotations$ID

mat <- as.matrix(meth[,c(-1,-2,-3)])


datTraits <- datTraits[,c(-2)]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p", method = "spearman")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

ageSig.df <- data.frame(moduleTraitCor, moduleTraitPvalue)
ageSig.df$ages.pAdj <- p.adjust(ageSig.df$ages.1, method = "fdr")
ageSig.df$sig <- ifelse(ageSig.df$ages.pAdj < 1e-5, colnames(MEs), NA)


## plot the outer portion of Fig 2A

e <- ggplot(data = ageSig.df, aes(x = rownames(ageSig.df), y = -log10(ages.pAdj), size = abs(ages))) + 
    xlab("WGCNA Module") + 
    geom_point(color = stringr::str_remove_all(rownames(moduleTraitCor), "ME")) + ylab("-log10(p) for Spearman Correlation") +
    theme(axis.text.x = element_blank(), legend.text = element_text(size = 11),
              legend.title = element_text(size = 14, face = "bold"), axis.ticks.x = element_blank()) + 
    labs(size = "Abs. r") + xlab("")
e

ggsave(e, filename = "~/Parrott_Lab/AnoleAge/WGCNA/images/ModuleAgeSig.svg", device = "svg",
       width = 8, height = 4)


### Plot Blue eigenvalue versus age

annot_MEs <- cbind(annotations, MEs)

e.inlay <- ggplot(data = annot_MEs, aes(x = ages, y = MEblue)) + 
                   geom_point() + geom_smooth(method = "lm", formula = y ~ log(x), alpha = 0.2) +
                   xlab("Age (months)") + ylab("MEblue") + theme_bw() + 
                   theme(axis.text=element_text(size=6),
                   axis.title=element_text(size=8,face="bold"))

e.inlay

ggsave(e.inlay, filename = "~/Parrott_Lab/AnoleAge/WGCNA/images/blueMod_eigenVSage.svg", device = "svg",
       width = 2, height = 2)

logLM <- lm(data = annot_MEs, formula = MEblue ~ log(ages))
summary(logLM) ## p = 6e-14, R2 = 0.82


