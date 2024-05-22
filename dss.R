library(bsseq)
library(DSS)
library(stringr)
library(data.table)

design = data.frame("Sex" = c("F","F","F","M","M","F","M","F","M","F","M","F","M","F","F","F","M","M",
                              "F","M","F","M","M","F","F","F","F","M","M","M","M","M","M","F","M","F","F"),
                    "Age" = c(36,36,36,36,36,36,48,48,48,48,48,48,18,18,18,18,
                              18,18,1,1,1,1,7,7,7,7,7,7,7,7,7,60,60,60,60,60,60))

X = model.matrix(~Sex + Age + Sex*Age, design)

path = file.path("/scratch/es88065/store_temp/methExtract/dss")

fileList <- list.files(path = path, pattern = ".dss")
fileList <- str_sort(fileList, numeric = TRUE)

sampleIDs <- c("S10", "S15", "S19", "S23", "S28", "S33", "S41", "S46", "S49", "S51", "S53", "S61",
               "S62", "S64", "S66", "S67", "S70", "S76", "S84", "S85", "S95", "S96", "S107", "S110",
               "S111", "S116", "S117", "S120", "S122", "S125", "S126", "S127", "S128", "S129", "S131", "S132", "S133")


for (i in 1:length(sampleIDs)) {
    assign(paste0(sampleIDs[i]), fread(file.path(path, fileList[i]), 
                                       header=FALSE, sep = "\t", 
                                       col.names = c("chr", "pos", "N", "X")))
}

 
BSobj = makeBSseqData(list(S10, S15, S19, S23, S28, S33, S41, S46, S49, S51, S53, S61,
                            S62, S64, S66, S67, S70, S76, S84, S85, S95, S96, S107, S110,
                            S111, S116, S117, S120, S122, S125, S126, S127, S128, S129, S131, S132, S133),
                       sampleNames = c("S10", "S15", "S19", "S23", "S28", "S33", "S41", "S46", "S49", "S51", "S53", "S61",
                                       "S62", "S64", "S66", "S67", "S70", "S76", "S84", "S85", "S95", "S96", "S107", "S110",
                                       "S111", "S116", "S117", "S120", "S122", "S125", "S126", "S127", "S128", "S129", "S131", "S132", "S133"))

rm(list=ls(S10, S15, S19, S23, S28, S33, S41, S46, S49, S51, S53, S61,
           S62, S64, S66, S67, S70, S76, S84, S85, S95, S96, S107, S110,
           S111, S116, S117, S120, S122, S125, S126, S127, S128, S129, S131, S132, S133))

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~Sex+Age+Sex:Age)

save(DMLfit, "/scratch/es88065/AnoleAge/DMLfit_logAge_smoothed200.fit")


load(file = "E:/anole/CpG/cov/DMLfit_logAge_smoothed200.fit")

saveToBed <- function(DMLframe, file){
    DMLframe$start <- DMLframe$pos - 1
    columns <- DMLframe[,c("chr", "start", "pos", "stat", "fdrs")]
    write.table(columns, file = file, quote = F, sep = "\t", row.names = F, col.names = F)
}

options(scipen = 999)

DMLtest.age <- DMLtest.multiFactor(DMLfit, coef=3)
sig.age.sites <- DMLtest.age[which(DMLtest.age$fdrs < 0.05),]
saveToBed(sig.age.sites, file = "C:/Users/sheal/Desktop/AnoleAge/misc/DSSdata/LogAgeSites_smoothed200.bed")
rm(DMLtest.age)

DMLtest.sex <- DMLtest.multiFactor(DMLfit, coef=2)
sig.sex.sites <- DMLtest.sex[which(DMLtest.sex$fdrs < 0.05),]
saveToBed(sig.sex.sites, file = "C:/Users/sheal/Desktop/AnoleAge/misc/DSSdata/LogSexSites_smoothed200.bed")
rm(DMLtest.sex)

DMLtest.sexAge <- DMLtest.multiFactor(DMLfit, coef=4)
sig.interact.sites <- DMLtest.sexAge[which(DMLtest.sexAge$fdrs < 0.05),]
saveToBed(sig.interact.sites, file = "C:/Users/sheal/Desktop/AnoleAge/misc/DSSdata/LogInteractSites_smoothed200.bed")
rm(DMLtest.sexAge)

### Call DMRs from fit DML models

saveDMRToBed <- function(DMRframe, file){
    columns <- DMRframe[,c("chr", "start", "end", "areaStat", "nCG")]
    write.table(columns, file = file, quote = F, sep = "\t", row.names = F, col.names = F)
}


age.DMR <- callDMR(sig.age.sites, p.threshold=0.05)
saveDMRToBed(age.DMR, file = "C:/Users/sheal/Desktop/AnoleAge/Fig1-DSSage/data/LogAge_DMRs_smoothed200.bed")


sex.DMR <- callDMR(sig.sex.sites, p.threshold=0.05)
saveDMRToBed(sex.DMR, file = "C:/Users/sheal/Desktop/AnoleAge/Fig2/data/LogSex_DMRs_smoothed200.bed")


sexAge.DMR <- callDMR(sig.interact.sites, p.threshold=0.05)
saveDMRToBed(sex.DMR, file = "C:/Users/sheal/Desktop/AnoleAge/Fig4-SexDiff/data/LogSex_DMRs_smoothed200.bed")




