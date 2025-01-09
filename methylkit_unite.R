library(methylKit)
library(dplyr)
library(data.table)

setDTthreads(threads=2)

##create list of files which will contribute to the datset
setwd("~/Parrott_Lab/Anole_Data/methylkit")

file.list <- list("10",
                  "15",
                  "19",
			        	  "23",
                  "28",
                  "33",
                  "41",
                  "46",
                  "49",
                  "51",
                  "53",
                  "61",
                  "62",
                  "64",
                  "66",
                  "67",
                  "70",
                  "76",
                  "84",
                  "85",
                  "95",
                  "96",
                  "107",
                  "110",
                  "111",
                  "116",
                  "117",
                  "120",
                  "122",
                  "125",
                  "126",
                  "127",
                  "128",
                  "129",
                  "131",
                  "132",
				          "133")

fileLOC <- paste("~/Parrott_Lab/Anole_Data/methylkit/1x/", file.list, "_CpG.txt", sep = "")

fileLOC.list <- as.list(fileLOC)

#create methylRawList object which will consist of methylation calls from each sample


methOBJ <- methRead(location = fileLOC.list,
                    sample.id = as.list(file.list),
                    assembly = "AnoSag2.1", 
                    context = "CpG",
                    mincov = 5,
                    treatment = rep(0, length(fileLOC.list)))  
					


# create unfiltered MethylBase file 
methFrame.total <- unite(methOBJ, destrand = FALSE, mc.cores = 2, min.per.group = 1L,
                         save.db = TRUE, suffix = "total_allSites_", dbdir = getwd())


## Create methylation matrix
mat.methFrame.total <- percMethylation(methFrame.total)


## Save data tables
methFrame.total.raw <- as(methFrame.total, "methylBase")

methFrame.coords <- methFrame.total.raw[,c(1,2,3,4)]

fwrite(cbind(methFrame.coords, mat.methFrame.total), 
       file = "~/Parrott_Lab/Anole_Data/methylkit/processed/labeled_allSites_stranded.tab", sep = "\t", quote = F)

## remove large tables from memory
rm("methFrame.total", "mat.methFrame.total", "methFrame.total.raw", "methFrame.coords")



#### Do same process for filtered dataset (5x coverage, 30 individuals)

methFrame.filt <- unite(methOBJ, destrand = TRUE, mc.cores = 2, min.per.group = 30L,
                         save.db = TRUE, suffix = "total_30L_", dbdir = getwd())

methFrame.filt <- filterByCoverage(methFrame.filt,lo.count = 5,lo.perc = NULL,
                            hi.count = NULL,hi.perc = 99.5)

methFrame.filt <- normalizeCoverage(methFrame.filt)

mat.methFrame.filt <- percMethylation(methFrame.filt)

methFrame.filt.raw <- as(methFrame.filt, "methylBase")

methFrame.filt.coords <- methFrame.filt.raw[,c(1,2,3,4)]

fwrite(cbind(methFrame.coords, mat.methFrame.total), 
       file = "~/Parrott_Lab/Anole_Data/methylkit/processed/labeled_total80perc.tab", sep = "\t", quote = F)
