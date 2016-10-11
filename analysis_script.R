load("~/Documents/prad.Rda")
eset <- extract(prad, "rnaseq2genenorm")

mydata <- exprs(eset)
clinical <- pData(eset)

## check that sample_code is primary tumor (after loading RTCGAToolbox)
samptab

## only looking at patients with tumors in the sample code
tum_clin <- clinical[clinical$sample_code == "01", ]
dim(tum_clin)

## take participant IDs from those who have tumors
participantIDS <- rownames(tum_clin)
mydata[ , participantIDS]

## Subset original data for those with tumors
mysubset <- mydata[ , participantIDS]

## load and clean race variable
racevar <- read.csv("racevariable.csv", stringsAsFactors = FALSE)
racevar[, 1] <- paste0(racevar[, 1], "-01")
racevar[, 1] <- gsub(".", "-", racevar[, 1], fixed=TRUE)
names(racevar) <- c("patientID", "race")

## inspect racevar and mysubset to understand how to merge
head(racevar)
nrow(racevar)

head(tum_clin)
nrow(tum_clin)

tum_clin$patientID <- rownames(tum_clin)
newClin <- merge(x = tum_clin, y = racevar, by = "patientID")
# newClin2 <- merge(x = tum_clin, y = racevar, by.x = "patientID", by.y = "patientID")

## Replace clinical data and expression data in the original eset object
pData(eset) <- newClin
exprs(eset) <- mysubset

## Create a fake variable
newClin$fakeTreatment <- rbinom(n = nrow(newClin), size = 1, 0.5)

## Create a new race variable that is a factor
newClin$newRace <- factor(newClin$race.y)

## Manually change the factor levels of the race variable
levels(newClin$newRace) <- c("other", "other", "black", "white")

## Relevel the factor variable to make "white" the reference group
newClin$newRace <- relevel(newClin$newRace, ref = "white")

## Check reference group is "white"
contrasts(newClin$newRace)

## look at all the variable names in the dataset/dataframe
names(newClin)

## search through the variables in the dataset 
grep("score", names(newClin), value = TRUE, ignore.case = TRUE)
grep("clinic", names(newClin), value = TRUE, ignore.case = TRUE)

table(newClin$gleason_score_combined)
class(newClin$gleason_score_combined)

newClin$gleason_score_combined <- as.numeric(newClin$gleason_score_combined)

table(newClin$gleason_score_combined)
factor(newClin$gleason_score_combined)

## Formula notation
myFormula <- fakeTreatment ~ newRace + factor(patient.clinical_cqcf.histological_type)
class(myFormula)

model1 <- glm(fakeTreatment ~ newRace + factor(gleason_score_combined), data = newClin, family = binomial(link = "logit"))
model2 <- glm(myFormula, data = newClin, family = binomial(link = "logit"))

summary(model1)
summary(model2)

exp(coef(model1))
exp(cbind(OR = coef(model1), confint(model1)))

## Cross tab of country and race
xtabs(~ patient.clinical_cqcf.country + newRace, data = newClin)
prop.table(table(newClin$newRace, newClin$patient.clinical_cqcf.country), 1)

## proportions of each category in the data
table(newClin$newRace)
prop.table(table(newClin$newRace))

## get all IDs for those who are black
blackIDS <- newClin$patientID[newClin$newRace == "black"]

blackExpressionLevels <- mysubset[, colnames(mysubset) %in% blackIDS]
mean(blackExpressionLevels[rownames(blackExpressionLevels) == "ACPP",])

whiteIDS <- newClin$patientID[newClin$newRace == "white"]

whiteExpressionLevels <- mysubset[, colnames(mysubset) %in% whiteIDS]
mean(whiteExpressionLevels[rownames(whiteExpressionLevels) == "ACPP",])

geneList <- sample(rownames(mysubset), 10)

## Get means of all the genes in the data by race
# result <- sapply(rownames(mysubset), function(gene) {
#     c(mean(whiteExpressionLevels[rownames(whiteExpressionLevels) == gene]),
#       mean(blackExpressionLevels[rownames(blackExpressionLevels) == gene,]))
# })
# rownames(result) <- c("white", "black")
# result <- t(result)
# write.csv(result, "meanExpressionWB.csv")
# meanExp <- read.csv("meanExpressionWB.csv")

head(meanExp[order(meanExp$perChange, decreasing = TRUE),], 40)
meanExp2 <- meanExp[!is.infinite(meanExp$perChange), ]

psadat <- data.frame(psa = as.numeric(as.character(newClin$patient.stage_event.psa.psa_value)),
                     race = newClin$newRace,
                     PVT1=mysubset["PVT1", ])

# Add race variable to psa value dataset (alternative)
# psadat$race <- newClin$newRace[match(rownames(psadat), newClin$patientID)]

group_by(psadat, race) %>% summarise(avg_psa = mean(psa, na.rm = TRUE),
                                     n = n(),
                                     stddev = sd(psa, na.rm = TRUE))

# Explore "stage" variables
grep("stage", names(newClin), value =TRUE)

psadat.complete <- psadat[complete.cases(psadat), ]
plot(PVT1 ~ psa, data=psadat.complete, xlab="clinical PSA", ylab="PVT1 tumor expression", log="xy")
fit <- lowess(x=psadat.complete$psa, y=psadat.complete$PVT1)
lines(fit, col="red", lw=3)

cor.test(x=psadat$PVT1, y=psadat$psa, method="spearman")

boxplot(psadat$PVT1 ~ psadat$race, ylab = "PSA")

racevar = factor(mergeVecs(as.character(newClin$patient.clinical_cqcf.race), as.character(newClin$patient.race)))

library(DESeq2)
countData <- exprs(eset)
colData <- pData(eset)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ gleason_score + batch_number + racevar)
dds <- DESeq(dds)

# Exon Analysis -----------------------------------------------------------

library(TCGAbiolinks)
prad <- TCGAquery(tumor = "prad", platform = "IlluminaHiSeq_RNASeqV2",
                  level = 3)
TCGAdownload(prad, path = path.expand(dataFolder), type = "bt.exon_quantification")

## Install packages
# BiocInstaller::biocLite("GenomicFeatures")
# BiocInstaller::biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("GenomicFeatures")
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # database package 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene #shorthand (for convenience)
txdb

seqlevels(txdb) <- "chr8"
# seqlevels(txdb) <- seqlevels0(txdb)
exon <- exonsBy(txdb, by="gene")

## Entrez GeneID for PVT1 = 5820
allExons <- exon$`5820`

## RefSeq Genes
exon9 <- allExons[10]
exon8 <- allExons[9]

select(txdb, keys = "5820", columns = c("EXONID", "EXONRANK"), keytype = "GENEID")

dataFolder <- "~/Documents/data/"
exonFiles <- dir(dataFolder, pattern = ".txt", full.names = TRUE)

exonFiles

# BiocInstaller::biocLite("waldronlab/TCGAmisc")
library(TCGAmisc)

myGRL <- readExonFiles(exonFiles)

Exon9s <- lapply(myGRL, function(sample) {
    # sample <- sample[seqnames(sample) == "chr8"]
    sample[overlapsAny(sample, exon9, type = "any")]
    # sample[sample %within% exon9]
})
Exon8s <- lapply(myGRL, function(sample) {
    # sample <- sample[seqnames(sample) == "chr8"]
    sample[overlapsAny(sample, exon8, type = "any")]
    # sample[sample %within% exon9]
})

byExons <- lapply(myGRL, function(sample){
    sample[overlapsAny(sample, allExons, type = "any")]
})

sapply(1:9, function(x, i) {x$RPKM[i]}, x = byExons)
exonsByRPK <- do.call(cbind, lapply(byExons, function(x) x$RPKM))
rownames(exonsByRPK) <- paste0("exon", 1:9)
apply(exonsByRPK, 1, mean)

rpkms9 <- sapply(Exon9s, function(x) {x$RPKM})
allRPK9 <- as.numeric(rpkms9)
e9rpk <- mean(allRPK9, na.rm = TRUE)

rpkms8 <- sapply(Exon8s, function(x) {x$RPKM})
allRPK8 <- as.numeric(rpkms8)
e8rpk <- mean(allRPK8, na.rm = TRUE)

t.test(allRPK9, allRPK8)



library(readr)
racevar <- read_csv("https://raw.githubusercontent.com/lwaldron/tcga_prad/master/racevariable.csv")
racevar$patientID <- barcode(racevar$X)

exonpatient <- barcode(names(Exon9s))

races <- racevar$race[match(exonpatient, racevar$patientID)]

idtorace <- data.frame(cbind(sampID =names(Exon9s), races), stringsAsFactors = FALSE)
