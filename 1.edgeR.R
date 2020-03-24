library(data.table)
library(magrittr)
library(edgeR)
library(stringr)
library(clusterProfiler)
library(openxlsx)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library(org.Mm.eg.db)
library(SummarizedExperiment)

projectNumber <- "033"
orgDb <- "org.Mm.eg.db"

makeSumExp <- function(file, projectNumber)
{
  counts <- fread(file)
  countsMat <- counts[, -(1:6)] %>% as.matrix
  rownames(countsMat) <- str_remove(counts$Geneid, "\\.[[:digit:]]+$")
  geneInfo <- counts[, 1:6]
  
  coveredGenes <- rowSums(countsMat) != 0
  
  countsMat <- countsMat[coveredGenes, ]
  geneInfo <- geneInfo[coveredGenes, ]
  
  expSetup <- read.xlsx(paste0(projectNumber, "_Metadata.xlsx")) %>% data.table
  setkey(expSetup, Sequencing.ID)
  expSetup <- expSetup[, .(Sequencing.ID, Genotype, Treatment)]
  setnames(expSetup, c("Sequencing.ID"), c("id"))
  expSetup[, Genotype:=factor(Genotype, c("WT", "KO"))]
  expSetup[, Treatment:=factor(Treatment, levels = c("vechicle", "treatment"))]
  expSetup[, Group:=paste0(Genotype, "_", Treatment)]
  
  cDat <- as.data.frame(expSetup)
  rownames(cDat) <- cDat$id
  cDat <- cDat[colnames(countsMat), ]
  
  SummarizedExperiment(assays = list(counts = countsMat), colData = cDat, rowData = as(geneInfo, "GRanges"))
}


colScale <- scale_colour_manual(
  values = c("WT_vechicle" = "#a6cee3",
             "WT_treatment" = "#1f78b4",
             "KO_vechicle" = "#b2df8a",
             "KO_treatment" = "#33a02c"
             ),
  labels = c("WT_vechicle" = "WT, vehicle",
             "WT_treatment" = "WT, treatment",
             "KO_vechicle" = "KO, vehicle",
             "KO_treatment" = "KO, treatment"
  )
)

sumExp <- makeSumExp("featureCounts/counts.csv.gz", projectNumber)

# Sample 13 clusters on its own and is heavily enriched for genes related to inflammation and response to infection
# For this reason sample 13 is excluded
badSamples <- c("033_13")
sumExp <- sumExp[, !(colnames(sumExp) %in% badSamples)]
colnames(sumExp) <- sumExp$sampleId
metaData <- colData(sumExp) %>% as.data.table(keep.rownames = TRUE)
setnames(metaData, "rn", "Samples")
fwrite(metaData, file = "metaData.csv")

design <- model.matrix(~ 0 + Group, data = colData(sumExp))
colnames(design) <- str_replace_all(colnames(design), "Group|\\(|\\)", "")
colnames(design) <- str_replace(colnames(design), ":", "_")

sumExp$sampleId <- NULL
# Prepare the data: Filter and perform QC plotting
source("supportFunctions.R")
orgDb <- "org.Mm.eg.db"
y <- rnaPrepare(sumExp = sumExp,
             folder = "QC", 
             group = "Group",
             mdsLabel = "id", 
             mdsColour = "Group", 
             design = design,
             orgDb = orgDb,
             colScale = colScale)


# Fit the model
y <- estimateDisp(y, design = design) %T>%
  redirectToFile(file.path("QC", "BCV.pdf"), plotBCV)

efit <- glmQLFit(y, design = design) %T>%
  redirectToFile(file.path("QC", "QLDispersion.pdf"), plotQLDisp)

colnames(design)
ctrsts <- makeContrasts(
  TreatmentInWT = WT_treatment - WT_vechicle,
  TreatmentInKO = KO_treatment - KO_vechicle,
  KOInVehicle = KO_vechicle - WT_vechicle,
  KOInTreatment = KO_treatment - WT_treatment,
  Interaction = (KO_treatment - KO_vechicle) - (WT_treatment - WT_vechicle),
  levels = design
)

dgeResults <- apply(ctrsts, 2, . %>% 
                      glmQLFTest(glmfit = efit, contrast = .) %>% 
                      topTags(n = Inf, p.value = 1) %>% 
                      extract2("table") %>% 
                      as.data.table(keep.rownames = TRUE))

annotated <- lapply(dgeResults, annotateFromEns, orgDb = orgDb)

dir.create("edgeR_results", showWarnings = FALSE)
write.xlsx(annotated, file = paste0("edgeR_results/", projectNumber, "_RNA.xlsx"), asTable = TRUE)

save(y, efit, design, ctrsts, annotated, file = "edgerDataObjects.RData")
