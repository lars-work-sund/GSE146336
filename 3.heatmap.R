library(vegan)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(edgeR)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(clipr)
library(openxlsx)
library(PoiClaClu)
library(factoextra)

load("edgerDataObjects.RData")

cpmsBySample <- cpm(y, log = TRUE) %>% data.table(keep.rownames = TRUE)
setnames(cpmsBySample, "rn", "ENSEMBL")
expSetup <- fread("metaData.csv")

all(colnames(cpmsBySample)[-1] == expSetup$id)
expSetup <- expSetup[order(Genotype, Treatment, decreasing = TRUE)]

prepareForPlotting <- function(x)
{
  exprMat <- as.matrix(x[, -c(1)])
  rownames(exprMat) <- x$ENSEMBL
  visData <- exprMat
  visData <- t(scale(t(visData), center = T, scale = T))
  visData
}

annotDF <- as.data.frame(expSetup[, .(Treatment, Genotype)])
rownames(annotDF) <- expSetup$sampleId
annotDF$Treatment <- factor(annotDF$Treatment)
levels(annotDF$Treatment) <- c("treatment" = "Treatment", "vechicle" = "Vehicle")[levels(annotDF$Treatment)]
annotDF$Treatment <- relevel(annotDF$Treatment, "Vehicle")

annotDF$Genotype <- factor(annotDF$Genotype)
annotDF$Genotype <- relevel(annotDF$Genotype, "WT")

annotColours <- list(
  Treatment = c(
    "Vehicle" = "#e41a1c",
    "Treatment" = "#377eb8"
  ),
  Genotype = c(
    "WT" = "#4daf4a",
    "KO" = "#984ea3"
  )
)

myHeatmap <- function(x, showRow = FALSE)
{
  colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
  breaks <- c(min(x), seq(from = -2, to = 2, length.out = 253), max(x))
  
  pheatmap(x, 
           show_colnames = FALSE,
           show_rownames = showRow,
           breaks = breaks,
           clustering_method = "ward.D2",
           clustering_distance_rows="correlation",
           treeheight_row = 0,
           cluster_cols=FALSE, 
           col=colors,
           silent = TRUE
  )
}

plotData <- prepareForPlotting(cpmsBySample)
i <- "TreatmentInWT"
subData <- plotData[annotated[[i]][1:30, ENSEMBL], expSetup$Samples]
rownames(subData) <- annotated[[i]][1:30, fifelse(is.na(SYMBOL), "No Gene Symbol", SYMBOL)]
png(paste0("edgeR_results/heatmap_", i, "_limited_colour_scale_no_dendrogram.png"), width = 6.5, height = 5, units = "in", res = 1200)
myHeatmap(x = subData, showRow = TRUE) %>% print
dev.off()
