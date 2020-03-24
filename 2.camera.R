library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(SummarizedExperiment)
library(reactome.db)
library(GO.db)
library(clusterProfiler)
library(openxlsx)
library(edgeR)

projectNumber <- "033"
load("edgerDataObjects.RData")
# Camera and fry enrichment analysis

cameraFn <- function(contrast, index, y, design)
{
  out <- camera(y, index = index, design = design, contrast = contrast) %>%
    data.table(keep.rownames = TRUE) %>% 
    setnames("rn", "ID") %>%
    extract(!is.na(PValue))
  out[order(PValue)]
}

#### Convert from ENSEMBL to ENTREZID
conv <- bitr(rownames(y$counts), fromType='ENSEMBL', toType='ENTREZID', OrgDb = "org.Mm.eg.db") %>%
  data.table(key = "ENSEMBL")

# Closer inspection reveales that few genes have multiple ENTREZ IDs
convNames <- conv[rownames(y), ENTREZID, mult = "first"]
y <- y[!duplicated(convNames), ]
rnams <- conv[rownames(y), ENTREZID, mult = "first"]
rnams[is.na(rnams)] <- "-1000"
rownames(y) <- rnams

#### Gene Ontology
goFn <- function(selectedTerms)
{
  keysGO <- keys(GO.db)
  termGO <- select(GO.db, keys=keysGO, columns=c("TERM", "ONTOLOGY")) %>% data.table
  termGO <- termGO[ONTOLOGY == selectedTerms]
  termGO[, ONTOLOGY:=NULL]
  setnames(termGO, "GOID", "ID")
  
  cyt.go.genes <- as.list(org.Mm.egGO2ALLEGS)
  cyt.go.genes <- cyt.go.genes[names(cyt.go.genes) %in% termGO$ID]
  cyt.go.genes <- Filter(. %>% length %>% is_greater_than(4), cyt.go.genes) # Remove small categories
  cyt.go.genes <- Filter(. %>% length %>% is_less_than(501), cyt.go.genes) # Remove large categories
  
  cameraGO <- apply(ctrsts, 2, cameraFn, index = cyt.go.genes, y = y, design = design)
  cameraGO <- lapply(cameraGO, merge, termGO, by = "ID", nomatch = FALSE)
  cameraGO <- lapply(cameraGO, extract, order(PValue, decreasing = FALSE)) 
  cameraGO
}
cameraGO_CC <- goFn("CC")
cameraGO_MF <- goFn("MF")
cameraGO_BP <- goFn("BP")

# Annotate with genes related to GO terms
annotateWithGenes <- function(tests, termList, fromType){
  conv <- bitr(unique(unlist(termList)), fromType = fromType, toType = 'SYMBOL', OrgDb = "org.Mm.eg.db") %>%
    data.table(key = fromType)
  pasteSymbols <- function(x) conv[termList[[x]], paste0(unique(SYMBOL), collapse = ", ")]
  termSymbols <- lapply(names(termList), pasteSymbols) %>% data.table
  termSymbols[, id:=names(termList)]
  setnames(termSymbols, c("symbolList", "id"))
  setkey(termSymbols, id)
  
  tests <- lapply(tests, copy)
  for (i in names(tests)){
    tests[[i]][, genesInTerm:=termSymbols[ID, symbolList]]
  }
  tests
}

cameraGO_CC_Annotated <- annotateWithGenes(cameraGO_CC, as.list(org.Mm.egGO2ALLEGS), fromType = 'ENTREZID')
cameraGO_MF_Annotated <- annotateWithGenes(cameraGO_MF, as.list(org.Mm.egGO2ALLEGS), fromType = 'ENTREZID')
cameraGO_BP_Annotated <- annotateWithGenes(cameraGO_BP, as.list(org.Mm.egGO2ALLEGS), fromType = 'ENTREZID')

save(cameraGO_CC_Annotated, cameraGO_BP_Annotated, cameraGO_MF, file = "cameraResults_annotated.RData")

write.xlsx(cameraGO_BP_Annotated, paste0("edgeR_results/", projectNumber, "_cameraGO_BP_with_genes.xlsx"), asTable = TRUE)
write.xlsx(cameraGO_CC_Annotated, paste0("edgeR_results/", projectNumber, "_cameraGO_CC_with_genes.xlsx"), asTable = TRUE)
write.xlsx(cameraGO_MF_Annotated, paste0("edgeR_results/", projectNumber, "_cameraGO_MF_with_genes.xlsx"), asTable = TRUE)
