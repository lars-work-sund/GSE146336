redirectToFile <- function(y, file, fun, last = TRUE, ...)
{
  pdf(file, width = 4, height = 4) 
  fun(y, ...) 
  if(last) dev.off()
}

plotDensity <- function(y, name)
{
  pD <- reshape2::melt(cpm(y, normalized.lib.sizes = TRUE, log = TRUE))
  p <- ggplot(pD, aes(value)) + 
    geom_density() + 
    facet_wrap(~Var2)
  
  print(p)
}

plotMDfun <- function(y, name)
{
  oldpar <- par()$mfrow
  par(mfrow = c(2, 2))
  for (i in seq_len(ncol(y))){
    plotMD(y, column = i)
    abline(h = 0)
  }
  par(mfrow = oldpar)
}

writeCPM <- function(y, file, orgDb = orgDb)
{
  cpmData <- cpm(y, log = TRUE) %>% data.table(keep.rownames = TRUE, key = "rn")
  setnames(cpmData, "rn", "ENSEMBL")
  conv <- bitr(geneID = rownames(y), 
               fromType = "ENSEMBL", 
               toType = c("SYMBOL", "GENENAME"), 
               OrgDb = orgDb) %>%
    data.table(key = "ENSEMBL")
  cpmData <- conv[cpmData]
  write.xlsx(cpmData, file = file, asTable = TRUE)
  invisible(TRUE)
}

mdsFormat <- function(y, columns)
{
  mdsData <- plotMDS(y, ndim = 3, plot = FALSE)
  mdsData <- mdsData$cmdscale.out %>% data.table(keep.rownames = TRUE)
  mdsData <- mdsData[as.data.table(y$samples[, columns], keep.rownames = TRUE), on = "rn"]
  setnames(mdsData, c("rn", "V1", "V2", "V3"), c("ID", "dim1", "dim2", "dim3"))
  mdsData
}

mdsPlot <- function(mdsData, colScale, label, colour)
{
  pBase <- ggplot(mdsData, aes_string(label = label, colour = colour)) + 
    geom_point(size = 5) +
    geom_label(show.legend = FALSE, size = 2, position = position_dodge(0.5)) +
    theme_bw() +
    colScale
  
  plotter <- function(d1, d2)
  {
    p <- pBase %+% 
      aes_string(x = d1, y = d2) +
      xlab(paste("Leading logFC dim", str_remove(d1, "dim"))) +
      ylab(paste("Leading logFC dim", str_remove(d2, "dim")))
    print(p)
  }
  plotter("dim1", "dim2")
  plotter("dim1", "dim3")
  plotter("dim2", "dim3")
}

poisHeatmap <- function(y, outFile, sampleData)
{
  poisd <- PoissonDistance(t(y$counts))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- colnames(cpm(y))
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  annot <- as.data.frame(sampleData)
  annot$id <- NULL
  
  pals <- c("Accent", "Dark2", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
  pals <- sample(pals, size = ncol(annot), replace = ifelse(ncol(annot) > length(pals), TRUE, FALSE))
  makeColours <- function(x, pal){
    nams <- unique(x)
    n <- length(nams)
    if (n > 8){
      out <- colors[floor(seq(1, length(colors), length.out = n))]
    } else {
      out <- brewer.pal(max(n, 3), pal)
      if (n < 3) {
        out <- out[seq_len(n)]
      }
    }
    names(out) <- nams
    out
  }
  annotCols <- Map(makeColours, annot, pals)
  
  pheatmap(mat = samplePoisDistMatrix, 
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           color = colors, 
           filename = outFile, 
           annotation_row = annot, 
           annotation_colors = annotCols)
}

combinedMDS <- function(y, label, colour, colScale)
{
  y %>%
    mdsFormat(c(colour, label)) %>%
    mdsPlot(colScale = colScale, 
            label = label, 
            colour = colour)
}

annotateFromEns <- function(rnaRes, orgDb = orgDb)
{
  dat <- copy(rnaRes)
  if (nrow(dat) == 0) return(dat)
  setnames(dat, names(dat)[1], "ENSEMBL")
  ens2symbol <- bitr(dat$ENSEMBL, fromType='ENSEMBL', toType='SYMBOL', OrgDb = orgDb)
  ens2symbol %<>% data.table %>% setkey(ENSEMBL)
  ens2description <- bitr(dat$ENSEMBL, fromType='ENSEMBL', toType='GENENAME', OrgDb = orgDb)
  ens2description %<>% data.table %>% setkey(ENSEMBL)
  dat[, SYMBOL:=ens2symbol[ENSEMBL, SYMBOL, mult = 'first']]
  dat[, GENENAME:=ens2description[ENSEMBL, GENENAME, mult = 'first']]
  dat
}

rnaPrepare <- function(sumExp, folder = "qcFigures", group, mdsLabel, mdsColour, colScale, design, orgDb) # Add support for min.count and design in filterByExpr
{
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  dgeList <- DGEList(counts = assay(sumExp, "counts"),
                     group = colData(sumExp)[, group], 
                     genes = as.data.frame(rowRanges(sumExp)),
                     samples = colData(sumExp))
  dgeList %T>%
    redirectToFile(file.path(folder, "beforeFiltering_MD.pdf"), plotMDfun) %T>%
    redirectToFile(file.path(folder, "beforeFiltering_density.pdf"), plotDensity) %>%
    extract(filterByExpr(., design = design), , keep.lib.sizes=FALSE) %>%
    calcNormFactors %T>%
    redirectToFile(file.path(folder, "afterFiltering_MD.pdf"), plotMDfun) %T>%
    redirectToFile(file.path(folder, "afterFiltering_density.pdf"), plotDensity) %T>%
    writeCPM(file.path(folder, "cpmData.xlsx"), orgDb = orgDb) %T>%
    redirectToFile(file.path(folder, "mdsPlot.pdf"), combinedMDS, label = mdsLabel, colour = mdsColour, colScale = colScale) %T>%
    poisHeatmap(outFile = file.path(folder, "sampleHeatmap.pdf"), sampleData = colData(sumExp))
}