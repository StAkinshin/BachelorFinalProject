library(tximport)
library(rhdf5)
library(DESeq2)
library(clusterProfiler)
library(gprofiler2)
library(ggplot2)
library(VennDiagram)
library(mixOmics)

#results - res
#significant results - sign.res

#METADATA
metadata.gamma <-  read.csv('C:/Users/ПК/Documents/bachelor/data/reference_data/metadata_gamma.csv', 
                            row.names = 1,
                            sep = ';')
metadata.proton <- read.csv('C:/Users/ПК/Documents/bachelor/data/reference_data/metadata_proton.csv', 
                            row.names = 1,
                            sep = ';')

#FUNCTIONS
STARtoDDS <- function(STARCountsFolderPath, MetaData) {
  
  CountsFiles <- list.files(STARCountsFolderPath,
                            pattern = '*.txt',
                            full.names = T
                            )
  CountsList <- lapply(CountsFiles, function(file) read.table(file,
                                   header = T,
                                   row.names = 1,
                                   na.strings = 'NA',
                                   dec = '.'))
  Matrix <- do.call(cbind, CountsList)
  colnames(Matrix) <- rownames(MetaData)
  DDS <-  DESeqDataSetFromMatrix(Matrix,
                                 colData = MetaData,
                                 design = ~ condition)
  DDS <- DESeq(DDS)
  return(DDS)
}

KALLISTOtoDDS <- function(KallistoFiles, MetaData) {
    names(KallistoFiles) <- rownames(MetaData)
    TXI <- tximport(KallistoFiles, type = "kallisto", txOut = TRUE)
    DDS <- DESeqDataSetFromTximport(TXI, colData = MetaData, design = ~condition)
    DDS$condition <- relevel(DDS$condition, ref = "control")
    DDS <- DESeq(DDS)
    return(DDS)
}

###MAIN###
##UPLOAD & DDS ESTIMATION##
dds.star.gamma <- STARtoDDS('C://Users/ПК/Documents/bachelor/data/2_Aligning/star+featureCounts/counts_gamma/', 
                            metadata.gamma)

dds.star.proton <- STARtoDDS('C://Users/ПК/Documents/bachelor/data/2_Aligning/star+featureCounts/counts_proton',
                             metadata.proton)

##RESULTS COMPILATION##
res.star.gamma <- results(dds.star.gamma)
res.star.proton <- results(dds.star.proton)
res.kallisto.gamma <- results(dds.kallisto.gamma)
res.kallisto.proton <- results(dds.kallisto.proton)

summary(res.star.gamma, alpha = 0.05)
summary(res.star.proton, alpha = 0.05)
summary(res.kallisto.gamma, alpha = 0.05)
summary(res.kallisto.proton, alpha = 0.05)

sign.res.star.gamma <- subset(res.star.gamma, padj < 0.05 & abs(log2FoldChange) > 1)
sign.res.star.proton <- subset(res.star.proton, padj < 0.05 & abs(log2FoldChange) > 1)
sign.res.kallisto.gamma <- subset(res.kallisto.gamma, padj < 0.05 & abs(log2FoldChange) > 1)
sign.res.kallisto.proton <- subset(res.kallisto.proton, padj < 0.05 & abs(log2FoldChange) > 1)

upreg.sign.res.star.gamma <- subset(sign.res.star.gamma, log2FoldChange > 0)
downreg.sign.res.star.gamma <- subset(sign.res.star.gamma, log2FoldChange < 0)
upreg.sign.res.star.proton <- subset(sign.res.star.proton, log2FoldChange > 0)
downreg.sign.res.star.proton <- subset(sign.res.star.proton, log2FoldChange < 0)

#НАДО ДОСТАТЬ ИНФУ ПРО ОБА ВОЗДЕЙСТВИЯ
common.upreg.sign.star <- intersect(rownames(upreg.sign.res.star.gamma), rownames(upreg.sign.res.star.proton))
common.downreg.sign.star <- intersect(rownames(downreg.sign.res.star.gamma), rownames(downreg.sign.res.star.proton))

unique.upreg.sign.star.gamma <- sign.res.star.gamma[setdiff(rownames(upreg.sign.res.star.gamma), rownames(upreg.sign.res.star.proton)), ]
unique.downreg.sign.star.gamma <- sign.res.star.gamma[setdiff(rownames(downreg.sign.res.star.gamma), rownames(downreg.sign.res.star.proton)), ]
unique.upreg.sign.star.proton <- sign.res.star.proton[setdiff(rownames(upreg.sign.res.star.proton), rownames(upreg.sign.res.star.gamma)), ]
unique.downreg.sign.star.proton <- sign.res.star.proton[setdiff(rownames(downreg.sign.res.star.proton), rownames(downreg.sign.res.star.gamma)), ]

##RESULTS VISUALIZATION##
#цвета презентации '#006668', '#003d3e'
#для протонов '#e59872'
#для гамма '#7fbf5a'

plotMA(
       res.star.gamma,
       alpha = 0.05,
       colNonSig = '#E2E3BE',
       colSig = '#7fbf5a',
       colLine = '#003d3e',
       ylim = c(-5, 9)
       )

plotMA(
  res.star.proton,
  alpha = 0.05,
  colNonSig = '#E2E3BE',
  colSig = '#e59872',
  colLine = '#003d3e',
  ylim = c(-5, 9)
)

vst.star.gamma <- vst(dds.star.gamma)
vst.star.proton <- vst(dds.star.proton)
vst.kallisto.gamma <- vst(dds.kallisto.gamma)
vst.kallisto.proton <- vst(dds.kallisto.proton)

plotPCA(vst.star.gamma)
plotPCA(vst.star.proton)

VennDiagram::venn.diagram(
  x = list(
    "Протоны" = rownames(upreg.sign.res.star.proton),
    "Гамма" = rownames(upreg.sign.res.star.gamma)
  ),
  category.names = c("Протоны", "Гамма"),
  filename = 'C://Users/ПК/Documents/bachelor/my_papers/Картиночки/STAR_UpRegVenn2.png',
  disable.logging = T,
  cat.cex = 0,
  lwd = 0,
  lty = 'blank',
  fill = c('#e59872','#7fbf5a'),
  cex = 3,
  fontface = "plain",
  rotation.degree = 180
)

VennDiagram::venn.diagram(
  x = list(
    "Гамма" = rownames(downreg.sign.res.star.gamma),
    "Протоны" = rownames(downreg.sign.res.star.proton)
  ),
  category.names = c("Гамма", "Протоны"),
  filename = 'C://Users/ПК/Documents/bachelor/my_papers/Картиночки/STAR_DownRegVenn2.png',
  disable.logging = T,
  cat.cex = 0,
  lwd = 2,
  lty = 'blank',
  fill = c('#7fbf5a', '#e59872'),
  cex = 3,
  fontface = "plain",
  rotation.degree = 180
)

## GENE CLUSTERING
rld.star.gamma = rlog(dds.star.gamma, blind = FALSE)
de.star.gamma = (res.star.gamma$padj <= 0.05)
de.star.gamma[is.na(de.star.gamma)] = FALSE

mixOmics::cim(t(assay(rld.star.gamma)[de.star.gamma, ]))

rld.star.proton = rlog(dds.star.proton, blind = FALSE)
de.star.proton = (res.star.proton$padj <= 0.05)
de.star.proton[is.na(de.star.proton)] = FALSE

mixOmics::cim(t(assay(rld.star.proton)[de.star.proton, ]))

##GO##
library(AnnotationHub)
library(biomaRt)
library(GOSemSim)
library(clusterProfiler)

genelist.star.gamma <- sign.res.star.gamma[,2]
names(genelist.star.gamma) = rownames(sign.res.star.gamma)
genelist.star.gamma = sort(genelist.star.gamma, decreasing = T)

genelist.star.proton <- sign.res.star.proton[,2]
names(genelist.star.proton) = rownames(sign.res.star.proton)
genelist.star.proton = sort(genelist.star.proton, decreasing = T)

ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")
searchDatasets(ensembl_plants, pattern = "Hordeum")
ensembl_hvulgare <- useEnsemblGenomes(biomart = "plants_mart", 
                                      dataset = "hvulgare_eg_gene")

hub <- AnnotationHub()
q <- query(hub, "Hordeum")
id <- q$ah_id[length(q)]
Hvulgare <- hub[[id]]

GSEA(
  genelist.star.gamma,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE,
)