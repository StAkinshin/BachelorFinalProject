library(DESeq2)

library(ggplot2)
library(gplots)
library(patchwork)

library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)

library(tximport)
library(rhdf5)
library(clusterProfiler)
library(gprofiler2)
library(ape)

#FUNCTIONS
tx2gene <- read.csv2("C:/Users/ПК/Documents/bachelor/data/reference_data/MappingIDs.csv") %>%
  dplyr::select(BaRT_ID, HORVU_ID) %>%
  dplyr::rename(transcript_id = BaRT_ID, gene_id = HORVU_ID)


KALLISTOtoDDS <- function(KallistoFilesPath, MetaData) {
  SampleNames <- rownames(MetaData)
  KallistoFiles <- file.path(KallistoFilesPath, SampleNames, "abundance.h5")
  names(KallistoFiles) <- rownames(MetaData)
  
  TXI <- tximport(KallistoFiles,
                  type = "kallisto",
                  tx2gene = tx2gene,
                  ignoreTxVersion = TRUE)
  
  DDS <- DESeqDataSetFromTximport(TXI, colData = MetaData, design = ~condition)
  DDS$condition <- relevel(DDS$condition, ref = "control")
  DDS <- DESeq(DDS)
  return(DDS)
}

FromResToGeneList <- function(Res){
  GeneList <- Res[,2]
  names(GeneList) = rownames(Res)
  GeneList = sort(GeneList, decreasing = T)
}

GOEA <- function(goSource, GeneList, BackgroundGenes) {
  
  gene2go <- goSource %>%
    filter(ensembl_gene_id %in% GeneList) %>%
    dplyr::select(ensembl_gene_id, go_id) %>%
    group_by(ensembl_gene_id) %>%
    summarise(go_ids = list(go_id)) %>%
    tibble::deframe()
  
  enricher(
    gene = names(gene2go),
    universe = if (!is.null(BackgroundGenes)) BackgroundGenes else unique(goSource$ensembl_gene_id),
    TERM2GENE = goSource[, c("go_id", "ensembl_gene_id")], 
    TERM2NAME = goSource[, c("go_id", "name_1006")],
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
}

clean_dotplot <- function(goea_obj, title = "", n = 10) {
  dotplot(goea_obj, showCategory = n) +
    coord_flip() +
    labs(title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45,    # Diagonal angle
        hjust = 1,     # Adjusts position so labels are readable
        vjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#METADATA
metadata.gamma <-  read.csv('C:/Users/ПК/Documents/bachelor/data/reference_data/metadata_gamma.csv', 
                            row.names = 1,
                            sep = ';')
metadata.proton <- read.csv('C:/Users/ПК/Documents/bachelor/data/reference_data/metadata_proton.csv', 
                            row.names = 1,
                            sep = ';')
metadata.ref <- read.csv('C:/Users/ПК/Documents/bachelor/data/reference_data/metadata_ref.csv', 
                         row.names = 1,
                         sep = ';')

dds.kallisto.ref <- KALLISTOtoDDS('C://Users/ПК/Documents/bachelor/data/2_Aligning/kallisto/', 
                                  metadata.ref)

background.goea.kallisto <- rownames(DDS)[rowMeans(counts(DDS, normalized=TRUE)) > 1]
background.goea.star <- rownames(read.csv('C:/Users/ПК/Documents/bachelor/data/reference_data/backgroundgenes_STAR.csv', 
                                          row.names = 1))


##GENE ONTOLOGY DATABASES
plants.ensembl <- useEnsemblGenomes(biomart = "plants_mart")
searchDatasets(plants.ensembl, pattern = "Hordeum")
hv.ensembl <- useEnsemblGenomes(biomart = "plants_mart", 
                                dataset = "hvulgare_eg_gene")

ant.go.hv <- getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name",
                                  "go_id",
                                  "name_1006",
                                  "namespace_1003"),
                   mart = hv.ensembl)
ant.go.hv <- ant.go.hv %>% filter(go_id != "") #ant - annotation
bp.ant.go.hv <- ant.go.hv %>% filter(namespace_1003 == "biological_process")
cc.ant.go.hv <- ant.go.hv %>% filter(namespace_1003 == "cellular_component")
mf.ant.go.hv <- ant.go.hv %>% filter(namespace_1003 == "molecular_function")


###MAIN###
##UPLOAD & DDS ESTIMATION##
dds.kallisto.gamma <- KALLISTOtoDDS('C://Users/ПК/Documents/bachelor/data/2_Aligning/kallisto/', 
                            metadata.gamma)
dds.kallisto.proton <- KALLISTOtoDDS('C://Users/ПК/Documents/bachelor/data/2_Aligning/kallisto/', 
                                    metadata.proton)

##RESULTS COMPILATION##
res.kallisto.gamma <- results(dds.kallisto.gamma)
res.kallisto.proton <- results(dds.kallisto.proton)

sign.res.kallisto.gamma <- subset(res.kallisto.gamma, padj < 0.05 & abs(log2FoldChange) > 1)
sign.res.kallisto.proton <- subset(res.kallisto.proton, padj < 0.05 & abs(log2FoldChange) > 1)

upreg.sign.res.kallisto.gamma <- subset(sign.res.kallisto.gamma, log2FoldChange > 0)
downreg.sign.res.kallisto.gamma <- subset(sign.res.kallisto.gamma, log2FoldChange < 0)
upreg.sign.res.kallisto.proton <- subset(sign.res.kallisto.proton, log2FoldChange > 0)
downreg.sign.res.kallisto.proton <- subset(sign.res.kallisto.proton, log2FoldChange < 0)

genelist.upreg.kallisto.gamma <- rownames(upreg.sign.res.kallisto.gamma)
genelist.downreg.kallisto.gamma <- rownames(downreg.sign.res.kallisto.gamma)
genelist.upreg.kallisto.proton <- rownames(upreg.sign.res.kallisto.proton)
genelist.downreg.kallisto.proton <- rownames(downreg.sign.res.kallisto.proton)

#write.table(x, file, append = FALSE, sep = " ", dec = ".",
#            row.names = TRUE, col.names = TRUE)

geneFC.upreg.kallisto.gamma <- FromResToGeneList(upreg.sign.res.kallisto.gamma)
geneFC.downreg.kallisto.gamma <- FromResToGeneList(downreg.sign.res.kallisto.gamma)
geneFC.upreg.kallisto.proton <- FromResToGeneList(upreg.sign.res.kallisto.proton)
geneFC.downreg.kallisto.proton <- FromResToGeneList(downreg.sign.res.kallisto.proton)


###
bp.go.upreg.kallisto.gamma <- bp.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.upreg.kallisto.gamma)
cc.go.upreg.kallisto.gamma <- cc.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.upreg.kallisto.gamma)
mf.go.upreg.kallisto.gamma <- mf.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.upreg.kallisto.gamma)

bp.go.downreg.kallisto.gamma <- bp.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.downreg.kallisto.gamma)
cc.go.downreg.kallisto.gamma <- cc.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.downreg.kallisto.gamma)
mf.go.downreg.kallisto.gamma <- mf.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.downreg.kallisto.gamma)

bp.go.upreg.kallisto.proton <- bp.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.upreg.kallisto.proton)
cc.go.upreg.kallisto.proton <- cc.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.upreg.kallisto.proton)
mf.go.upreg.kallisto.proton <- mf.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.upreg.kallisto.proton)

bp.go.downreg.kallisto.proton <- bp.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.downreg.kallisto.proton)
cc.go.downreg.kallisto.proton <- cc.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.downreg.kallisto.proton)
mf.go.downreg.kallisto.proton <- mf.ant.go.hv %>% filter(ensembl_gene_id %in% genelist.downreg.kallisto.proton)

####
bp.goea.upreg.kallisto.gamma <- GOEA(bp.ant.go.hv, genelist.upreg.kallisto.gamma, background.goea.star)
cc.goea.upreg.kallisto.gamma <- GOEA(cc.ant.go.hv, genelist.upreg.kallisto.gamma, background.goea.star)
mf.goea.upreg.kallisto.gamma <- GOEA(mf.ant.go.hv, genelist.upreg.kallisto.gamma, background.goea.star)

bp.goea.downreg.kallisto.gamma <- GOEA(bp.ant.go.hv, genelist.downreg.kallisto.gamma, background.goea.star)
cc.goea.downreg.kallisto.gamma <- GOEA(cc.ant.go.hv, genelist.downreg.kallisto.gamma, background.goea.star)
mf.goea.downreg.kallisto.gamma <- GOEA(mf.ant.go.hv, genelist.downreg.kallisto.gamma, background.goea.star)

bp.goea.upreg.kallisto.proton <- GOEA(bp.ant.go.hv, genelist.upreg.kallisto.proton, background.goea.star)
cc.goea.upreg.kallisto.proton <- GOEA(cc.ant.go.hv, genelist.upreg.kallisto.proton, background.goea.star)
mf.goea.upreg.kallisto.proton <- GOEA(mf.ant.go.hv, genelist.upreg.kallisto.proton, background.goea.star)

bp.goea.downreg.kallisto.proton <- GOEA(bp.ant.go.hv, genelist.downreg.kallisto.proton, background.goea.star)
cc.goea.downreg.kallisto.proton <- GOEA(cc.ant.go.hv, genelist.downreg.kallisto.proton, background.goea.star)
mf.goea.downreg.kallisto.proton <- GOEA(mf.ant.go.hv, genelist.downreg.kallisto.proton, background.goea.star)

####
dotplot(bp.goea.upreg.kallisto.gamma, showCategory=10)
dotplot(cc.goea.upreg.kallisto.gamma, showCategory=10)
dotplot(mf.goea.upreg.kallisto.gamma, showCategory=10)

dotplot(bp.goea.downreg.kallisto.gamma, showCategory=10)
#dotplot(cc.goea.downreg.kallisto.gamma, showCategory=10)
dotplot(mf.goea.downreg.kallisto.gamma, showCategory=10)

dotplot(bp.goea.upreg.kallisto.proton, showCategory=10)
dotplot(cc.goea.upreg.kallisto.proton, showCategory=10)
dotplot(mf.goea.upreg.kallisto.proton, showCategory=10)

#dotplot(bp.goea.downreg.kallisto.proton, showCategory=10)
#dotplot(cc.goea.downreg.kallisto.proton, showCategory=10) 
#dotplot(mf.goea.downreg.kallisto.proton, showCategory=10)

heatplot(bp.goea.upreg.kallisto.gamma, foldChange=geneFC.upreg.kallisto.gamma)
heatplot(cc.goea.upreg.kallisto.gamma, foldChange=geneFC.upreg.kallisto.gamma)
heatplot(mf.goea.upreg.kallisto.gamma, foldChange=geneFC.upreg.kallisto.gamma)

heatplot(bp.goea.downreg.kallisto.gamma, foldChange=geneFC.downreg.kallisto.gamma)
#heatplot(cc.goea.downreg.kallisto.gamma, foldChange=geneFC.downreg.kallisto.gamma)
heatplot(mf.goea.downreg.kallisto.gamma, foldChange=geneFC.downreg.kallisto.gamma)

heatplot(bp.goea.upreg.kallisto.proton, foldChange=geneFC.upreg.kallisto.proton)
heatplot(cc.goea.upreg.kallisto.proton, foldChange=geneFC.upreg.kallisto.proton)
heatplot(mf.goea.upreg.kallisto.proton, foldChange=geneFC.upreg.kallisto.proton)

#heatplot(bp.goea.downreg.kallisto.proton, foldChange=geneFC.downreg.kallisto.proton)
#heatplot(cc.goea.downreg.kallisto.proton, foldChange=geneFC.downreg.kallisto.proton)
#heatplot(mf.goea.downreg.kallisto.proton, foldChange=geneFC.downreg.kallisto.proton)

#####
p_bp <- clean_dotplot(bp.goea.upreg.kallisto.proton, "Biological Process", 5)
p_cc <- clean_dotplot(cc.goea.upreg.kallisto.proton, "Cellular Component", 5)
p_mf <- clean_dotplot(mf.goea.upreg.kallisto.proton, "Molecular Function", 5)

p_bp + p_cc + p_mf
