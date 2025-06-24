library(ggplot2, ggdendro)

ExprDataImport <- function(STARCountsFolderPath, MetaData) {
  
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
  return(Matrix)
  }
 
exprdata.star.gamma <- ExprDataImport('C://Users/ПК/Documents/bachelor/data/2_Aligning/star+featureCounts/counts_gamma/', 
                                      metadata.gamma)
exprdata.star.proton <- ExprDataImport('C://Users/ПК/Documents/bachelor/data/2_Aligning/star+featureCounts/counts_proton',
                                       metadata.proton)

exprdata.star.merged <- merge(exprdata.star.gamma, exprdata.star.proton, by = "row.names", all = T)
exprdata.star.merged <- exprdata.star.merged[, -1] 

cor.star.gamma <- cor(exprdata.star.gamma)
cor.star.proton <- cor(exprdata.star.proton)
cor.star.merged <- cor(exprdata.star.merged)

distance.star.gamma <- as.dist(1 - cor.star.gamma)
distance.star.proton <- as.dist(1 - cor.star.proton)
distance.star.merged <- as.dist(1 - cor.star.merged)

tree.star.gamma <- hclust(distance.star.gamma)
tree.star.proton <- hclust(distance.star.proton)
tree.star.merged <- hclust(distance.star.merged)

ggdendro::ggdendrogram(tree.star.gamma, rotate = T)
ggdendro::ggdendrogram(tree.star.proton, rotate = T)
ggdendro::ggdendrogram(tree.star.merged)



