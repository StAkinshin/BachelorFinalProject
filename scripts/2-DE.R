# Add libraries ------------------------------

library(tximport)
library(rhdf5)
library(DESeq2)
library(ggplot2)
library(VennDiagram)

# Define variables ---------------------------

project.metadata <- read.csv(
  file = "D:/ASD_BachelorFinalProject/BachelorFinalProject/project-metadata.csv", 
  row.names = 1,
  sep = ';',
  header = F)

conditions <- c("gamma", "proton")
aligners <- c("STAR", "kallisto")

# Check correlation in counts' table ---------

# Generate DE Matrix  ----------------

generate_matrix(aligner, condition) {
  
  if aligner == "STAR" {
    
    counts.files <- list.files(counts.folder.path,
                              pattern = '*.txt',
                              full.names = T
    )
    counts.list <- lapply(counts.files, 
                          function(file) read.table(file,
                                                    header = T,
                                                    row.names = 1,
                                                    na.strings = 'NA',
                                                    dec = '.'))
    matrix <- do.call(cbind, counts.list)
    colnames(matrix) <- rownames(metadata)
    
    return(matrix)
  }
  
  if aligner == "kallisto" {
    
   # SampleNames <- rownames(metadata.gamma)
   # KallistoFiles <- file.path('C://Users/ПК/Documents/bachelor/data/2_Aligning/kallisto/', rownames(metadata.gamma), "abundance.h5")
    # names(KallistoFiles) <- rownames(metadata.gamma)
    
    #TXI <- tximport(KallistoFiles,
    #                type = "kallisto",
    #                tx2gene = tx2gene,
    #                ignoreTxVersion = TRUE)
    
    #Matrix.kallisto <- as.data.frame(TXI$counts)
    
  }
  
  else print "There's no aligner with such a name, sorry!"
  
}


# Check PCA for Matrix -----------------------

# Perform DE analysis ------------------------

# Plot DE results ----------------------------
