library(ngsReports)
library(ggplot2)

model_genome <- list.files("D:/", pattern = ".*fasta", full.names = TRUE)
model_transcriptome <- list.files("C:/Users/ПК/Documents/bachelor/data/reference_data/", pattern = ".*fasta", full.names = TRUE)

refDir_raw <- list.files("C:/Users/ПК/Documents/bachelor/data/1_QC/raw_qc/", pattern = "HV_Rep[1-3]_Ref.*fastqc.zip", full.names = TRUE)
gammaDir_raw <- list.files("C:/Users/ПК/Documents/bachelor/data/1_QC/raw_qc/", pattern = "HV_Rep[1-3]_Gamma.*fastqc.zip", full.names = TRUE)
protonDir_raw <- list.files("C:/Users/ПК/Documents/bachelor/data/1_QC/raw_qc/", pattern = "HV_Rep[1-3]_Proton.*fastqc.zip", full.names = TRUE)

refDir_fastp <- list.files("C:/Users/ПК/Documents/bachelor/data/1_QC/fastp_qc/", pattern = "HV_Rep[1-3]_Ref.*fastqc.zip", full.names = TRUE)
gammaDir_fastp <- list.files("C:/Users/ПК/Documents/bachelor/data/1_QC/fastp_qc/", pattern = "HV_Rep[1-3]_Gamma.*fastqc.zip", full.names = TRUE)
protonDir_fastp <- list.files("C:/Users/ПК/Documents/bachelor/data/1_QC/fastp_qc/", pattern = "HV_Rep[1-3]_Proton.*fastqc.zip", full.names = TRUE)

data <- FastqcDataList(protonDir_fastp)

plotSummary(data)
plotSeqContent(data[1], plotType = "line")
plotBaseQuals(data[1], plotType = "boxplot")
plotSeqQuals(data[1], plotType = "line")
plotAdapterContent(data[1], plotType = "line")
plotDupLevels(data[1], plotType = "line")
plotGcContent(data, Fastafile = model_genome, plotType = "line")
plotGcContent(data, Fastafile = model_transcriptome, plotType = "line")
plotOverrep(data)
