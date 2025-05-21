# Add libraries ------------------------------

library(ngsReports)
library(ggplot2)

# Define variables -------------

project.metadata <- read.csv("D:/ASD_BachelorFinalProject/BachelorFinalProject/project-metadata.csv", 
                             row.names = 1,
                             sep = ';',
                             header = F)

ref.genome_fasta <- project.metadata[1, 1]
ref.transcriptome_fasta <- project.metadata[3, 1]

qc.data.folder_raw <- project.metadata[9, 1]
qc.data.folder_trim <- project.metadata[10, 1]

conditions <- c("Ref", "Gamma", "Proton") 
  
# QC -----------------------------------------

condition <- "Proton"

sample.pattern <- paste0("HV_Rep[1-3]_", condition, ".*fastqc.zip")
print(sample.pattern)

qc.data.name_raw <- paste0(project.metadata[14, 1], "HV_Rep[1-3]_", condition, "_QC-data_raw.RData")
qc.data.name_trim <- paste0(project.metadata[14, 1], "HV_Rep[1-3]_", condition, "_QC-data_trim.RData")

for (condition in conditions) {
  sample.pattern <- paste0("HV_Rep[1-3]_", condition, ".*fastqc.zip")
  
  qc.data.filename_raw <- paste0(project.metadata[14, 1], "HV_Rep[1-3]_", condition, "_QC-data_raw.RData")
  qc.data.filename_trim <- paste0(project.metadata[14, 1], "HV_Rep[1-3]_", condition, "_QC-data_trim.RData")
  
  data_raw <- list.files(
    qc.data.folder_raw,
    pattern = sample.pattern,
    full.names = TRUE
  )
  
  qc.data_raw <- FastqcDataList(data_raw)
  save(qc.data_raw, file = qc.data.filename_raw)
  
  data_trim <- list.files(
    qc.data.folder_trim,
    pattern = sample.pattern,
    full.names = TRUE
  )
  
  qc.data_trim <- FastqcDataList(data_trim)
  save(qc.data_trim, file = qc.data.filename_trim)
}

# Plot QC results ----------------------------

#for (condition in conditions) {
  qc.data_raw <- paste0(project.metadata[14, 1], "HV_Rep[1-3]_", condition, "_QC-data_raw.RData")
  qc.data_trim <- paste0(project.metadata[14, 1], "HV_Rep[1-3]_", condition, "_QC-data_trim.RData")
  
  load(qc.data_raw)
  plotSummary(qc.data_raw)
  plotSeqContent(qc.data_raw[1:6], plotType = "line", )
  plotBaseQuals(qc.data_raw[1:6], plotType = "boxplot")
  plotSeqQuals(qc.data_raw[1:6], plotType = "line")
  plotAdapterContent(qc.data_raw[1:6], plotType = "line")
  plotDupLevels(qc.data_raw[1:6], plotType = "line")
  plotGcContent(qc.data_raw, Fastafile = ref.genome_fasta, plotType = "line")
  plotGcContent(qc.data_raw, Fastafile = ref.transcriptome_fasta, plotType = "line")
  plotOverrep(qc.data_raw)
  
  load(qc.data_trim)
  plotSummary(qc.data_trim)
  plotSeqContent(qc.data_trim[1:6], plotType = "line", )
  plotBaseQuals(qc.data_trim[1:6], plotType = "boxplot")
  plotSeqQuals(qc.data_trim[1:6], plotType = "line")
  plotAdapterContent(qc.data_trim[1:6], plotType = "line")
  plotDupLevels(qc.data_trim[1:6], plotType = "line")
  plotGcContent(qc.data_trim, Fastafile = ref.genome_fasta, plotType = "line")
  plotGcContent(qc.data_trim, Fastafile = ref.transcriptome_fasta, plotType = "line")
  plotOverrep(qc.data_trim)
}

