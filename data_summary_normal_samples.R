startTime <- Sys.time()

library(foreach)
library(doMC)

cat("> START: data_summary.R\n")
# Rscript data_summary_normal_samples.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/TAD_DE_pipeline_v2_TCGAdata")

registerDoMC(ifelse(SSHFS, 2, 40))

source("analysis_utils.R")

outFold <- "DATA_SUMMARY_NORMAL"
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "data_summary_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""


# /mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain/[...].rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
exprFolder <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain")

annotFile <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_sample_annotation/tcga.pancanAtlas_sample_type_annotation_v4.txt")

# firstly, run using v4 (_GAMv4 folders)
#gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v4.RData")
# then, update: e-mail from Marco 29.11.2018
gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v5.RData")

txt <- paste0("exprFolder\t=\t", exprFolder, "\n" )
printAndLog(txt, logFile)
txt <- paste0("annotFile\t=\t", annotFile, "\n" )
printAndLog(txt, logFile)
txt <- paste0("gamFile\t=\t", gamFile, "\n" )
printAndLog(txt, logFile)

stopifnot(file.exists(exprFolder))
stopifnot(file.exists(annotFile))
stopifnot(file.exists(gamFile))

##### EXPRESSION DATA:
all_exprFiles <- list.files(exprFolder, full.names = TRUE, pattern=paste0(".+illuminahiseq.+rnaseqv2.+RSEM.+"))
expr_cohorts <- gsub("^(.+)\\.rna.+", "\\1", basename(all_exprFiles))



nSamp_DT <- foreach(i = seq_along(all_exprFiles), .combine='rbind') %dopar% {
  
exprFile <- all_exprFiles[i]
cohort <- expr_cohorts[i]

cat("... start ", cohort, "\n")
  
stopifnot(file.exists(exprFile))

dt <- read.delim(exprFile, header=T)
samples <- colnames(dt)

stopifnot(samples[1] == "Hybridization.REF")
samples <- samples[-1]
stopifnot(samples[1] != "Hybridization.REF")

sampTypes <- substr(x=samples, start=14, stop=15)

nNormal <- sum(as.character(sampTypes) == "11")

data.frame(cohort = cohort, 
           nNormal = nNormal,
           stringsAsFactors = FALSE)
  
}

outFile <- file.path(outFold, "nSamp_DT.RData")
save(nSamp_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "nSamp_DT.txt")
write.table(nSamp_DT, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

#dt = read.delim("/mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain/STAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=T)
#samples = colnames(dt)
#unique(substr(x=samples, start=14, stop=15))




