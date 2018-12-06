startTime <- Sys.time()

library(foreach)
library(doMC)

cat("> START: data_summary.R\n")
# Rscript data_summary.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/TAD_DE_pipeline_v2_TCGAdata")

registerDoMC(ifelse(SSHFS, 2, 40))

source("analysis_utils.R")

outFold <- "DATA_SUMMARY"
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

#diff DATA_SUMMARY/tcga_datasets_DT.txt DATA_SUMMARY_GAMv4/tcga_datasets_DT.txt
# -> no diff

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

txt <- paste0("# cohorts expression data\t=\t", length(expr_cohorts), "\n")
printAndLog(txt, logFile)


##### ANNOTATION DATA
cat("... load annotation data\n")
annotDT <- read.delim(annotFile, header=T, stringsAsFactors = FALSE)
table(annotDT$type)

annot_cohorts <- unique(annotDT$type)
txt <- paste0("# cohorts annotation data\t=\t", length(annot_cohorts), "\n")
printAndLog(txt, logFile)

expr_cohorts[!expr_cohorts %in% annot_cohorts] # "COADREAD" "KIPAN"    "STES"    
annot_cohorts[!annot_cohorts %in% expr_cohorts] # "COAD" "READ"


txt <- paste0("!!! remove rows with missing information !!!")
printAndLog(txt, logFile)
txt <- paste0("... nrow before na.omit\t=\t", nrow(annotDT),"\n")
printAndLog(txt, logFile)
annotDT <- na.omit(annotDT)
txt <- paste0("... nrow after na.omit\t=\t", nrow(annotDT),"\n")
printAndLog(txt, logFile)

#take only tumor samples
sampleNbr <- as.numeric(gsub(".+-.+-.+-(.+)", "\\1", annotDT$sample_id))
stopifnot(!is.na(sampleNbr))
stopifnot(sampleNbr < 10)

annotDT$type_subtype <- ifelse(annotDT$type == annotDT$subtype, annotDT$type, paste0(annotDT$type, "_", annotDT$subtype))

stopifnot(all(annotDT$type_subtype == annotDT$class))

annotDT$type_subtype <- NULL
table(annotDT$class)

stopifnot(!duplicated(annotDT$sample_id))
stopifnot(!duplicated(annotDT$patient_id))

# trial0 <- setNames(aggregate(sample_id ~ type, FUN=length, data=annotDT), c("type", "nSamples"))
# trial1 <- aggregate(list(nSamples=annotDT$sample_id), list(type=annotDT$type), FUN= length) 
dt1 <- setNames(aggregate(sample_id ~ type, FUN=length, data=annotDT), c("type", "nSamples"))
dt2 <- setNames(aggregate(sample_id ~ type+subtype, FUN=length, data=annotDT), c("type", "subtype", "nSamples"))
dt3 <- setNames(aggregate(sample_id ~ type+subtype+class, FUN=length, data=annotDT), c("type", "subtype", "class", "nSamples"))
dt4 <- setNames(aggregate(sample_id ~ class, FUN=length, data=annotDT), c("type", "class"))
stopifnot(nrow(dt3) == nrow(dt4))
stopifnot(nrow(dt3) == nrow(dt2))

summaryDT <- setNames(aggregate(sample_id ~ type+subtype+class, FUN=length, data=annotDT), c("type", "subtype", "class", "nSamples"))

all_types <- unique(annotDT$type)

##### GAM MATRIX
cat("... load GAM data\n")
x <- load(gamFile)
stopifnot("pancan_gam" %in% x)
gamDT <- pancan_gam
gam_samples <- gsub("-", ".", rownames(gamDT))

################################################################################################################################################

i=1
tcga_datasets_DT <- foreach(i = seq_len(nrow(summaryDT)), .combine='rbind') %dopar% {
  
  curr_type <- summaryDT$type[i]
  curr_subtype <- summaryDT$subtype[i]
  curr_class <- summaryDT$class[i]

  cat("start for\t", curr_type, " - ", curr_subtype, "\n")
    
  curr_annotDT <- annotDT[annotDT$type == curr_type & annotDT$subtype == curr_subtype & annotDT$class == curr_class,]
  
  ### check  expression data available
  cat("... load expression data\n")
  exprFile <- all_exprFiles[grepl(curr_type, basename(all_exprFiles))]
  stopifnot(length(exprFile) == 1)
  exprDT <- read.delim(exprFile)
  annot_samples <- gsub("-", ".", curr_annotDT$sample_id)
  stopifnot(colnames(exprDT)[1] == "Hybridization.REF")
  expr_samples <- colnames(exprDT)[2:ncol(exprDT)]
  expr_samples <- substr(x = expr_samples, start = 1, stop = 15)
  all(annot_samples %in% expr_samples)
  
  ### check mutation data available
  sum(annot_samples %in% gam_samples)
  
  nAnnot <- summaryDT$nSamples[i]
  nAnnotExpr <- sum(annot_samples %in% expr_samples)
  nAnnotExprMut <- sum(annot_samples %in% expr_samples & annot_samples %in% gam_samples)

  stopifnot(length(annot_samples) == nAnnot)
    
  data.frame(
    type = curr_type,
    subtype = curr_subtype,
    class = curr_class,
    nSamplesAnnot = nAnnot,
    nSamplesAnnotExpr = nAnnotExpr,
    nSamplesAnnotExprMut = nAnnotExprMut,
    stringsAsFactors = FALSE
  )
  
}

outFile <- file.path(outFold, "tcga_datasets_DT.Rdata")
save(tcga_datasets_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "tcga_datasets_DT.txt")
write.table(tcga_datasets_DT, col.names=T, row.names=F, sep="\t", quote=F, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

################################################################################################################################################



######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

#dt = read.delim("/mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain/STAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=T)
#samples = colnames(dt)
#unique(substr(x=samples, start=14, stop=15))




