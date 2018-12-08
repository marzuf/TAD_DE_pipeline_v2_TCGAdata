startTime <- Sys.time()
cat(paste0("> Rscript prep_settingFile.R\n"))

# Rscript prep_settingFile.R

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/TAD_DE_pipeline_v2_TCGAdata")
source("analysis_utils.R")

outFold <- file.path("PREP_DATA")
system(paste0("mkdir -p ", outFold))

cancerType <- "CESC"
cond1 <- "AdenoCarcinoma"
cond2 <- "SquamousCarcinoma"
type_out <- tolower(cancerType)
cond1_out <- "adeno"
cond2_out <- "squamous"

mutType <- "symbol"

# Rscript prep_data.R <cancer type> <cond1> <cond1_out> <cond2> <cond2_out>
# Rscript prep_data.R CESC AdenoCarcinoma adeno SquamousCarcinoma squamous
# Rscript prep_data.R PAAD paad wt mutAMP mut
# Rscript prep_data.R STAD norm norm GS gs  

checkFile <- "DATA_SUMMARY/tcga_datasets_DT.Rdata"
checkDT <- eval(parse(text = load(checkFile)))
#type subtype      class nSamplesAnnot nSamplesAnnotExpr nSamplesAnnotExprMut
#1  ACC     ACC        ACC            76                76                   76

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
cancerType <- args[1]
cond1 <- args[2]
cond1_out <- args[3]
cond2 <- args[4]
cond2_out <- args[5]
type_out <- tolower(cancerType)

cmp_name <- paste0("TCGA", type_out, "_", cond1_out, "_", cond2_out)

logFile <- file.path(outFold, paste0(cmp_name, "_prep_data_logFile.txt"))
system(paste0("rm -f ", logFile))


## if mut cond2 <- "mutKRAAS"

if(cond1 == "norm") {
  nCheck1 <- NULL
  nCheck2 <- checkDT$nSamplesAnnotExpr[checkDT$type == cancerType & checkDT$subtype == cond2]
} else if(grepl("^mut.+", cond2)) {
  mymut <- gsub("^mut(.+)$", "\\1", cond2)
  stopifnot(length(mymut) == 1)
  mymut <- as.character(mymut)
  nCheck1 <- checkDT$nSamplesAnnotExprMut[checkDT$type == cancerType & checkDT$subtype == cond1]
  nCheck2 <- checkDT$nSamplesAnnotExprMut[checkDT$type == cancerType & checkDT$subtype == cond2]
} else {
  nCheck1 <- checkDT$nSamplesAnnotExpr[checkDT$type == cancerType & checkDT$subtype == cond1]
  nCheck2 <- checkDT$nSamplesAnnotExpr[checkDT$type == cancerType & checkDT$subtype == cond2]
}

# /mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain/[...].rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
exprFolder <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain")
annotFile <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_sample_annotation/tcga.pancanAtlas_sample_type_annotation_v4.txt")
# firstly, run using v4 (_GAMv4 folders)
#gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v4.RData")
# then, update: e-mail from Marco 29.11.2018
gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v5.RData")

gamEntrez2symb <- eval(parse(text = load(file.path("GAM_IDs_ENTREZ", "entrezID_2_symbol.Rdata"))))
# 57801   388585     
# "HES4"   "HES5" 
gamAlterationsByEntrez <- eval(parse(text = load(file.path("GAM_IDs_ENTREZ", "alterations_by_entrezID.Rdata"))))
# $`57801`
# [1] "DEL.HES4"

txt <- paste0("... using exprFolder\t=\t", exprFolder, "\n")
printAndLog(txt, logFile)
txt <- paste0("... using annotFile\t=\t", annotFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... using gamFile\t=\t", gamFile, "\n")
printAndLog(txt, logFile)

outputExprFile <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name, "rnaseqDT_v2.Rdata")
outputSamp1File <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name,  paste0(cond1_out, "_ID.Rdata"))
outputSamp2File <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name, paste0(cond2_out, "_ID.Rdata"))

#cat(outputExprFile, "\n")
#stopifnot(!file.exists(outputExprFile))
#stopifnot(!file.exists(outputSamp1File))
#stopifnot(!file.exists(outputSamp2File))

system(paste0("mkdir -p ", dirname(outputExprFile)))
system(paste0("mkdir -p ", dirname(outputSamp1File)))
system(paste0("mkdir -p ", dirname(outputSamp2File)))

all_exprFiles <- list.files(exprFolder, full.names = TRUE, pattern=paste0(".+illuminahiseq.+rnaseqv2.+RSEM.+"))

annotDT <- read.delim(annotFile, stringsAsFactors = FALSE)


if(grepl("_mut", cmp_name)) {
  
  cat("... load GAM data\n")
  x <- load(gamFile)
  stopifnot("pancan_gam" %in% x)
  gamDT <- pancan_gam
  gamDT[1:3,1:3]
  gam_samples <- rownames(gamDT)

  all_samp_ID <- as.character(annotDT$sample_id[annotDT$type == cancerType])
  
  all_gam_samp_ID <- all_samp_ID[all_samp_ID %in% gam_samples]
  stopifnot(length(all_gam_samp_ID) > 0)
  
  txt <- paste0("# annotated samples for ", cancerType, " ", cond1, " with gamData\t=\t", length(all_gam_samp_ID), "\n")
  printAndLog(txt, logFile)

  curr_gamDT <- gamDT[rownames(gamDT) %in% all_gam_samp_ID,]
  stopifnot(nrow(curr_gamDT) == length(all_gam_samp_ID))
  
  curr_gamDT <- t(curr_gamDT) # now put the samples as columns
  curr_gamDT[1:3,1:3]
  
  # select the alterations for the current mutation
  txt <- paste0("... for alterations in\t=\t", mymut, "\n" )
  printAndLog(txt, logFile)
  
  if(mutType == "symbol") {
    stopifnot(mymut %in% gamEntrez2symb)
    mymut_entrez <- names(gamEntrez2symb[gamEntrez2symb == mymut])
    stopifnot(length(mymut_entrez) == 1)
  } else if(mutType == "entrez") { 
    stopifnot(mymut %in% names(gamEntrez2symb))
    mymut_entrez <- mymut
  } else {
    stop("-- error -- \n")
  }
  stopifnot(mymut_entrez %in% names(gamAlterationsByEntrez))
  myalterations <- gamAlterationsByEntrez[[paste0(mymut_entrez)]]
  
  txt <- paste0("... found following alterations:\t=\t", paste0(myalterations, collapse=","), "\n")
  printAndLog(txt, logFile)
  
  stopifnot(myalterations %in% rownames(curr_gamDT))
  
  curr_gamDT <- curr_gamDT[myalterations,]
  stopifnot(nrow(curr_gamDT) == length(myalterations))
  
  altBySamp <- colSums(curr_gamDT)
  
  # wt (cond1) or mut (cond2)
  all_cond1_ID <- names(altBySamp[altBySamp == 0])
  all_cond2_ID <- names(altBySamp[altBySamp > 0])
    
  stopifnot(all_cond1_ID %in% all_gam_samp_ID)
  stopifnot(all_cond2_ID %in% all_gam_samp_ID)
  
  ### /// select all_cond1_ID and all_cond2_ID according to the mutation
  
  ### // select the row with the given mutation in gamDT
  
  ### by definition: cond1 = wt and cond2 = mut
  
  
  
} else if(cond1 == "norm") {
  all_cond2_ID <- as.character(annotDT$sample_id[annotDT$type == cancerType & annotDT$subtype %in% c(cond2)])
  txt <- paste0("# annotated samples for ", cancerType, " ", cond2, "\t=\t", length(all_cond2_ID), "\n")
  printAndLog(txt, logFile)
} else {
  all_cond1_ID <- as.character(annotDT$sample_id[annotDT$type == cancerType & annotDT$subtype %in% c(cond1)])
  all_cond2_ID <- as.character(annotDT$sample_id[annotDT$type == cancerType & annotDT$subtype %in% c(cond2)])
  
  txt <- paste0("# annotated samples for ", cancerType, " ", cond1, "\t=\t", length(all_cond1_ID), "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("# annotated samples for ", cancerType, " ", cond2, "\t=\t", length(all_cond2_ID), "\n")
  printAndLog(txt, logFile)
}

exprFile <- all_exprFiles[grepl(cancerType, basename(all_exprFiles))]
stopifnot(length(exprFile) == 1)
cat("... load expression data\n")
exprDT <- read.delim(exprFile, stringsAsFactors = FALSE)
colnames(exprDT) <- gsub("\\.", "-", colnames(exprDT))
exprDT[1:3,1:3]

if(cond1 == "norm") {
  samples <- colnames(exprDT)
  stopifnot(samples[1] == "Hybridization-REF")
  samples <- samples[-1]
  stopifnot(samples[1] != "Hybridization-REF")
  sampTypes <- substr(x=samples, start=14, stop=15)
  all_cond1_ID <- samples[which(sampTypes == "11")]
  stopifnot(grepl("TCGA-..-....-11.+", all_cond1_ID))
  all_cond1_ID <- substr(x = all_cond1_ID, start = 1, stop = 15)
}
stopifnot(length(all_cond1_ID) > 0)
stopifnot(length(all_cond2_ID) > 0)

stopifnot(colnames(exprDT)[1] == "Hybridization-REF")
stopifnot(exprDT[which(exprDT[,1] == "gene_id"), 2:ncol(exprDT)] == "normalized_count")
exprDT <- exprDT[-which(exprDT[,1] == "gene_id"),]
exprDT[1:3,1:3]
geneIDs <- strsplit(as.character(exprDT[,"Hybridization-REF"]), split="\\|")
geneIDs <- sapply(geneIDs, function(x) x[[2]])
stopifnot(!duplicated(geneIDs))
# rownames(exprDT) <- exprDT[,"Hybridization-REF"]
rownames(exprDT) <- geneIDs
exprDT <- exprDT[,-which(colnames(exprDT) == "Hybridization-REF")]
exprDT[1:3,1:3]

tmpDT <- as.data.frame(data.matrix(exprDT))
stopifnot(as.character(tmpDT[4,5]) == as.character(exprDT[4,5]))
stopifnot(is.numeric(tmpDT[4,5]))
exprDT[1:3,1:3]
tmpDT[1:3,1:3]
exprDT <- tmpDT

# expr_samples <- colnames(exprDT)
# expr_samples <- substr(x = expr_samples, start = 1, stop = 15)
colnames(exprDT) <- substr(x = colnames(exprDT), start = 1, stop = 15)
expr_samples <- colnames(exprDT)

stopifnot(!duplicated(expr_samples))
stopifnot(any(all_cond1_ID %in% expr_samples))
stopifnot(any(all_cond2_ID %in% expr_samples))

cond1_ID <- all_cond1_ID[all_cond1_ID %in% expr_samples]
stopifnot(cond1_ID %in% colnames(exprDT))

txt <- paste0("# annotated samples for ", cancerType, " ", cond1, " with exprData\t=\t", length(cond1_ID), "\n")
printAndLog(txt, logFile)

cond2_ID <- all_cond2_ID[all_cond2_ID %in% expr_samples]
stopifnot(cond2_ID %in% colnames(exprDT))

txt <- paste0("# annotated samples for ", cancerType, " ", cond2, " with exprData\t=\t", length(cond2_ID), "\n")
printAndLog(txt, logFile)

exprDT <- exprDT[, c(cond1_ID, cond2_ID)]
stopifnot(ncol(exprDT) == (length(cond1_ID) + length(cond2_ID)))

exprDT[1:3,1:3]

#stopifnot(!file.exists(outputExprFile))
#stopifnot(!file.exists(outputSamp1File))
#stopifnot(!file.exists(outputSamp2File))


save(exprDT, file = outputExprFile)
txt <- paste0("... written:", outputExprFile, "\n")
printAndLog(txt, logFile)

save(cond1_ID, file = outputSamp1File)
txt <- paste0("... written:", outputSamp1File, "\n")
printAndLog(txt, logFile)

save(cond2_ID, file = outputSamp2File)
txt <- paste0("... written:", outputSamp2File, "\n")
printAndLog(txt, logFile)

stopifnot(length(cond2_ID) == nCheck2)

if(cond1 != "norm"){
  stopifnot(length(cond1_ID) == nCheck1)  
}



######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



