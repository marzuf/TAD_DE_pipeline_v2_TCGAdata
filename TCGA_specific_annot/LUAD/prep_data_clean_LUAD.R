startTime <- Sys.time()
cat(paste0("> Rscript prep_data_clean_LUAD.R\n"))

# Rscript prep_data_clean_LUAD.R



options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/TAD_DE_pipeline_v2_TCGAdata")
source("../../analysis_utils.R")

outFold <- file.path("PREP_DATA_VALL")
system(paste0("mkdir -p ", outFold))

cancerType <- "LUAD"
cond1 <- "nonsmoker"
cond2 <- "smoker"
type_out <- tolower(cancerType)
cond1_out <- "nonsmoker"
cond2_out <- "smoker"


all_cond1_ID_file <- paste0(cond1_out, "_ID.Rdata")
stopifnot(file.exists(all_cond1_ID_file))
all_cond1_ID <- eval(parse(text = load(all_cond1_ID_file)))

all_cond2_ID_file <- paste0(cond2_out, "_ID.Rdata")
stopifnot(file.exists(all_cond2_ID_file))
all_cond2_ID <- eval(parse(text = load(all_cond2_ID_file)))



cmp_name <- paste0("TCGA", type_out, "_", cond1_out, "_", cond2_out)

logFile <- file.path(outFold, paste0(cmp_name, "_prep_data_logFile.txt"))
system(paste0("rm -f ", logFile))



### ADDED 07.12.2018

# scaled estimates folder
seFolderMain <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_gdac_full_snapshots", cancerType)
seSubFolders <- list.files(seFolderMain, pattern="^2")
# if("20160715" %in% seSubFolders) {
#   data_release <- "20160715"
# } else if("20160128" %in% seSubFolders) {
#   data_release <- "20160128"
# } else{
#   stop("*** ERROR ***\n")
# }
data_release <- "20160128"
seFolder <- file.path(seFolderMain, data_release)
seFile <- list.files(seFolder, full.names = TRUE, pattern=paste0(".+illuminahiseq.+rnaseqv2.+RSEM_genes__data.Level_3.+gz$"))
stopifnot(length(seFile) == 1)

# /mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain/[...].rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
exprFolder <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain")
annotFile <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_sample_annotation/tcga.pancanAtlas_sample_type_annotation_v4.txt")
# firstly, run using v4 (_GAMv4 folders)
#gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v4.RData")
# then, update: e-mail from Marco 29.11.2018
gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v5.RData")

gamEntrez2symb <- eval(parse(text = load(file.path("../../GAM_IDs_ENTREZ", "entrezID_2_symbol.Rdata"))))
# 57801   388585     
# "HES4"   "HES5" 
gamAlterationsByEntrez <- eval(parse(text = load(file.path("../../GAM_IDs_ENTREZ", "alterations_by_entrezID.Rdata"))))
# $`57801`
# [1] "DEL.HES4"

txt <- paste0("... using exprFolder\t=\t", exprFolder, "\n")
printAndLog(txt, logFile)
txt <- paste0("... using annotFile\t=\t", annotFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... using gamFile\t=\t", gamFile, "\n")
printAndLog(txt, logFile)

outputFPKMFile <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name, "fpkmDT.Rdata")
outputExprFile <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name, "rnaseqDT_v2.Rdata")
outputSamp1File <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name,  paste0(cond1_out, "_ID.Rdata"))
outputSamp2File <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name, paste0(cond2_out, "_ID.Rdata"))
cat(outputExprFile, "\n")
stopifnot(!file.exists(outputExprFile))
stopifnot(!file.exists(outputSamp1File))
stopifnot(!file.exists(outputSamp2File))


#cat(outputExprFile, "\n")
#stopifnot(!file.exists(outputExprFile))
#stopifnot(!file.exists(outputSamp1File))
#stopifnot(!file.exists(outputSamp2File))

system(paste0("mkdir -p ", dirname(outputExprFile)))
system(paste0("mkdir -p ", dirname(outputSamp1File)))
system(paste0("mkdir -p ", dirname(outputSamp2File)))

all_exprFiles <- list.files(exprFolder, full.names = TRUE, pattern=paste0(".+illuminahiseq.+rnaseqv2.+RSEM.+"))

annotDT <- read.delim(annotFile, stringsAsFactors = FALSE)


exprFile <- all_exprFiles[grepl(cancerType, basename(all_exprFiles))]
stopifnot(length(exprFile) == 1)



cat("... load expression data\n")
exprDT <- read.delim(exprFile, stringsAsFactors = FALSE)
colnames(exprDT) <- gsub("\\.", "-", colnames(exprDT))
exprDT[1:3,1:3]

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
cat("exprDT[4,5] = ", exprDT[4,5], "\n")
cat("tmpDT[4,5] = ", tmpDT[4,5], "\n")
stopifnot(as.numeric(as.character(tmpDT[4,5])) == as.numeric(as.character(exprDT[4,5])))
stopifnot(is.numeric(tmpDT[4,5]))
exprDT[1:3,1:3]
tmpDT[1:3,1:3]
exprDT <- tmpDT

# expr_samples <- colnames(exprDT)
# expr_samples <- substr(x = expr_samples, start = 1, stop = 15)
colnames(exprDT) <- substr(x = colnames(exprDT), start = 1, stop = 15)
expr_samples <- colnames(exprDT)

expr_samples_short <- substr(x = expr_samples, start = 1, stop = 12)

#stopifnot(!duplicated(expr_samples_short)) # not true
stopifnot(any(all_cond1_ID %in% expr_samples_short))
stopifnot(any(all_cond2_ID %in% expr_samples_short))


cond1_ID <- expr_samples[expr_samples_short %in% all_cond1_ID]
stopifnot(cond1_ID %in% colnames(exprDT))

txt <- paste0("# annotated samples for ", cancerType, " ", cond1, " with exprData\t=\t", length(cond1_ID), "\n")
printAndLog(txt, logFile)

cond2_ID <- expr_samples[expr_samples_short %in% all_cond2_ID]
stopifnot(cond2_ID %in% colnames(exprDT))

txt <- paste0("# annotated samples for ", cancerType, " ", cond2, " with exprData\t=\t", length(cond2_ID), "\n")
printAndLog(txt, logFile)

exprDT <- exprDT[, c(cond1_ID, cond2_ID)]
stopifnot(ncol(exprDT) == (length(cond1_ID) + length(cond2_ID)))

exprDT[1:3,1:3]

save(exprDT, file = outputExprFile)
txt <- paste0("... written:", outputExprFile, "\n")
printAndLog(txt, logFile)

save(cond1_ID, file = outputSamp1File)
txt <- paste0("... written:", outputSamp1File, "\n")
printAndLog(txt, logFile)

save(cond2_ID, file = outputSamp2File)
txt <- paste0("... written:", outputSamp2File, "\n")
printAndLog(txt, logFile)


cat("... length(cond1_ID) = ", length(cond1_ID), "\n")
#cat("... nCheck1 = ", nCheck1, "\n")
cat("... length(cond2_ID) = ", length(cond2_ID), "\n")
#cat("... nCheck2 = ", nCheck2, "\n")

#######################################################################################  PREP THE FPKM DATA
# discussion with Marco 07.12.2018:
# to have similar data like FPKM:
# take the TCGA scaled estimates file
# the TPM = EC/l * 10⁶/# reads [EC = expected counts, l = length]
# the scaled estimates = EC/l * 1/reads
# the EC are the raw counts divided by the 75% value, not corrected for gene length
# for limma: better to use the EC
# for comparability with FPKM:
# scaled estimates * 10^6
# (TPM: sum of all TPMs in each sample can be different; RPKM: sum of normalized reads may be different)

cat("... start preparing TPM data for RSEM\n")

cmd <- paste0("tar -xzvf ", seFile)
cat(paste0(cmd, "\n"))
system(cmd)
stopifnot(file.exists(seFile))

extractFolder <- gsub(".tar.gz$", "", basename(seFile))
seInFile <- file.path(extractFolder, gsub(paste0("^.+", cancerType, ".Merge_"), paste0(cancerType, "."), basename(seFile)))
seInFile <- file.path(extractFolder, gsub("(.+__data\\.)Level_3.+", "\\1data.txt", basename(seInFile)))

seDT <- read.delim(seInFile, stringsAsFactors = FALSE)
colnames(seDT) <- gsub("\\.", "-", colnames(seDT))
stopifnot(colnames(seDT)[1] == "Hybridization-REF")
stopifnot(seDT[1,1] == "gene_id")

colsToKeep <- seDT[1,] %in% c("gene_id", "scaled_estimate")
stopifnot(length(colsToKeep) == ncol(seDT))
seDT <- seDT[, colsToKeep]
stopifnot(ncol(seDT) == sum(colsToKeep))
stopifnot(seDT[which(seDT[,1] == "gene_id"), 2:ncol(seDT)] == "scaled_estimate")
seDT <- seDT[-which(seDT[,1] == "gene_id"),]
seDT[1:3,1:3]
geneIDs <- strsplit(as.character(seDT[,"Hybridization-REF"]), split="\\|")
geneIDs <- sapply(geneIDs, function(x) x[[2]])
stopifnot(!duplicated(geneIDs))
# rownames(seDT) <- seDT[,"Hybridization-REF"]
rownames(seDT) <- geneIDs
seDT <- seDT[,-which(colnames(seDT) == "Hybridization-REF")]
seDT[1:3,1:3]

tmpDT <- as.data.frame(data.matrix(seDT))
cat("seDT[4,5] = ", seDT[4,5], "\n")
cat("tmpDT[4,5] = ", tmpDT[4,5], "\n")
stopifnot(as.numeric(as.character(tmpDT[4,5])) == as.numeric(as.character(seDT[4,5])))
stopifnot(is.numeric(tmpDT[4,5]))
seDT[1:3,1:3]
tmpDT[1:3,1:3]
seDT <- tmpDT

colnames(seDT) <- substr(x = colnames(seDT), start = 1, stop = 15)
stopifnot(cond1_ID %in% colnames(seDT))
stopifnot(cond2_ID %in% colnames(seDT))

stopifnot(rownames(exprDT) %in% rownames(seDT))
seDT <- seDT[rownames(exprDT),c(cond1_ID, cond2_ID)]

stopifnot(dim(seDT) == dim(exprDT))

stopifnot(file.exists(seInFile))
stopifnot(file.exists(seFile))
stopifnot(grepl(".tar.gz$", seFile))  # the file should always exist !
stopifnot(grepl(".txt$", seInFile))  # the file should always exist !
folderToRemove <- dirname(seInFile)
stopifnot(folderToRemove != seFile)

cmd <- paste0("rm -rf ", folderToRemove)
cat(paste0(cmd, "\n"))
system(cmd)

stopifnot(file.exists(seFile))
stopifnot(! file.exists(seInFile))
stopifnot(! file.exists(folderToRemove))

save(seDT, file = outputFPKMFile)
txt <- paste0("... written:", outputFPKMFile, "\n")
printAndLog(txt, logFile)

#######################################################################################


######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



