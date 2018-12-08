startTime <- Sys.time()
cat(paste0("> Rscript prep_SKCM_infiltrated_noninfiltrated.R\n"))

# Rscript prep_SKCM_infiltrated_noninfiltrated.R

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/TAD_DE_pipeline_v2_TCGAdata")
source("../../analysis_utils.R")

# outFold <- file.path("PREP_DATA_VALL")
# system(paste0("mkdir -p ", outFold))

cancerType <- "SKCM"
cond1 <- "infiltrated"
cond2 <- "noninf"
type_out <- tolower(cancerType)
cond1_out <- "infiltrated"
cond2_out <- "noninf"

#' Run ssGSEA to calculate an enrichment score for each sample
#'
#' @param expression gene expression. genes on rows, samples on columns
#' @param signatures list of signatures. Each signature is a vector of genes. No direction is supported
#' @return A matrix with the results of the analysis
#'
ssGSEA <- function(expression, signatures) {
  require(GSVA)
  gsva_analysis <- gsva(expression, signatures, mx.diff=FALSE, verbose=TRUE, method='ssgsea', 
                        tau=0.25, parallel.sz=0, rnaseq = TRUE, ssgsea.norm=TRUE)
  return(gsva_analysis)
}

infiltSignatures <- eval(parse(text = load("signatures.RData")))
names(infiltSignatures)[1:5]


Tcell_sign_genes <-  c(union(infiltSignatures[["Tcell_infiltration_1"]][["gene"]],
            infiltSignatures[["Tcell_infiltration_2"]][["gene"]]))



# all_cond1_ID_file <- paste0(cond1_out, "_ID.Rdata")
# stopifnot(file.exists(all_cond1_ID_file))
# all_cond1_ID <- eval(parse(text = load(all_cond1_ID_file)))
# 
# all_cond2_ID_file <- paste0(cond2_out, "_ID.Rdata")
# stopifnot(file.exists(all_cond2_ID_file))
# all_cond2_ID <- eval(parse(text = load(all_cond2_ID_file)))
# 
cmp_name <- paste0("TCGA", type_out, "_", cond1_out, "_", cond2_out)
# 

outFold="."

logFile <- file.path(outFold, paste0(cmp_name, "_prep_data_logFile.txt"))
system(paste0("rm -f ", logFile))


# /mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain/[...].rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt
exprFolder <- file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_gdac/TCGA_transcriptome_Oct2016/plain")

txt <- paste0("... using exprFolder\t=\t", exprFolder, "\n")
printAndLog(txt, logFile)

# outputSamp1File <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name,  paste0(cond1_out, "_ID.Rdata"))
# outputSamp2File <- file.path(setDir, "/mnt/ed4/marie/other_datasets", cmp_name, paste0(cond2_out, "_ID.Rdata"))
# cat(outputExprFile, "\n")
# stopifnot(!file.exists(outputExprFile))
# stopifnot(!file.exists(outputSamp1File))
# stopifnot(!file.exists(outputSamp2File))
# 
# system(paste0("mkdir -p ", dirname(outputSamp1File)))
# system(paste0("mkdir -p ", dirname(outputSamp2File)))

all_exprFiles <- list.files(exprFolder, full.names = TRUE, pattern=paste0(".+illuminahiseq.+rnaseqv2.+RSEM.+"))

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
geneIDs <- sapply(geneIDs, function(x) x[[1]])

stopifnot(all(Tcell_sign_genes %in% geneIDs))


signatureRows <- unlist(sapply( Tcell_sign_genes, function(x){
  grep(paste0("^", x, "\\|.+"), exprDT$`Hybridization-REF`) }))

stopifnot(length(signatureRows) == length(Tcell_sign_genes))


Tcell_sign_genes_ID <- exprDT$`Hybridization-REF`[signatureRows]

Tcell_sign_genes_ID

Tcell_sign_set <- list(Tcell_sign = Tcell_sign_genes_ID)


rownames(exprDT) <- exprDT[,"Hybridization-REF"]

exprDT[,"Hybridization-REF"] <- NULL

exprDT[1:3,1:3]

# expr_samples <- colnames(exprDT)
# expr_samples <- substr(x = expr_samples, start = 1, stop = 15)
colnames(exprDT) <- substr(x = colnames(exprDT), start = 1, stop = 15)
expr_samples <- colnames(exprDT)

exprDT[1:5,1:5]

exprMat <- data.matrix(exprDT)
exprMat[1:5,1:5]
stopifnot(is.numeric(exprMat[1,1]))

SKCM_signScores <- ssGSEA(expression = exprMat , 
                     signatures = Tcell_sign_set) 

head(SKCM_signScores)

# save(SKCM_signScores, file = "SKCM_signScores.Rdata")

### scores: 0-25% = low, 25-50 = mid, 75-100 = high infiltration

SKCM_signScores <- setNames(SKCM_signScores[1,], colnames(SKCM_signScores))

# <= 25
# >= 75
thresh25 <- quantile(SKCM_signScores, probs = 0.25)
thresh75 <- quantile(SKCM_signScores, probs = 0.75)

lowInf_ID <- names(SKCM_signScores[SKCM_signScores <=  thresh25])
length(lowInf_ID)

highInf_ID <- names(SKCM_signScores[SKCM_signScores >=  thresh75])
length(highInf_ID)

head(highInf_ID)
head(lowInf_ID)

save(lowInf_ID, file = "lowInf_ID.Rdata")
save(highInf_ID, file = "highInf_ID.Rdata")

