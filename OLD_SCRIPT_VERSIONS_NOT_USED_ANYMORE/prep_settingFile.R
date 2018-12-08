startTime <- Sys.time()
cat(paste0("> Rscript prep_settingFile.R\n"))

# Rscript prep_settingFile.R

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/TAD_DE_pipeline_v2_TCGAdata")
source("analysis_utils.R")

outFold <- file.path("PREP_SETTINGFILE")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "prep_settingfile_logFile.txt")
system(paste0("rm -f ", logFile))

settingFileFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")

cancerType <- "CESC"
cond1 <- "AdenoCarcinoma"
cond2 <- "SquamousCarcinoma"
type_out <- tolower(cancerType)
cond1_out <- "adeno"
cond2_out <- "squamous"


# Rscript prep_settingFile.R <cancer type> <cond1> <cond1_out> <cond2> <cond2_out>
# Rscript prep_settingFile.R CESC AdenoCarcinoma adeno SquamousCarcinoma squamous
# Rscript prep_settingFile.R PAAD paad wt mutAMP mut

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
cancerType <- args[1]
cond1 <- args[2]
cond1_out <- args[3]
cond2 <- args[4]
cond2_out <- args[5]
type_out <- tolower(cancerType)



cmp_name <- paste0("TCGA", type_out, "_", cond1_out, "_", cond2_out)

settingFile <- file.path(settingFileFolder, paste0("run_settings_", cmp_name, ".R"))
cat("... settingFile = ", settingFile, "\n")
stopifnot(!file.exists(settingFile))

settingFileTxt <- paste0("

# in this file, settings that are specific for a run on a dataset
             
# gives path to output folder
pipOutFold <- file.path(\"OUTPUT_FOLDER\", \"",cmp_name, "\")
\n
# full path (starting with /mnt/...)
# following format expected for the input
# colnames = samplesID
# rownames = geneID
# !!! geneID are expected not difficulted
\n
# *************************************************************************************************************************\n
# ************************************ SETTINGS FOR 0_prepGeneData\n
# *************************************************************************************************************************\n
\n
rnaseqDT_file <- \"/mnt/ed4/marie/other_datasets/", cmp_name, "/rnaseqDT_v2.Rdata\"
my_sep <- \"\t\"
# input is Rdata or txt file ?
# TRUE if the input is Rdata
inRdata <- TRUE

# can be ensemblID, entrezID, geneSymbol
geneID_format <- \"entrezID\"
stopifnot(geneID_format %in% c(\"ensemblID\", \"entrezID\", \"geneSymbol\"))

# are geneID rownames ? -> \"rn\" or numeric giving the column
geneID_loc <- \"rn\"
stopifnot(geneID_loc == \"rn\" | is.numeric(geneID_loc))

removeDupGeneID <- TRUE

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 1_runGeneDE
# *************************************************************************************************************************

# labels for conditions
cond1 <- \"", cond1_out, "\"
cond2 <- \"", cond2_out, "\"

# path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
sample1_file <- \"/mnt/ed4/marie/other_datasets/", cmp_name, "/", cond1_out, "ID_rnaseq_v2.Rdata\"
sample2_file <- \"/mnt/ed4/marie/other_datasets/", cmp_name, "/", cond2_out, "ID_rnaseq_v2.Rdata\"

minCpmRatio <- 20/888

inputDataType <- \"RSEM\"

")

writeLines(text = settingFileTxt, con = settingFile)
cat(paste0("... written: ", settingFile, "\n"))




######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




