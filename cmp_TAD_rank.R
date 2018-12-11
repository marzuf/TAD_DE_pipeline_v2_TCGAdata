startTime <- Sys.time()

cat(paste0("> Rscript cmp_TAD_rank.R\n"))

# Rscript cmp_TAD_rank.R OUTPUT_FOLDER/TCGAstad_msi_gs OUTPUT_FOLDER_TCGA_RUN3_WRONG_FPKM/TCGAstad_msi_gs
# Rscript cmp_TAD_rank.R OUTPUT_FOLDER/TCGAcoad_msi_mss OUTPUT_FOLDER_TCGA_RUN1/TCGAcrc_msi_mss 
# Rscript cmp_TAD_rank.R OUTPUT_FOLDER/TCGAbrca_lum_bas OUTPUT_FOLDER_TCGA_RUN1/TCGAbrca_lum_bas 
setDir <- ""

source("analysis_utils.R")
source(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18/analysis_utils.R"))

outFold <- "CMP_TAD_RANK"
system(paste0("mkdir -p ", outFold))

infold0 <- "OUTPUT_FOLDER/TCGAstad_msi_gs"
infold1 <- "OUTPUT_FOLDER_TCGA_RUN3_WRONG_FPKM/TCGAstad_msi_gs"

plotType <- "svg"
myHeight <- 7
myWidth <- 7

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
infold0 <- args[1]
infold1 <- args[2]

script11_name <- "11_runEmpPvalCombined"

outFolder0 <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", infold0)
outFolder1 <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", infold1)

ds0 <- basename(infold0)
ds1 <- basename(infold1)

v0_pval_file <- file.path(outFolder0, script11_name, "emp_pval_combined.Rdata")
stopifnot(file.exists(v0_pval_file))
v0_pvals <- eval(parse(text = load(v0_pval_file)))
v0_pvals <- p.adjust(v0_pvals, method="BH")

v1_pval_file <- file.path(outFolder1, script11_name, "emp_pval_combined.Rdata")
stopifnot(file.exists(v1_pval_file))
v1_pvals <- eval(parse(text = load(v1_pval_file)))
v1_pvals <- p.adjust(v1_pvals, method="BH")

intersect_tads <- intersect(names(v0_pvals), names(v1_pvals))

mySub <- paste0("(# TADs = ", length(intersect_tads), ")")

myx <- v0_pvals[intersect_tads]
myy <- v1_pvals[intersect_tads]

outFile <- file.path(outFold, paste0(ds1, "_vs_", ds0, "_cmp_TADpvals.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot(
#   xlab = paste0(ds0, " - TAD adj. pval"),
#   ylab = paste0(ds1, " - TAD adj. pval"),
#   main = paste0(ds1, " vs. ", ds0, "results"),
#   x = v0_pvals[intersect_tads],
#   y = v1_pvals[intersect_tads],
#   pch =16, cex =0.7
# )
densplot(
  xlab = paste0(ds0, " - TAD adj. pval"),
  ylab = paste0(ds1, " - TAD adj. pval"),
  main = paste0(ds1, " vs. ", ds0, "results"),
  x = myx,
  y = myy,
         pch = 16, cex = 0.7
)
mtext(side=3, text = mySub)
add_curv_fit(x = myx,
             y= myy, withR2 = FALSE, lty=2, col="darkgray")

addCorr(x=myx,
        y = myy,
        bty="n",
        legPos="bottomright")

foo <- dev.off()
cat(paste0("written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


