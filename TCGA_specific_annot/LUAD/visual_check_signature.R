curr_ds <-"TCGAluad_nonsmoker_smoker"

SSHFS <- TRUE


outFold <- "VISUAL_CHECK"
system(paste0("mkdir -p ", outFold))

setDir <- ifelse(SSHFS, "/media/electron", "")

inputFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")

source(file.path(inputFolder, paste0("run_settings_", curr_ds,".R")))

samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))

rnaDT <- eval(parse(text = load(file.path(setDir, rnaseqDT_file))))


rnaDT[1:5,1:5]

dim(rnaDT)

# remove the columns with only 0 (otherwise cannot scale)
tmpDT <- t(rnaDT)
tmpDT <- tmpDT[,-which(colSums(tmpDT)==0)]
rna_pca <- prcomp(tmpDT, scale=TRUE, center=TRUE)

plot(rna_pca)

plot_pca <- function(prcomp_obj, PC_x = 1, PC_y = 2, ...) {
  var_exp <- prcomp_obj$sdev^2/sum(prcomp_obj$sdev^2)
  PC_x <- as.numeric(as.character(PC_x))
  PC_y <- as.numeric(as.character(PC_y))
  stopifnot(!is.na(PC_x))
  stopifnot(!is.na(PC_y))
  xlab <- paste0("PC", PC_x, " (", round(var_exp[PC_x]*100,2), "%)")
  ylab <- paste0("PC", PC_y, " (", round(var_exp[PC_y]*100,2), "%)")
  stopifnot(length(var_exp)  >= PC_x)
  stopifnot(length(var_exp)  >= PC_y)
  xpos <- prcomp_obj$x[,PC_x]
  ypos <- prcomp_obj$x[,PC_y]
  stopifnot(length(xpos) == length(ypos))
  
  plot(x=xpos,
       y=ypos,
       xlab=xlab,
       ylab=ylab,
       ...)
}

mycols <- ifelse(rownames(rna_pca$x) %in% samp1, "red", 
                 ifelse(rownames(rna_pca$x) %in% samp2, "blue", 
                 NA))
stopifnot(!is.na(mycols))

outFile <- file.path(outFold, "pca_smoker_smoker.png" )
png(outFile, height=400, width=400)
plot_pca(rna_pca, pch=16, cex=0.7, col=mycols, main=paste0(curr_ds))
legend("topleft", legend = c("nonsmoker", "smoker"), pch=16, col = c("red", "blue"), bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



