curr_ds <-"TCGAskcm_lowInf_highInf"

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")

outFold <- "VISUAL_CHECK"
system(paste0("mkdir -p ", outFold))

inputFolder <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")

source(file.path(inputFolder, paste0("run_settings_", curr_ds,".R")))

samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))

rnaDT <- eval(parse(text = load(file.path(setDir, rnaseqDT_file))))

rnaDT[1:5,1:5]

plotType <- "png"
myWidth <- 400
myHeight <- 400

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

plot_pca(rna_pca, pch=16, cex=0.7, col=mycols, main=paste0(curr_ds))
legend("topleft", legend = c("lowInf", "highInf"), pch=16, col = c("red", "blue"), bty="n")


entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)

load("signatures.RData")

signature_genes <- c(signatures[["Tcell_infiltration_1"]][["gene"]], signatures[["Tcell_infiltration_2"]][["gene"]])
# retrieve entrez ID
signature_entrez <- sapply(signature_genes, function(x)as.character(entrezDT$entrezID[entrezDT$symbol == x]))

DE_DT <- eval(parse(text = load(file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER", 
                                          curr_ds, "1_runGeneDE", "DE_topTable.Rdata"))))

DE_DT <- DE_DT[order(DE_DT$adj.P.Val, decreasing=F),]

outFile <- file.path(outFold, paste0("pca_infiltration_pval_all.", plotType))
do.call(plotType, list(outFile, height=myHeight, width = myWidth))
plot(-log10(DE_DT$adj.P.Val), type="h",
     main = paste0(curr_ds, " - infiltration signature - pval"),
      xlab="",     
     ylab="-log10 adj. pval",
     col=ifelse(DE_DT$genes %in% signature_entrez, "red", "gray"))
legend("topright", 
       legend=c(paste0("n=", length(DE_DT$genes) ),
                # paste0(sum(signature_entrez %in% DE_DT$genes), "/", length(DE_DT$genes) ),
                paste0(sum(signature_entrez %in% DE_DT$genes), "/", length(signature_genes) )),
       bty="n"
       )
mtext(side=3, text = paste0("(all genes)"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("pca_infiltration_pval_top1000.", plotType))
do.call(plotType, list(outFile, height=myHeight, width = myWidth))
plot(-log10(DE_DT$adj.P.Val[1:1000]), type="h",
     main = paste0(curr_ds, " - infiltration signature - pval"),
     ylab="-log10 adj. pval",
     col=ifelse(DE_DT$genes[1:1000] %in% signature_entrez, "red", "gray"))
legend("topright", 
       legend=c(paste0("n=", length(DE_DT$genes[1:1000]) ),
                # paste0(sum(signature_entrez %in% DE_DT$genes[1:1000]), "/", length(DE_DT$genes[1:1000]) ),
                paste0(sum(signature_entrez %in% DE_DT$genes[1:1000]), "/", length(signature_genes) )),
       bty="n"
)
mtext(side=3, text = paste0("(top 1000 genes)"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

DE_DT <- DE_DT[order(DE_DT$logFC, decreasing = T),]

outFile <- file.path(outFold, paste0("pca_infiltration_logfc_all.", plotType))
do.call(plotType, list(outFile, height=myHeight, width = myWidth))
plot(DE_DT$logFC, type="h",
     main = paste0(curr_ds, " - infiltration signature - logFC"),
     xlab="",     
     ylab="log FC",
     col=ifelse(DE_DT$genes %in% signature_entrez, "red", "gray"))
legend("topright", 
       legend=c(paste0("n=", length(DE_DT$genes) ),
                # paste0(sum(signature_entrez %in% DE_DT$genes), "/", length(DE_DT$genes) ),
                paste0(sum(signature_entrez %in% DE_DT$genes), "/", length(signature_genes) )),
       bty="n"
)
mtext(side=3, text = paste0("(all genes)"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("pca_infiltration_logfc_top1000.", plotType))
do.call(plotType, list(outFile, height=myHeight, width = myWidth))
plot(DE_DT$logFC[1:1000], type="h",
     main = paste0(curr_ds, " - infiltration signature - logFC"), 
     ylab="log FC",
     col=ifelse(DE_DT$genes[1:1000] %in% signature_entrez, "red", "gray"))
legend("topright", 
       legend=c(paste0("n=", length(DE_DT$genes[1:1000]) ),
                # paste0(sum(signature_entrez %in% DE_DT$genes[1:1000]), "/", length(DE_DT$genes[1:1000]) ),
                paste0(sum(signature_entrez %in% DE_DT$genes[1:1000]), "/", length(signature_genes) )),
       bty="n"
)
mtext(side=3, text = paste0("(top 1000 genes)"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

