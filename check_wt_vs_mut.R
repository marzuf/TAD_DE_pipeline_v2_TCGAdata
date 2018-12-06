all_datasets <- c("blca", "lusc", "stad", "kich", "sarc", "sarc", "tgct", "lgg", "cesc", "hnsc", "sarc")
all_cond1 <- c("norm", "norm", "norm", "norm", "lms", "ddlps", "sem", "IDHwt", "adeno", "HPVneg", "ddlps")
all_cond2 <- c("blca", "lusc", "gs", "kich", "mfs", "mfs", "nonsem", "IDHmutnc", "squam", "HPVpos", "lms")


all_datasets<-c()
all_cond1<-c()
all_cond2<-c()
all_datasets <- c(all_datasets, "stad", "ucec", "thca", "stad", "skcm", "skcm", "paad", "luad", "lihc", "laml", "acc", "luad")
all_cond1 <- c(all_cond1, "EBVpos", "msi", "wt", "msi", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "mutKRAS")
all_cond2 <- c(all_cond2, "gs", "cnl", "mutBRAF", "gs", "mutCTNNB1", "mutBRAF", "mutKRAS", "mutKRAS", "mutCTNNB1", "mutFLT3", "mutCTNNB1", "mutEGFR")
# => ok



checkFile <- "DATA_SUMMARY/tcga_datasets_DT.Rdata"
checkDT <- eval(parse(text = load(checkFile)))
#type subtype      class nSamplesAnnot nSamplesAnnotExpr nSamplesAnnotExprMut
#1  ACC     ACC        ACC            76                76                   76

gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v5.RData")
x <- load(gamFile)
stopifnot("pancan_gam" %in% x)
gamDT <- pancan_gam
gamDT[1:3,1:3]


for(i in seq_along(all_datasets)){
  
  cancerType <- all_datasets[i]
  cond1 <- all_cond1[i]
  cond2 <- all_cond2[i]
  
  if( ! grepl("^mut.+$", cond2)) next
  
  dataset <- paste0("TCGA", cancerType, "_", cond1, "_", cond2)
  samp1_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", dataset, paste0(cond1, "_ID.Rdata"))
  samp2_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", dataset, paste0(cond2, "_ID.Rdata"))
  samp1 <- eval(parse(text = load(samp1_file)))
  samp2 <- eval(parse(text = load(samp2_file)))
  
  # checkType <- toupper(cancerType)
  # if(grepl("^mut.+$", cond2)) {
  # } else {
  #   if(cond1 != "norm" & cond1 %in% checkDT$subtype[checkDT$type == checkType]) {
  #     nCond1 <- checkDT$nSamplesAnnotExpr[toupper(checkDT$subtype) == toupper(cond1) & [checkDT$type == checkType]]
  #     stopifnot(nCond1 == length(samp1))
  #   }
  # }
  
  stopifnot(samp1 %in% rownames(gamDT))
  stopifnot(samp2 %in% rownames(gamDT))
  
  mut2 <- gsub("^mut", "", cond2)
  mutCol2 <- grep(paste0("MUT\\.", mut2), colnames(gamDT))
  
  stopifnot(all(gamDT[samp2, mutCol2] == 1))
  
  if(grepl("^mut.+$", cond1)) {
    mut1 <- gsub("^mut", "", cond1)
    mutCol1 <- grep(paste0("MUT\\.", mut1), colnames(gamDT))
    stopifnot(all(gamDT[samp1, mutCol1] == 1))
  } else {
    stopifnot(all(gamDT[samp1, mutCol2] == 0))
  }
}