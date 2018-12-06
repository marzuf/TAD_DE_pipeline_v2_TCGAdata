setDir <- "/media/electron"

dataset <- "TCGAblca_norm_blca" # => ok
cond1 <- "norm" # =>  ok
cond2 <- "blca" # => ok
# RUN 2 DATASETS:
# Rscript prep_data_clean.R BLCA norm norm BLCA blca
# Rscript prep_data_clean.R LUSC norm norm LUSC lusc
# Rscript prep_data_clean.R STAD norm norm GS gs
# Rscript prep_data_clean.R KICH norm norm KICH kich
# 
# Rscript prep_data_clean.R SARC LMS lms MFS/UPS mfs
# Rscript prep_data_clean.R SARC DDLPS ddlps MFS/UPS mfs
# Rscript prep_data_clean.R TGCT seminoma sem non-seminoma nonsem
# 
# Rscript prep_data_clean.R LGG IDHwt IDHwt IDHmut-non-codel IDHmutnc
# Rscript prep_data_clean.R CESC AdenoCarcinoma adeno SquamousCarcinoma squam
# Rscript prep_data_clean.R HNSC HPV- HPVneg HPV+ HPVpos
# Rscript prep_data_clean.R SARC DDLPS ddlps LMS lms
# 
# Rscript prep_data_clean.R THCA mut.RAS mut.RAS mutBRAF mutBRAF


### FROM RUN 2
all_datasets <- c("blca", "lusc", "stad", "kich", "sarc", "sarc", "tgct", "lgg", "cesc", "hnsc", "sarc")
all_cond1 <- c("norm", "norm", "norm", "norm", "lms", "ddlps", "sem", "IDHwt", "adeno", "HPVneg", "ddlps")
all_cond2 <- c("blca", "lusc", "gs", "kich", "mfs", "mfs", "nonsem", "IDHmutnc", "squam", "HPVpos", "lms")


for(i in seq_along(all_datasets)){
  
  cancerType <- all_datasets[i]
  cond1 <- all_cond1[i]
  cond2 <- all_cond2[i]
  
  dataset <- paste0("TCGA", cancerType, "_", cond1, "_", cond2)
  samp1_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", dataset, paste0(cond1, "_ID.Rdata"))
  samp2_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", dataset, paste0(cond2, "_ID.Rdata"))
  
  samp1v0_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets_TCGA_RUN2_WRONG_PREP_MUT", dataset, paste0(cond1, "_ID.Rdata"))
  samp2v0_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets_TCGA_RUN2_WRONG_PREP_MUT", dataset, paste0(cond2, "_ID.Rdata"))
  
  
  samp1 <- eval(parse(text = load(samp1_file)))
  samp1v0 <- eval(parse(text = load(samp1v0_file)))
  samp2 <- eval(parse(text = load(samp2_file)))
  samp2v0 <- eval(parse(text = load(samp2v0_file)))
  
  stopifnot(setequal(samp1, samp1v0))
  stopifnot(setequal(samp2, samp2v0))
}

# => ok

### FROM RUN 1
# Rscript prep_data_clean.R STAD EBV EBVpos GS gs
# Rscript prep_data_clean.R UCEC MSI msi CN_LOW cnl
# Rscript prep_data_clean.R THCA THCA thca mutBRAF mutBRAF
# Rscript prep_data_clean.R STAD MSI msi GS gs
# Rscript prep_data_clean.R SKCM skcm mutCTNNB1 mutCTNNB1
# Rscript prep_data_clean.R SKCM skcm mutBRAF mutBRAF
# Rscript prep_data_clean.R PAAD paad mutKRAS mutKRAS
# Rscript prep_data_clean.R LUAD luad mutKRAS mutKRAS
# Rscript prep_data_clean.R LIHC lihc mutCTNNB1 mutCTNNB1
# Rscript prep_data_clean.R LAML laml mutFLT3 mutFLT3
# Rscript prep_data_clean.R ACC acc mutCTNNB1 mutCTNNB1

all_datasets <- c("stad", "ucec", "thca", "stad", "skcm", "skcm", "paad", "luad", "lihc", "laml", "acc")
all_cond1 <- c("EBVpos", "msi", "wt", "msi", "wt", "wt", "wt", "wt", "wt", "wt", "wt")
all_cond2 <- c("gs", "cnl", "mutBRAF", "gs", "mutCTNNB1", "mutBRAF", "mutKRAS", "mutKRAS", "mutCTNNB1", "mutFLT3", "mutCTNNB1")

stopifnot(length(all_datasets) == length(all_cond1))
stopifnot(length(all_cond2) == length(all_cond1))
i=11
for(i in seq_along(all_datasets)){
  
  cancerType <- all_datasets[i]
  cond1 <- all_cond1[i]
  cond2 <- all_cond2[i]
  
  dataset <- paste0("TCGA", cancerType, "_", cond1, "_", cond2)
  samp1_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", dataset, paste0(cond1, "_ID.Rdata"))
  samp2_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", dataset, paste0(cond2, "_ID.Rdata"))
  
  if(i == 1) {
    samp1v0_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets_TCGA_RUN1/TCGA_stad_EBVneg_EBVpos/17.08_prepData/EBVneg_ID.Rdata")
    samp2v0_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets_TCGA_RUN1/TCGA_stad_EBVneg_EBVpos/17.08_prepData/EBVpos_ID.Rdata")
    samp1 <- eval(parse(text = load(samp1_file)))
    samp1v0 <- eval(parse(text = load(samp2v0_file)))
    samp2 <- eval(parse(text = load(samp2_file)))
    samp2v0 <- eval(parse(text = load(samp1v0_file)))
    

    
  }
   
  
  
  if(i == 11) {
    
    samp1v0_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets_TCGA_RUN1/TCGA_acc_acc_mutCTNNB1/wt_ID.Rdata")
    samp2v0_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets_TCGA_RUN1/TCGA_acc_acc_mutCTNNB1/mut_ID.Rdata")
    samp1 <- eval(parse(text = load(samp1_file)))
    samp1v0 <- eval(parse(text = load(samp1v0_file)))
    samp2 <- eval(parse(text = load(samp2_file)))
    samp2v0 <- eval(parse(text = load(samp2v0_file)))
    # => 1 sample wrongly annotated in wt  
    

  }
  
  
  if(length(samp1) == length(samp1v0)) {
    stopifnot(setequal(samp1, samp1v0))  
  } else {
    stopifnot(samp1v0 %in% samp1)
  }
  if(length(samp2) == length(samp2v0)) {
    stopifnot(setequal(samp2, samp2v0))
  } else {
    stopifnot(samp1v0 %in% samp1)
  }
  
  
}