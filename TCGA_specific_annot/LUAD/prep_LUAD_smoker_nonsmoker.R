suppTable <- "LUAD_tableS7.csv"
suppDT <- read.delim(suppTable, sep=",", stringsAsFactors = FALSE)
head(suppDT)
table(suppDT$Smoking.Status)
# Current reformed smoker for > 15 years Current reformed smoker for < or = 15 years 
# 69                                          73 
# Current smoker                         Lifelong Non-smoker 
# 45                                          32 
# [Not Available] 
# 11 

suppDT <- suppDT[!is.na(suppDT$Smoking.Status),]

suppDT_smoker_ID <- suppDT$Tumor.ID[grepl("smoker", suppDT$Smoking.Status) &
                                      !grepl("Non-smoker", suppDT$Smoking.Status)]
length(suppDT_smoker_ID)
# 187

suppDT_nonsmoker_ID <- suppDT$Tumor.ID[grepl("Non-smoker", suppDT$Smoking.Status)]
length(suppDT_nonsmoker_ID)
# 32

stopifnot(length(intersect(suppDT_smoker_ID, suppDT_nonsmoker_ID)) == 0)

stopifnot(length(suppDT_nonsmoker_ID) + length(suppDT_smoker_ID) + 11 == nrow(suppDT))


tcgaTable <- "clinical.project-TCGA-LUAD.2018-12-07/LUAD_exposure.tsv"
tcgaDT <- read.delim(tcgaTable, sep="\t", stringsAsFactors = F)

tcgaDT <- tcgaDT[!is.na(tcgaDT$years_smoked),]


tcgaDT_smoker_ID <- tcgaDT$submitter_id[tcgaDT$years_smoked != "--"]
length(tcgaDT_smoker_ID)
# 194

tcgaDT_nonsmoker_ID <- tcgaDT$submitter_id[tcgaDT$years_smoked == "--"]
length(tcgaDT_nonsmoker_ID)
# 391


sum(suppDT_smoker_ID %in% tcgaDT_smoker_ID)
# 88
sum(suppDT_smoker_ID %in% tcgaDT_nonsmoker_ID)  # => some of the tcgaDT nonsmoker should indeed be NA data ?!
# 99


sum(suppDT_nonsmoker_ID %in% tcgaDT_nonsmoker_ID)
# 32
sum(suppDT_nonsmoker_ID %in% tcgaDT_smoker_ID)
# 0

# don't use the tcgaDT_nonsmoker -> because I don't know if not smoker or NA data !

smoker_ID <- union(suppDT_smoker_ID, tcgaDT_smoker_ID)
stopifnot(!smoker_ID %in% suppDT_nonsmoker_ID)

nonsmoker_ID <- suppDT_nonsmoker_ID
stopifnot(!nonsmoker_ID %in% suppDT_smoker_ID)
stopifnot(!nonsmoker_ID %in% tcgaDT_smoker_ID)

length(smoker_ID)
length(nonsmoker_ID)

save(smoker_ID, file = "smoker_ID.Rdata")
save(nonsmoker_ID, file = "nonsmoker_ID.Rdata")