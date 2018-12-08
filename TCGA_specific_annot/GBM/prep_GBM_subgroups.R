gbm_annotDT <- read.delim("GBM_subtype_classif.csv", stringsAsFactors = FALSE, sep=",")

gbm_annotDT <- gbm_annotDT[!is.na(gbm_annotDT$Cluster),]

table(gbm_annotDT$Cluster)
# Classical      G-CIMP Mesenchymal      Neural   Proneural 
# 147          39         160          96         102 

classical_ID <- gbm_annotDT$sample.id[gbm_annotDT$Cluster == "Classical"]
length(classical_ID)
neural_ID <- gbm_annotDT$sample.id[gbm_annotDT$Cluster == "Neural"]
length(neural_ID)
proneural_ID <- gbm_annotDT$sample.id[gbm_annotDT$Cluster == "Proneural"]
length(proneural_ID)
mesenchymal_ID <- gbm_annotDT$sample.id[gbm_annotDT$Cluster == "Mesenchymal"]
length(mesenchymal_ID)

save(classical_ID, file ="classical_ID.Rdata")
save(neural_ID, file ="neural_ID.Rdata")
save(proneural_ID, file ="proneural_ID.Rdata")
save(mesenchymal_ID, file ="mesenchymal_ID.Rdata")