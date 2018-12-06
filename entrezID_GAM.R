setDir <- "/media/electron"
setDir <- ""

outFold <- "GAM_IDs_ENTREZ"
system(paste0("mkdir -p ", outFold))

# firstly, run using v4 (_GAMv4 folders)
#gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v4.RData")
# then, update: e-mail from Marco 29.11.2018
gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v5.RData")
x=load(gamFile)
x

all_genes_symbol <- unique(unlist(alterations[["alteration_gene"]]))
length(all_genes_symbol)

gffFile <- file.path(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gffDT <- read.delim(gffFile, header=T, stringsAsFactors = FALSE)

sum(all_genes_symbol %in% gffDT$symbol)

all_genes_symbol[!all_genes_symbol %in% gffDT$symbol]
# "IFNL1,IFNL2" "ERG,ETS2" 


alt_to_modif <- c("AMP.pc32_CNA.chr19:39752820-40157464", "DEL.pc32_CNA.chr21:39750686-40196987")

for(alt in alt_to_modif) {
  genes <- alterations[["alteration_gene"]][[alt]]
  alterations[["alteration_gene"]][[alt]] <- unlist(strsplit(x=genes, split = ","))
}

all_genes_symbol <- unique(unlist(alterations[["alteration_gene"]]))
stopifnot(length(all_genes_symbol[!all_genes_symbol %in% gffDT$symbol]) == 0)

# -> ok

entrezID_2_symbol <- setNames(
  gffDT$symbol[gffDT$symbol %in% all_genes_symbol],
  as.character(gffDT$entrezID[gffDT$symbol %in% all_genes_symbol])
)
stopifnot(length(entrezID_2_symbol) == length(all_genes_symbol))
stopifnot(entrezID_2_symbol %in% all_genes_symbol)
stopifnot(all_genes_symbol %in% entrezID_2_symbol)

stopifnot(!duplicated(entrezID_2_symbol))
stopifnot(!duplicated(names(entrezID_2_symbol)))

head(entrezID_2_symbol)
# 57801   388585     6146   390992    54626    54206 
# "HES4"   "HES5"  "RPL22"   "HES3"   "HES2" "ERRFI1" 

alterations_by_entrezID <- lapply(names(entrezID_2_symbol), function(x) {
  geneEntrez <- as.character(x)
  geneSymbol <- as.character(entrezID_2_symbol[geneEntrez])
  stopifnot(!is.na(geneSymbol))
  matchingAlterations <- Filter(function(x) {geneSymbol %in% x}, alterations[["alteration_gene"]])
  names(matchingAlterations)
})
names(alterations_by_entrezID) <- as.character(names(entrezID_2_symbol))

## check for few genes:
geneSymb <- "BRAF"
geneEntrez <- names(entrezID_2_symbol[entrezID_2_symbol == geneSymb])
stopifnot(length(geneEntrez) == 1)
geneAlterations <- alterations_by_entrezID[[paste0(geneEntrez)]]
alterations[["alteration_gene"]][ names(alterations[["alteration_gene"]]) %in% geneAlterations]

outFile <- file.path(outFold, "entrezID_2_symbol.Rdata")
save(entrezID_2_symbol, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "alterations_by_entrezID.Rdata")
save(alterations_by_entrezID, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


