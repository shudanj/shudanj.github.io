
# Project 2: Gene Ontology Analysis

#using the same EdgeR pipeline to get top differentially expressed genes
# as the previous project

library(readr)
library(edgeR)
library(tidyr)
library(dplyr)
library(clusterProfiler)
library(org.Dr.eg.db)   #using zebrafish database, db must match organism gene ids
library(enrichplot)

slc_tissue <- read.csv("zebrafish_slc_tissue.csv", row.names = 1) #load in csv file

wt_slc <- grep("^WT", colnames(slc_tissue), value = TRUE)

py_slc <- grep("^PY", colnames(slc_tissue), value = TRUE)  #Grep pulls every column beginning with WT (control) or PY (treatment)

#Telling edgeR the counts and reps for each sample

counts <- slc_tissue[, c(wt_slc, py_slc)]
group  <- factor(c(rep("WT", 4), rep("PY", 5)))

#Normalization, filters out counts that are too low 
dge  <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)

#Statistics 
design <- model.matrix(~group)  #sets WT as the control/baseline
dge    <- estimateDisp(dge, design)  #accounts for RNA-seq overdispersion
fit    <- glmQLFit(dge, design) #applies a model
qlf    <- glmQLFTest(fit, coef = 2)   #runs the hypothesis test (significance test)

results      <- topTags(qlf, n = Inf)$table
results$gene <- rownames(results)

#Assigning significance/P-values
results$significance <- "Not Significant"
results$significance[results$FDR < 0.05 & results$logFC >  1] <- "Up in PY"
results$significance[results$FDR < 0.05 & results$logFC < -1] <- "Down in PY"

top_genes <- results[results$FDR < 0.05 & abs(results$logFC) > 3, ]

print(top_genes) #checkpoint to make sure the pipeline is working 

sig_ids <- bitr(
  rownames(top_genes),         # your significant gene symbols
  fromType = "SYMBOL",          
  toType   = "ENTREZID",
  OrgDb    = org.Dr.eg.db       #zebrafish
)

#Convert the full background (all tested genes)
bg_ids <- bitr(
  rownames(results),            #all genes that passed filtering
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Dr.eg.db
)

#Run GO enrichment
go_res <- enrichGO(
  gene          = sig_ids$ENTREZID,
  universe      = bg_ids$ENTREZID,
  OrgDb         = org.Dr.eg.db,
  ont           = "BP",           #bp = biological process, there are other functions that can be substituted for
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE            #converts Entrez IDs back to gene symbols in output
)

#checkpoint
head(as.data.frame(go_res)[, c("Description", "GeneRatio", "p.adjust", "geneID")])

#Plot
dotplot(go_res, showCategory = 15, title = "GO Biological Process in Treatment Group")

