############################################################
## COMPLETE STAT3 ChIP-seq PROMOTER ANALYSIS
## INPUT: BED FILE ONLY
############################################################

############################################################
## 0. Load libraries
############################################################
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(ChIPpeakAnno)
  library(rtracklayer)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(UpSetR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)  # for mapIds etc.
  library(ggplot2) 
  library(ChIPseeker)# for plotting (dotplot uses ggplot2)
})
############################################################
## 1. Set paths
############################################################
setwd("PATH")
list.files()
peak_file <- "STAT3_chipseq_peaks.narrowPeak"
############################################################
## 2. Import BED/NARROW_PEAKS peaks
############################################################
peaks_gr1 <- toGRanges(peak_file, format = "narrowPeak", header = FALSE)
cat("Total peaks:", length(peaks_gr1), "\n")
############################################################
## 3. Gene annotation (EnsDb v86)
############################################################
ensdb1 <- EnsDb.Hsapiens.v86
ensembl_genes1 <- toGRanges(ensdb1, feature = "gene")
ensembl_genes1 <- keepStandardChromosomes(ensembl_genes1, pruning.mode = "coarse")
suppressWarnings(seqlevelsStyle(ensembl_genes1) <- "UCSC")
############################################################
## 4. Annotate peaks to nearest TSS
############################################################
peak_anno1 <- annotatePeakInBatch(
  peaks_gr1,
  AnnotationData = ensembl_genes1,
  output = "nearestLocation",
  FeatureLocForDistance = "TSS"
)
peak_anno_df1 <- as.data.frame(peak_anno1)

# Map gene symbols (may return list-column)
peak_anno_df1$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = peak_anno_df1$feature,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
# 'select()' returned 1:many mapping warning is normal if some keys map multiple symbols

# --- FIX: Flatten gene_symbol list-column to character vector ---
peak_anno_df1$gene_symbol <- sapply(peak_anno_df1$gene_symbol, function(x) {
  if (length(x) == 0) return(NA_character_)  # handle empty lists
  else return(as.character(x[1]))
})

# Write to CSV now that gene_symbol column is flattened
write.csv(peak_anno_df1, "STAT3_peak_annotation_from_NP.csv", row.names = FALSE)
############################################################
## 5. Promoter definition (-2 kb / +1 kb)
############################################################
txdb1 <- TxDb.Hsapiens.UCSC.hg38.knownGene

genes_gr1 <- genes(txdb1, single.strand.genes.only = FALSE)
genes_gr1 <- unlist(genes_gr1)
genes_gr1 <- keepStandardChromosomes(genes_gr1, pruning.mode = "coarse")
suppressWarnings(seqlevelsStyle(genes_gr1) <- "UCSC")

promoters_gr1 <- promoters(
  genes_gr1,
  upstream = 2000,
  downstream = 1000
)

############################################################
## 6. Overlap peaks with promoters
############################################################
ov_prom1 <- findOverlaps(peaks_gr1, promoters_gr1)

promoter_peaks1 <- peaks_gr1[queryHits(ov_prom1)]
promoter_genes1 <- promoters_gr1[subjectHits(ov_prom1)]

############################################################
## 7. Distance to TSS (strand-aware)
############################################################
gene_strand1 <- as.character(strand(promoter_genes1))
gene_strand1[gene_strand1 == "*"] <- "+"

gene_TSS1 <- ifelse(
  gene_strand1 == "+",
  start(promoter_genes1),
  end(promoter_genes1)
)

peak_center1 <- start(promoter_peaks1) + width(promoter_peaks1) %/% 2
distance_to_TSS1 <- peak_center1 - gene_TSS1

############################################################
## 8. Promoter dataframe (PEAK + PROMOTER + GENE ID)
############################################################
promoter_df1 <- data.frame(
  peak_id         = names(promoter_peaks1),
  peak_chr        = as.character(seqnames(promoter_peaks1)),
  peak_start      = start(promoter_peaks1),
  peak_end        = end(promoter_peaks1),
  
  promoter_chr    = as.character(seqnames(promoter_genes1)),
  promoter_start  = start(promoter_genes1),
  promoter_end    = end(promoter_genes1),
  
  gene_id         = names(promoter_genes1),  # ENTREZID
  gene_strand     = gene_strand1,
  gene_TSS        = gene_TSS1,
  distance_to_TSS = distance_to_TSS1,
  
  stringsAsFactors = FALSE
)

############################################################
## 9. TRUE gene body coordinates (TxDb)
############################################################
############################################################
## 9. TRUE gene body coordinates (TxDb) — TxDb-SAFE
############################################################

gene_body_gr <- genes(txdb1, single.strand.genes.only = FALSE)
gene_body_gr <- unlist(gene_body_gr)

## Extract ENTREZ IDs from names()
gene_body_df1 <- data.frame(
  gene_id = names(gene_body_gr),
  gene_body_chr   = as.character(seqnames(gene_body_gr)),
  gene_body_start = start(gene_body_gr),
  gene_body_end   = end(gene_body_gr),
  stringsAsFactors = FALSE
) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)


promoter_df1 <- promoter_df1 %>%
  dplyr::left_join(
    gene_body_df1,
    by = "gene_id"
  )

############################################################
## 10. Gene symbol + Ensembl ID mapping
############################################################
promoter_df1$symbol <- mapIds(
  org.Hs.eg.db,
  keys = promoter_df1$gene_id,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

promoter_df1$ensembl_gene_id <- mapIds(
  org.Hs.eg.db,
  keys = promoter_df1$gene_id,
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)

############################################################
## 11. Protein-coding gene filter (EnsDb)
############################################################
gene_info1 <- genes(
  ensdb1,
  columns = c("gene_id", "gene_biotype")
) %>%
  as.data.frame() %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

promoter_df_pc1 <- promoter_df1 %>%
  dplyr::left_join(
    gene_info1,
    by = c("ensembl_gene_id" = "gene_id")
  ) %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::filter(!grepl("^LOC", symbol))

############################################################
## 12. Save results
############################################################
write.csv(
  promoter_df1,
  "STAT3_promoter_peaks_ALL_from1_NP.csv",
  row.names = FALSE
)

write.csv(
  promoter_df_pc1,
  "STAT3_promoter_peaks_PROTEIN_CODING1_from_NP.csv",
  row.names = FALSE
)

############################################################
## 13. Genomic element UpSet plot
############################################################
res1 <- genomicElementUpSetR(
  peaks_gr1,
  TxDb.Hsapiens.UCSC.hg38.knownGene
)

upset(
  res1$plotData,
  nsets = ncol(res1$plotData),
  nintersects = NA
)

############################################################
## 14. Peak distribution relative to TSS
############################################################
if (is.null(mcols(peaks_gr1)$score)) {
  mcols(peaks_gr1)$score <- 1
}

binOverFeature(
  peaks_gr1,
  nbins = 20,
  annotationData = ensembl_genes1,
  xlab = "Peak distance from TSS (bp)",
  ylab = "Peak count",
  main = "Distribution of STAT3 peaks around TSS"
)

############################################################
## 15. Genomic feature distribution
############################################################
genomicElementDistribution(
  peaks_gr1,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
)

macs_peaks1 <- GRangesList(STAT3 = peaks_gr1)

genomicElementDistribution(
  macs_peaks1,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
)

############################################################
## 16. Final summary
############################################################
cat("\n===== FINAL SUMMARY =====\n")
cat("Total peaks: ", length(peaks_gr1), "\n")
cat("Promoter-overlapping peaks: ", nrow(promoter_df1), "\n")
cat("Unique promoter genes (all): ",
    length(unique(promoter_df1$symbol)), "\n")
cat("Protein-coding promoter genes: ",
    length(unique(promoter_df_pc1$symbol)), "\n")

###########################################################
#### GSEA Analysis
###########################################################
peaks <- readPeakFile("STAT3_chipseq_peaks.narrowPeak")

peakAnno <- annotatePeak(
  peaks,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

peakAnno1 <- annotatePeak(
  peaks,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  tssRegion = c(-2000, 1000),
  annoDb = "org.Hs.eg.db"
)

gene_df <- as.data.frame(peakAnno1)

gene_score <- gene_df %>%
  filter(!is.na(ENSEMBL)) %>%
  group_by(ENSEMBL) %>%
  summarise(score = max(V7, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(score))

gene_ids <- bitr(
  gene_score$ENSEMBL,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)


gene_score2 <- inner_join(
  gene_score,
  gene_ids,
  by = c("ENSEMBL" = "ENSEMBL")
)

gene_score2 <- gene_score2 %>%
  group_by(ENTREZID) %>%
  summarise(score = max(score)) %>%
  ungroup()

geneList <- gene_score2$score
names(geneList) <- gene_score2$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

any(duplicated(names(geneList)))   # MUST be FALSE
length(geneList)                   # ~300–800 is typical

entrez_genes <- unique(gene_score2$ENTREZID)

ego_bp <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.4,  # more lenient
  qvalueCutoff = 0.4
)

ego_cc <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.4,  # more lenient
  qvalueCutoff = 0.4
)


ego_mf <- enrichGO(
  gene = entrez_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.4,  # more lenient
  qvalueCutoff = 0.4
)

sum(ego_bp@result$p.adjust < 0.4)  # see how many pass now

dotplot(ego_bp, showCategory = 10)


