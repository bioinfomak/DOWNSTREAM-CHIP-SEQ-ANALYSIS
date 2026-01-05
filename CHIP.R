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
})

############################################################
## 1. Set paths
############################################################
setwd("C:")

peak_gff_file <- "STAT3_chipseq_peaks.gff3"
bed_file      <- "STAT3_chipseq_summits.bed"

############################################################
## 2. Import peaks
############################################################
peak_gr1 <- import(peak_gff_file)
peak_gr2 <- toGRanges(bed_file, format = "BED", header = FALSE)

############################################################
## 3. Merge overlapping peaks
############################################################
overlaps_peak <- findOverlapsOfPeaks(peak_gr1, peak_gr2)
overlaps_peak <- addMetadata(overlaps_peak, colNames = "score", FUN = mean)

merged_peaks <- overlaps_peak$mergedPeaks
merged_peaks <- keepStandardChromosomes(merged_peaks, pruning.mode = "coarse")
seqlevelsStyle(merged_peaks) <- "UCSC"

cat("Merged peaks:", length(merged_peaks), "\n")

############################################################
## 4. Gene annotation (EnsDb v86)
############################################################
ensdb <- EnsDb.Hsapiens.v86

ensembl_genes <- toGRanges(ensdb, feature = "gene")
ensembl_genes <- keepStandardChromosomes(ensembl_genes, pruning.mode = "coarse")
seqlevelsStyle(ensembl_genes) <- "UCSC"

############################################################
## 5. Annotate peaks to nearest TSS
############################################################
peak_anno <- annotatePeakInBatch(
  merged_peaks,
  AnnotationData = ensembl_genes,
  output = "nearestLocation",
  FeatureLocForDistance = "TSS"
)

peak_anno_df <- as.data.frame(peak_anno)

############################################################
## 6. Map gene symbols (ENSEMBL → SYMBOL)
############################################################
gene_ids <- unique(na.omit(peak_anno_df$feature))

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

peak_anno_df$gene_name <-
  gene_symbols[match(peak_anno_df$feature, names(gene_symbols))]

############################################################
## 7. Final peak annotation table
############################################################
peak_anno_df$peakNames <-
  sapply(peak_anno_df$peakNames, paste, collapse = ";")

peak_anno_df_final <- peak_anno_df %>%
  dplyr::select(
    seqnames, start, end, peakNames,
    feature, gene_name,
    distancetoFeature, shortestDistance,
    insideFeature, fromOverlappingOrNearest
  )

write.csv(peak_anno_df_final,
          "STAT3_peak_annotation_nearest_TSS.csv",
          row.names = FALSE)

############################################################
## 8. Promoter analysis (-2 kb / +1 kb)
############################################################
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

genes_gr <- genes(txdb, single.strand.genes.only = FALSE)
genes_gr <- unlist(genes_gr)
genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
seqlevelsStyle(genes_gr) <- "UCSC"

promoters_gr <- promoters(
  genes_gr,
  upstream = 2000,
  downstream = 1000
)

ov_prom <- findOverlaps(merged_peaks, promoters_gr)

promoter_peaks <- merged_peaks[queryHits(ov_prom)]
promoter_genes <- promoters_gr[subjectHits(ov_prom)]

############################################################
## 9. Distance to TSS
############################################################
gene_strand <- as.character(strand(promoter_genes))
gene_strand[gene_strand == "*"] <- "+"

gene_tss <- ifelse(
  gene_strand == "+",
  start(promoter_genes),
  end(promoter_genes)
)

peak_center <- start(promoter_peaks) + width(promoter_peaks) %/% 2
distance_to_TSS <- peak_center - gene_tss

############################################################
## 10. Promoter dataframe
############################################################
promoter_df <- data.frame(
  chr = as.character(seqnames(promoter_genes)),
  gene_id = names(promoter_genes),   # ENTREZID
  distance_to_TSS = distance_to_TSS,
  stringsAsFactors = FALSE
)

promoter_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = promoter_df$gene_id,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

promoter_df$ensembl_gene_id <- mapIds(
  org.Hs.eg.db,
  keys = promoter_df$gene_id,
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)

############################################################
## 11. Protein-coding filter (EnsDb)
############################################################
gene_info <- genes(
  ensdb,
  columns = c("gene_id", "gene_biotype")
) |> 
  as.data.frame() |> 
  dplyr::distinct(gene_id, .keep_all = TRUE)

promoter_df_pc <- promoter_df %>%
  dplyr::left_join(
    gene_info,
    by = c("ensembl_gene_id" = "gene_id")
  ) %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::filter(!grepl("^LOC", symbol))

############################################################
## 12. Save promoter results
############################################################
write.csv(promoter_df,
          "STAT3_promoter_peaks_ALL.csv",
          row.names = FALSE)

write.csv(promoter_df_pc,
          "STAT3_promoter_peaks_PROTEIN_CODING.csv",
          row.names = FALSE)

############################################################
## 13. Genomic element UpSet plot
############################################################
res <- genomicElementUpSetR(
  merged_peaks,
  TxDb.Hsapiens.UCSC.hg38.knownGene
)

upset(res$plotData,
      nsets = ncol(res$plotData),
      nintersects = NA)

############################################################
## 14. Final summary
############################################################
cat("\n===== FINAL SUMMARY =====\n")
cat("Total merged peaks: ", length(merged_peaks), "\n")
cat("Promoter-overlapping peaks: ", nrow(promoter_df), "\n")
cat("Unique promoter genes (all): ",
    length(unique(promoter_df$symbol)), "\n")
cat("Protein-coding promoter genes: ",
    length(unique(promoter_df_pc$symbol)), "\n")
############################################################
## 15. Peak distribution relative to TSS
############################################################

# binOverFeature requires a score column; set score = 1 if missing
if (is.null(mcols(merged_peaks)$score)) {
  mcols(merged_peaks)$score <- 1
}

binOverFeature(
  merged_peaks,
  nbins = 20,
  annotationData = ensembl_genes,
  xlab = "Peak distance from TSS (bp)",
  ylab = "Peak count",
  main = "Distribution of STAT3 peaks around TSS"
)

############################################################
## 16. Peak distribution over genomic features (single sample)
############################################################

genomicElementDistribution(
  merged_peaks,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
)

############################################################
## 17. Peak distribution over genomic features (UpSet-ready)
############################################################

macs_peaks <- GRangesList(STAT3 = merged_peaks)

genomicElementDistribution(
  macs_peaks,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
)

