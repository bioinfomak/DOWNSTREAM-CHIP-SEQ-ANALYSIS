# ============================================================
#      🧬 COMPLETE DOWNSTREAM CHIP-SEQ ANALYSIS IN R
# ============================================================
# Covers:
#   1️⃣ Load peak and signal files
#   2️⃣ Peak annotation
#   3️⃣ Gene extraction
#   4️⃣ GO enrichment
#   5️⃣ Average ChIP-seq signal profile (matrix.gz)
#   6️⃣ Heatmap visualization
#   7️⃣ Output saving
# ============================================================

# -------------------------------
# Load required libraries
# -------------------------------
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(EnsDb.Hsapiens.v86)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(zoo)
})

# -------------------------------
# 1️⃣ Define Input Files
# -------------------------------
peak_file   <- "sample_chipseq_peaks.narrowPeak"   # Peak file
bw_file     <- "sample_chipseq.bw"                 # BigWig (optional)
matrix_file <- "matrix.gz"                         # computeMatrix output file

# Output directories
dir.create("plots", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# -------------------------------
# 2️⃣ Peak Annotation
# -------------------------------
cat("📂 Loading peak file...\n")
peaks <- readPeakFile(peak_file)
cat("✅ Loaded", length(peaks), "peaks\n")

cat("🧭 Annotating peaks...\n")
peakAnno <- annotatePeak(peakS, tssRegion=c(-3000, 3000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
write.table(as.data.frame(peakAnno),
            "tables/peak_annotation_new.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("✅ Peak annotation saved to tables/peak_annotation_new.tsv\n")

# -------------------------------
# 3️⃣ Extract Unique Gene IDs
# -------------------------------
genes <- unique(na.omit(as.data.frame(peakAnno)$geneId))
cat("🧬 Number of unique Entrez gene IDs:", length(genes), "\n")

# -------------------------------
# 4️⃣ GO Enrichment Analysis
# -------------------------------
cat("📊 Running GO enrichment...\n")
ego <- enrichGO(
  gene          = genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

write.table(as.data.frame(ego@result),
            "tables/GO_enrichment_BP.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf("plots/GO_enrichment_dotplot.pdf", width = 8, height = 6)
dotplot(ego, showCategory = 15, title = "GO Enrichment: Biological Process")
dev.off()
cat("✅ GO enrichment saved to tables/GO_enrichment_BP.tsv and plots/GO_enrichment_dotplot.pdf\n")

# -------------------------------
# 5️⃣ Load Signal Matrix (computeMatrix output)
# -------------------------------
cat("📥 Loading signal matrix from computeMatrix...\n")

header_lines <- readLines(gzfile(matrix_file), n = 100)
skip_n <- length(grep("^@", header_lines))
cat("Skipping", skip_n, "header lines...\n")

mat_df <- read.table(gzfile(matrix_file), skip = skip_n, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
signal_mat <- as.matrix(apply(mat_df[, 4:ncol(mat_df)], 2, as.numeric))
cat("✅ Matrix loaded with", nrow(signal_mat), "regions and", ncol(signal_mat), "bins\n")

# -------------------------------
# 6️⃣ Parameter Setup
# -------------------------------
flank_bp      <- 2000
bin_size      <- 50
smooth_window <- 5
focus_region  <- 1000

coords <- seq(-flank_bp, flank_bp, length.out = ncol(signal_mat))
avg_signal <- colMeans(signal_mat, na.rm = TRUE)

avg_df <- data.frame(Position = coords, Signal = avg_signal)
avg_df <- avg_df %>% mutate(Smoothed = zoo::rollmean(Signal, k = smooth_window, fill = NA))

filtered_df <- avg_df %>%
  dplyr::filter(Position >= -focus_region & Position <= focus_region)


# -------------------------------
# 7️⃣ Plot Average ChIP-seq Signal
# -------------------------------
cat("🧩 Plotting average ChIP-seq signal profiles...\n")

p1 <- ggplot(avg_df, aes(x = Position, y = Smoothed)) +
  geom_line(color = "#0072B2", size = 1.2) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Average ChIP-seq Signal (±", flank_bp, " bp)"),
    subtitle = paste0("Bin size: ", bin_size, " bp | Smoothing: ", smooth_window, " bins"),
    x = "Distance from peak center (bp)",
    y = "Smoothed signal"
  )

p2 <- ggplot(filtered_df, aes(x = Position, y = Smoothed)) +
  geom_line(color = "#009E73", size = 1.2) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Zoomed Signal (±", focus_region, " bp around TSS/peak)"),
    x = "Distance from center (bp)",
    y = "Smoothed signal"
  )

ggsave("plots/avg_profile_full.pdf", p1, width = 7, height = 5)
ggsave("plots/avg_profile_zoom.pdf", p2, width = 7, height = 5)
cat("✅ Signal plots saved in ./plots\n")

# -------------------------------
# 8️⃣ Heatmap Visualization
# -------------------------------
cat("🔥 Generating heatmap of top signal regions...\n")

top_regions <- order(rowSums(signal_mat, na.rm = TRUE), decreasing = TRUE)[1:1000]
signal_sorted <- signal_mat[top_regions, ]

pheatmap(
  signal_sorted,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  main = paste("ChIP-seq Signal Heatmap (Top 1000 ±", flank_bp, "bp)")
)

cat("✅ Heatmap generated.\n")

# -------------------------------
# 9️⃣ Save Processed Data
# -------------------------------
write.table(avg_df, "tables/average_signal_full.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_df, "tables/average_signal_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ All analyses complete!\n")
cat("📁 Results saved in ./tables and ./plots\n")
