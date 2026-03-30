# CLEANED and DIAGNOSTIC PIPELINE
rm(list = ls())

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(edgeR)
library(biomaRt)
library(tibble)

# 1. Load expression matrix
expression_matrix <- read.table("expression_data.txt", header = TRUE, sep = "\t",
                                row.names = 1, check.names = FALSE)
colnames(expression_matrix) <- sub("^X", "", colnames(expression_matrix))

# 2. CPM normalization (logCPM)
dge <- DGEList(counts = as.matrix(expression_matrix))
dge <- calcNormFactors(dge)
df_cpm <- cpm(dge, log = TRUE, prior.count = 1)  # genes x samples (rows = genes)

# 3. Build meta_data (but DO NOT aggressively parse unless sure)
meta_data <- data.frame(sample_id = colnames(df_cpm), stringsAsFactors = FALSE)
# Diagnostic: show sample IDs
cat("Number of expression samples:", nrow(meta_data), "\n")
print(head(meta_data$sample_id, 20))

# If you *must* select samples whose name contains "_M" or end with "M", do an explicit filter:
# Example: keep samples that contain "_M" (adjust to your naming convention)
meta_data <- meta_data[grepl("_M$|_M_|^M", meta_data$sample_id) | grepl("M_", meta_data$sample_id), , drop = FALSE]
cat("Samples after pattern-based filter:", nrow(meta_data), "\n")
print(head(meta_data$sample_id, 20))

# Exclude any problematic obvious names (only if they exactly match)
exclude_samples <- c("7M_RCS", "19-2M_Pitt")
meta_data <- meta_data[!meta_data$sample_id %in% exclude_samples, , drop = FALSE]

# Subset expression to selected samples
df_cpm_sub <- df_cpm[, meta_data$sample_id, drop = FALSE]

# 4. Load clinical data and inspect
clinical_data <- read.csv("clinical_data.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# DO NOT change Case IDs blindly. First inspect them:
cat("Unique clinical Case examples:\n")
print(head(clinical_data$Case, 20))

# Trim whitespace only:
clinical_data$Case <- trimws(clinical_data$Case)

# 5. Ensure required clinical columns exist
required_cols <- c("Case", "Vital.Status", "OS")
if (!all(required_cols %in% colnames(clinical_data))) {
  stop(paste("Missing required clinical columns:", paste(setdiff(required_cols, colnames(clinical_data)), collapse = ", ")))
}
clinical_subset <- clinical_data[, required_cols, drop = FALSE]

# Convert OS numeric (diagnose non-numeric entries)
clinical_subset$OS_num <- suppressWarnings(as.numeric(as.character(clinical_subset$OS)))
if (any(is.na(clinical_subset$OS_num))) {
  warning("Some OS values could not be coerced to numeric. NAs present. Inspect formatting of OS column.")
}
clinical_subset$OS <- clinical_subset$OS_num
clinical_subset$Vital.Status <- ifelse(tolower(as.character(clinical_subset$Vital.Status)) %in% c("dead","deceased","1","yes","y","true","t"), 1, 0)

# 6. Prepare Ensembl IDs for biomaRt: strip version suffix if present
ensembl_ids <- rownames(df_cpm_sub)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)  # remove .X version if exists

# 7. Get gene symbols from biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids_clean,
  mart = mart
)

gene_annotations <- gene_annotations[gene_annotations$hgnc_symbol != "" & !is.na(gene_annotations$hgnc_symbol), ]
gene_annotations <- gene_annotations[!duplicated(gene_annotations$ensembl_gene_id), ]

# Map back: create mapping vector from original rownames -> hgnc
mapping <- setNames(gene_annotations$hgnc_symbol, gene_annotations$ensembl_gene_id)
# Create a vector of hgnc or NA for each original rowname
hgnc_for_rows <- mapping[ensembl_ids_clean]
# Keep only rows with mapped symbols
keep_idx <- which(!is.na(hgnc_for_rows))
if (length(keep_idx) == 0) stop("No Ensembl IDs matched in biomaRt. Check your Ensembl IDs (versions?).")

expr_annot <- df_cpm_sub[keep_idx, , drop = FALSE]
rownames(expr_annot) <- hgnc_for_rows[keep_idx]

# Collapse duplicate gene symbols by mean
expr_collapsed <- as.data.frame(expr_annot) %>%
  rownames_to_column("gene") %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("gene")
expression_final <- as.matrix(expr_collapsed)

# 7. Signature genes (your list)
signature_genes <- c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "CSPG4", "PTPRZ1", "GPR17")
found_signature <- signature_genes[signature_genes %in% rownames(expression_final)]
cat("Found signature genes:", paste(found_signature, collapse = ", "), "\n")
if (length(found_signature) < length(signature_genes)) {
  warning(paste("Missing signature genes:", paste(setdiff(signature_genes, found_signature), collapse = ", ")))
}
if (length(found_signature) == 0) stop("None of the signature genes were found in expression matrix.")

# 8. Compute signature score: standard approach
# z-score each gene across samples (standardize rows), then sample score = mean z across genes
signature_matrix <- expression_final[found_signature, , drop = FALSE]
signature_matrix_z <- t(scale(t(signature_matrix)))   # genes x samples (standardized per gene)
# If some genes have zero variance, scale() gives NA — handle:
signature_matrix_z[is.na(signature_matrix_z)] <- 0
signature_scores_per_sample <- colMeans(signature_matrix_z, na.rm = TRUE)

# Build signature_df with same Case/sample naming
signature_df <- data.frame(
  Case = names(signature_scores_per_sample),
  score = as.numeric(signature_scores_per_sample),
  stringsAsFactors = FALSE
)

# 9. Diagnostic: check overlap between clinical Case and expression sample names
cat("Number of signature samples:", nrow(signature_df), "\n")
cat("Number of clinical cases:", nrow(clinical_subset), "\n")
common <- intersect(clinical_subset$Case, signature_df$Case)
cat("Number of overlapping Case names:", length(common), "\n")
cat("Examples of overlap (up to 20):\n"); print(head(common, 20))

# If no overlap, show examples to help debug
if (length(common) == 0) {
  cat("No overlap between clinical Case and expression sample names. Inspect naming conventions:\n")
  cat("Sample names (expression):\n"); print(head(signature_df$Case, 30))
  cat("Clinical Case names:\n"); print(head(clinical_subset$Case, 30))
  stop("No overlap — fix Case/sample naming before merging.")
}

# 10. Merge and survival
merged <- merge(clinical_subset, signature_df, by = "Case", all = FALSE)
cat("Merged rows:", nrow(merged), "\n")
print(table(is.na(merged$OS), useNA = "ifany"))
# Remove rows with missing OS or NA score
merged <- merged[!is.na(merged$OS) & !is.na(merged$score), ]

# Median split & group sizes
merged$group <- ifelse(merged$score > median(merged$score, na.rm = TRUE), "high", "low")
cat("Group counts:\n"); print(table(merged$group))

# If group sizes < 3 per group, warn
if (any(table(merged$group) < 3)) warning("One or both groups have fewer than 3 samples — results will be unstable.")

surv_obj <- Surv(time = merged$OS, event = merged$Vital.Status)
fit <- survfit(surv_obj ~ group, data = merged)

# 11. Plot (with p-value)
p <- ggsurvplot(
  fit,
  data = merged,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  ggtheme = theme_minimal(),
  palette = c("#E41A1C", "#377EB8"),
  risk.table.height = 0.25
)
print(p)
