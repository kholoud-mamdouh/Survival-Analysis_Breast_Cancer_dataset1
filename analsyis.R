# ================================
# Survival Analysis using Gene Signature
# ================================

rm(list = ls())

library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(biomaRt)
library(dplyr)

# ---- 1. Load the data ----
clinical_data <- read.csv("clinical_data.csv", header = TRUE)
clinical_data$Case = gsub(pattern = "_", replacement = "M_", clinical_data$Case)
expression_matrix <- read.table("expression_data.txt", header = TRUE, sep = "\t", row.names = 1)
colnames(expression_matrix) <- sub("^X", "", colnames(expression_matrix))

# ---- 2. Subset clinical data (tumour type M only) ----
clincical_data_subset <- clinical_data[, c("Case", "Vital.Status", "SPBM", "OS")]


# Format survival variables
clincical_data_subset$OS <- as.numeric(clincical_data_subset$OS)
clincical_data_subset$Vital.Status <- ifelse(clincical_data_subset$Vital.Status == "Dead", 1, 0)

# ---- 3. Map Ensembl IDs → gene symbols ----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(expression_matrix)
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Remove blanks and duplicates
gene_annotations <- gene_annotations[gene_annotations$hgnc_symbol != "", ]
gene_annotations <- gene_annotations[!duplicated(gene_annotations$ensembl_gene_id), ]

# ---- 4. Collapse duplicates (mean expression per symbol) ----
expr_df <- cbind(ensembl_id = rownames(expression_matrix), expression_matrix)
expr_annotated <- merge(expr_df, gene_annotations, by.x = "ensembl_id", by.y = "ensembl_gene_id")

expr_collapsed <- expr_annotated %>%
  dplyr::select(-ensembl_id) %>%
  group_by(hgnc_symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE))


# Convert back to matrix
expression_matrix <- as.data.frame(expr_collapsed)
rownames(expression_matrix) <- expression_matrix$hgnc_symbol
expression_matrix$hgnc_symbol <- NULL
expression_matrix <- as.matrix(expression_matrix)

# ---- 5. Match samples between expression and clinical data ----

colnames(expression_matrix) <- gsub("\\.", "-", colnames(expression_matrix))

common_samples <- intersect(clincical_data_subset$Case, colnames(expression_matrix))

clincical_data_subset <- clincical_data_subset[clincical_data_subset$Case %in% common_samples, ]
expression_matrix <- expression_matrix[, common_samples]

# ---- 6. Define gene signature ----
signature_genes <- c("OLIG1", "OLIG2", "SOX10","PDGFRA","CSPG4","PTPRZ1","GPR17")   
# "OLIG1", "OLIG2", "SOX10","PDGFRA","CSPG4","PTPRZ1","GPR17"

# Keep only genes in the signature that exist in the expression matrix
signature_genes <- signature_genes[signature_genes %in% rownames(expression_matrix)]

if(length(signature_genes) == 0){
  stop("None of the signature genes were found in expression matrix!")
}

# ---- 7. Compute signature score (mean expression across genes) ----
signature_scores <- colMeans(expression_matrix[signature_genes, , drop = FALSE])

# Add to clinical data
clincical_data_subset$signature_score <- signature_scores[clincical_data_subset$Case]

# ---- 8. Split into high vs low (median cutoff) ----
clincical_data_subset$group <- ifelse(
  clincical_data_subset$signature_score >= median(clincical_data_subset$signature_score, na.rm = TRUE),
  "High", "Low"
)

# ---- 9. Kaplan–Meier survival analysis ----
surv_obj <- Surv(time = clincical_data_subset$OS, event = clincical_data_subset$Vital.Status)

fit <- survfit(surv_obj ~ group, data = clincical_data_subset)

# ---- 10. Plot ----
ggsurvplot(
  fit,
  data = clincical_data_subset,
  risk.table = TRUE,
  pval = TRUE,                # log-rank p-value
  conf.int = TRUE,            # CI
  xlab = "Time (days)",
  ylab = "SPBM Survival Probability",
  ggtheme = theme_minimal(),
  palette = c("#E41A1C", "#377EB8")
)

