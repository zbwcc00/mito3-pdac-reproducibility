options(stringsAsFactors = FALSE)

root_dir <- "D:/pdac_project"
input_dir <- file.path(root_dir, "incoming_data", "GSE240078")
result_dir <- file.path(root_dir, "data", "03_results", "feasibility", "gse240078")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(input_dir, "GSE240078_quantile_normalized_data_223samples.txt")
soft_file <- file.path(input_dir, "GSE240078_family.soft")
gene_audit_file <- file.path(root_dir, "data", "03_results", "feasibility", "pdac_mechanism_gene_set_audit.tsv")
stopifnot(file.exists(matrix_file), file.exists(soft_file), file.exists(gene_audit_file))

soft_lines <- readLines(soft_file, warn = FALSE)
starts <- grep("^!Sample_title = ", soft_lines)
ends <- c(starts[-1L] - 1L, length(soft_lines))
parse_value <- function(block, prefix) {
  x <- block[startsWith(block, prefix)]
  if (!length(x)) return(NA_character_)
  x <- sub(paste0("^", prefix), "", x[1L])
  x <- sub('^"', "", x)
  x <- sub('"\\s*$', "", x)
  x
}
meta_list <- lapply(seq_along(starts), function(i) {
  block <- soft_lines[starts[i]:ends[i]]
  title <- parse_value(block, "!Sample_title = ")
  geo <- parse_value(block, "!Sample_geo_accession = ")
  treatment <- parse_value(block, "!Sample_characteristics_ch1 = treatment: ")
  tissue_type <- parse_value(block, "!Sample_characteristics_ch1 = tissue type: ")
  matrix_id <- regmatches(title, regexpr("DSP-[^]]+", title))
  patient_id <- sub("\\..*$", "", sub(",.*$", "", title))
  data.frame(
    geo_accession = geo,
    title = title,
    matrix_id = paste0(matrix_id, ".dcc"),
    patient_id = patient_id,
    treatment = treatment,
    tissue_type = tissue_type,
    stringsAsFactors = FALSE
  )
})
metadata <- do.call(rbind, meta_list)

raw <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_symbol <- toupper(trimws(as.character(raw$GENE_ID)))
gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- toupper(trimws(as.character(raw$ID[is.na(gene_symbol) | gene_symbol == ""])))
expr <- as.matrix(raw[, -(1:2), drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- make.unique(gene_symbol)

matrix_samples <- colnames(expr)
metadata <- metadata[metadata$matrix_id %in% matrix_samples, , drop = FALSE]
metadata <- metadata[match(matrix_samples, metadata$matrix_id), , drop = FALSE]
metadata$matrix_id <- matrix_samples
metadata$group <- ifelse(tolower(metadata$treatment) == "naïve", "naive",
  ifelse(tolower(metadata$treatment) == "treated", "treated", "control"))
metadata$region <- tolower(metadata$tissue_type)

gene_audit <- read.delim(gene_audit_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_audit$gene_or_feature <- toupper(trimws(gene_audit$gene_or_feature))
gene_audit <- gene_audit[gene_audit$gene_or_feature != "" & !is.na(gene_audit$gene_or_feature), , drop = FALSE]
module_genes <- lapply(split(gene_audit$gene_or_feature, gene_audit$module), unique)

coverage <- do.call(rbind, lapply(names(module_genes), function(module_name) {
  genes <- module_genes[[module_name]]
  present <- intersect(genes, rownames(expr))
  data.frame(dataset = "GSE240078", module = module_name,
    n_prespecified = length(genes), n_present = length(present),
    coverage = length(present) / length(genes),
    present_genes = paste(present, collapse = ";"),
    missing_genes = paste(setdiff(genes, present), collapse = ";"),
    stringsAsFactors = FALSE)
}))
write.table(coverage, file.path(result_dir, "module_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

log_expr <- log2(expr + 1)
gene_z <- t(scale(t(log_expr)))
gene_z[is.na(gene_z)] <- 0
score_list <- list()
for (module_name in names(module_genes)) {
  present <- intersect(module_genes[[module_name]], rownames(gene_z))
  if (!length(present)) next
  score_list[[module_name]] <- data.frame(
    matrix_id = colnames(gene_z), module = module_name,
    score = colMeans(gene_z[present, , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE)
}
scores <- do.call(rbind, score_list)
scores <- merge(scores, metadata, by = "matrix_id", all.x = TRUE, sort = FALSE)
write.table(scores, file.path(result_dir, "roi_module_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

patient_scores <- aggregate(score ~ patient_id + treatment + group + region + module, data = scores, FUN = mean)
write.table(patient_scores, file.path(result_dir, "patient_region_module_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

test_rows <- list()
for (module_name in unique(patient_scores$module)) {
  for (region_name in unique(patient_scores$region)) {
    x <- patient_scores$score[patient_scores$module == module_name & patient_scores$region == region_name & patient_scores$group == "treated"]
    y <- patient_scores$score[patient_scores$module == module_name & patient_scores$region == region_name & patient_scores$group == "naive"]
    if (length(x) < 3 || length(y) < 3) next
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    test_rows[[length(test_rows) + 1L]] <- data.frame(
      dataset = "GSE240078", module = module_name, region = region_name,
      n_treated = length(x), n_naive = length(y),
      median_treated = median(x), median_naive = median(y),
      delta_median = median(x) - median(y), wilcox_p = wt$p.value,
      stringsAsFactors = FALSE)
  }
}
tests <- if (length(test_rows)) do.call(rbind, test_rows) else data.frame()
if (nrow(tests)) tests$fdr <- p.adjust(tests$wilcox_p, method = "BH")
write.table(tests, file.path(result_dir, "patient_level_treatment_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(metadata, file.path(result_dir, "roi_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  dataset = "GSE240078", n_matrix_rois = ncol(expr), n_metadata_rois = nrow(metadata),
  n_patients = length(unique(metadata$patient_id)),
  n_naive_patients = length(unique(metadata$patient_id[metadata$group == "naive"])),
  n_treated_patients = length(unique(metadata$patient_id[metadata$group == "treated"])),
  n_tumor_rois = sum(metadata$region == "tumor", na.rm = TRUE),
  n_stroma_rois = sum(metadata$region == "stroma", na.rm = TRUE),
  n_genes = nrow(expr), stringsAsFactors = FALSE)
write.table(summary, file.path(result_dir, "dataset_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(result_dir, "sessionInfo.txt"))
