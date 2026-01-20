# ==============================================================================
# LIPID GENOMICS PIPELINE
# Description: QC, Association Testing, and Visualization for Lipid Traits
# System: Windows/Mac/Linux compatible
# ==============================================================================

# ==============================================================================
# PART 1: DATA PREPARATION & POPULATION STRUCTURE
# ==============================================================================

# --- Setup & Dependencies ---
rm(list=ls())
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ---------------- USER-DEFINED PATHS ----------------
# Update these paths before running
plink_path      <- "/path/to/plink2"           # Path to PLINK2 executable
vcf_file        <- "data/raw_data.vcf"         # Input VCF file (if starting from VCF)
genotype_prefix <- "lipgen_data"               # Base name for generated PLINK files
phenotype_file  <- "data/phenotypes.csv"       # Clinical phenotype file
output_dir      <- "results/"                  # Directory for all outputs
# ----------------------------------------------------

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Detect PLINK2 (Fallback if user path is invalid)
if (!file.exists(plink_path)) {
  if (Sys.info()['sysname'] == "Windows") {
    plink_path <- "C:/Plink/plink2.exe"
  } else {
    plink_path <- "/opt/homebrew/bin/plink2"
  }
}
plink_path <- path.expand(plink_path)
if (!file.exists(plink_path)) stop("PLINK2 not found at: ", plink_path)

# --- Helper Functions ---

# Standardize IDs
normalize_id3 <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- sub("^#", "N", s)
  d <- suppressWarnings(as.integer(sub("^N", "", s)))
  ifelse(is.na(d), NA_character_, sprintf("N%03d", d %% 100L))
}

# Calculate Delta and Percent Change
chg <- function(post, pre) {
  list(delta = post - pre, pct = ifelse(is.finite(pre) & pre != 0, 100 * (post - pre) / pre, NA_real_))
}

# Format Mean (SD)
msd <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
}

# Paired T-test P-value
get_pval <- function(pre, post) {
  valid <- complete.cases(pre, post)
  if (sum(valid) < 3) return("NA")
  p <- t.test(pre[valid], post[valid], paired = TRUE)$p.value
  if (p < 0.001) return("<0.001") else return(sprintf("%.3f", p))
}

# Linear Model Wrapper
run_lm_clean <- function(delta_col, baseline_col, data, outcome_label) {
  covariates <- c(baseline_col, "Age", "PC1")
  if ("Sex" %in% names(data)) covariates <- c(covariates, "Sex")
  
  f <- as.formula(paste(delta_col, "~", paste(covariates, collapse=" + ")))
  fit <- lm(f, data = data)
  
  res <- as.data.table(summary(fit)$coefficients, keep.rownames = "Predictor")
  setnames(res, c("Predictor", "Beta", "SE", "t_stat", "P_Value"))
  res[, Outcome := outcome_label]
  res[, `:=`(Beta = round(Beta, 3), SE = round(SE, 3), P_Value = round(P_Value, 4))]
  return(res[, .(Outcome, Predictor, Beta, SE, P_Value)])
}

# Correlation Stats for Plotting
get_stats <- function(x, y, type = "corr") {
  keep <- complete.cases(x, y)
  if (type == "corr") {
    res <- cor.test(x[keep], y[keep], method = "spearman")
    return(sprintf("r = %.2f, p = %.3f", res$estimate, res$p.value))
  } else {
    fit <- lm(y[keep] ~ x[keep])
    return(sprintf("R² = %.2f", summary(fit)$r.squared))
  }
}

# --- Pipeline Execution ---

# 1. VCF Import
message("--- Step 1: VCF to PGEN ---")
out_raw <- file.path(output_dir, paste0(genotype_prefix, "_raw_auto"))

if (!file.exists(vcf_file)) warning("VCF file not found at: ", vcf_file)
system2(plink_path, args = c("--vcf", shQuote(vcf_file), "--autosome", "--make-pgen", "--out", out_raw, "--threads", "4"))

# 2. Phenotype Matching
message("--- Step 2: Load Phenotypes & Filter ---")
psam_raw <- fread(paste0(out_raw, ".psam"))
setnames(psam_raw, names(psam_raw)[1], "#IID") 

pheno_raw <- fread(phenotype_file, check.names = FALSE)
ph_clean <- copy(pheno_raw)
ph_clean[, IID := normalize_id3(ph_clean[[1]])]
psam_raw[, IID_norm := normalize_id3(`#IID`)]

keep_ids_norm <- intersect(ph_clean$IID, psam_raw$IID_norm)
keep_ids_raw <- psam_raw[IID_norm %chin% keep_ids_norm, `#IID`]

keep_file <- file.path(output_dir, "keep_raw_ids.txt")
fwrite(data.table(IID = keep_ids_raw), keep_file, sep = "\t", col.names = FALSE, quote = FALSE)
fwrite(psam_raw[IID_norm %chin% keep_ids_norm, .(Raw_ID = `#IID`, Fixed_ID = IID_norm)], file.path(output_dir, "MASTER_sample_list.csv"))

out_fixed <- file.path(output_dir, paste0(genotype_prefix, "_ids_fixed"))
system2(plink_path, args = c("--pfile", out_raw, "--keep", keep_file, "--make-pgen", "--out", out_fixed, "--threads", "4"))

# 3. Variant Cleaning
message("--- Step 3: Variant Cleaning ---")
pvar_file <- paste0(out_fixed, ".pvar")
all_lines <- readLines(pvar_file, warn = FALSE)
clean_lines <- all_lines[grep("^#CHROM", all_lines):length(all_lines)]

out_clean_base <- paste0(out_fixed, ".clean")
writeLines(clean_lines, paste0(out_clean_base, ".pvar"))
file.copy(paste0(out_fixed, ".pgen"), paste0(out_clean_base, ".pgen"), overwrite = TRUE)
file.copy(paste0(out_fixed, ".psam"), paste0(out_clean_base, ".psam"), overwrite = TRUE)

out_renamed <- file.path(output_dir, paste0(genotype_prefix, "_renamed"))
system2(plink_path, args = c("--pfile", out_clean_base, "--set-all-var-ids", "@:#:$r:$a", "--make-pgen", "--out", out_renamed, "--threads", "4"))

# 4. LD Pruning
message("--- Step 4: LD Pruning ---")
out_prune_list <- file.path(output_dir, "lipgen_prune_list")
system2(plink_path, args = c("--pfile", out_renamed, "--indep-pairwise", "200", "50", "0.1", "--out", out_prune_list, "--threads", "4"))

out_pruned <- file.path(output_dir, "lipgen_pruned")
system2(plink_path, args = c("--pfile", out_renamed, "--extract", paste0(out_prune_list, ".prune.in"), "--make-pgen", "--out", out_pruned, "--threads", "4"))

# 5. PCA Generation
message("--- Step 5: Generate PCA ---")
out_pca <- file.path(output_dir, "lipgen_pca_fixed")
system2(plink_path, args = c("--pfile", out_pruned, "--pca", "10", "--out", out_pca, "--threads", "4"))

# 6. Merge Data
message("--- Step 6: Merge Data ---")
PCs <- fread(paste0(out_pca, ".eigenvec"), header = FALSE)
if (ncol(PCs) == 12) setnames(PCs, c("#IID", "IID_raw", paste0("PC", 1:10))) else {
  setnames(PCs, c("IID_raw", paste0("PC", 1:10)))
  PCs[, `#IID` := IID_raw]
}
setcolorder(PCs, c("#IID", "IID_raw", paste0("PC", 1:10)))
PCs[, IID := normalize_id3(IID_raw)]

merged <- merge(PCs, ph_clean, by = "IID", all = FALSE)
merged <- merged[!(IID %in% c("#IID", "", NA))]
merged_dt <- as.data.table(copy(merged))

# 7. Variable Calculations
message("--- Step 7: Variable Calculations ---")
num_cols <- c("LDL1", "LDL2", "TAG1", "TAG2", "Age", "PC1")
for (cc in intersect(num_cols, names(merged_dt))) suppressWarnings(merged_dt[[cc]] <- as.numeric(merged_dt[[cc]]))

if ("Sex" %in% names(merged_dt)) {
  merged_dt[, Sex := tolower(trimws(as.character(Sex)))]
  merged_dt[Sex %in% c("male", "m"), Sex := "m"]
  merged_dt[Sex %in% c("female", "f"), Sex := "f"]
  merged_dt[, Sex := factor(Sex, levels = c("m", "f"))]
}

if (all(c("LDL1", "LDL2") %in% names(merged_dt))) { 
  z <- chg(merged_dt$LDL2, merged_dt$LDL1)
  merged_dt[, `:=`(dLDL = z$delta, pLDL = z$pct)] 
}
if (all(c("TAG1", "TAG2") %in% names(merged_dt))) { 
  z <- chg(merged_dt$TAG2, merged_dt$TAG1)
  merged_dt[, `:=`(dTAG = z$delta, pTAG = z$pct)] 
}

# Outlier Filter
if ("dTAG" %in% names(merged_dt)) merged_dt[dTAG > 200, dTAG := NA]

# 8. Kinship & Exclusion
message("--- Step 8: Kinship Analysis & Exclusion ---")
out_king <- file.path(output_dir, "lipgen_king_fixed")
system2(plink_path, args = c("--pfile", out_pruned, "--make-king-table", "--out", out_king, "--threads", "4"))

remove_id <- "N028" 
merged_dt_clean <- merged_dt[IID != remove_id]

out_final_csv <- file.path(output_dir, "lipgen_merged_deltas_final.csv")
fwrite(merged_dt_clean, out_final_csv)
message(sprintf("✅ Final Analytical N = %d", nrow(merged_dt_clean)))

# --- Output Generation (Descriptive) ---

# Table 1: Participant Characteristics
message("--- Generating Table 1 ---")
n_total <- nrow(merged_dt_clean)
n_female <- sum(merged_dt_clean$Sex == "f", na.rm=TRUE)
pct_female <- round(100 * n_female / n_total, 1)

demographics <- rbind(
  data.table(Variable = "Total N", Baseline = as.character(n_total), Post = "-", Delta = "-", P_Value = "-"),
  data.table(Variable = "Age (Years)", Baseline = msd(merged_dt_clean$Age), Post = "-", Delta = "-", P_Value = "-"),
  data.table(Variable = "Female Sex, n (%)", Baseline = sprintf("%d (%s%%)", n_female, pct_female), Post = "-", Delta = "-", P_Value = "-")
)

clin_rows <- list()
vars <- list(c("LDL-C (mg/dL)", "LDL1", "LDL2", "dLDL"), c("Triglycerides", "TAG1", "TAG2", "dTAG"))

for (v in vars) {
  name <- v[1]; pre <- merged_dt_clean[[v[2]]]; post <- merged_dt_clean[[v[3]]]; delta <- merged_dt_clean[[v[4]]]
  clin_rows[[length(clin_rows)+1]] <- data.table(
    Variable = name, Baseline = msd(pre), Post = msd(post), Delta = msd(delta), P_Value = get_pval(pre, post)
  )
}

table1_final <- rbind(demographics, rbindlist(clin_rows))
fwrite(table1_final, file.path(output_dir, "Table1_Descriptive_Stats.csv"))

# Table S1: Linear Models
message("--- Generating Table S1 (Linear Models) ---")
lm_results <- rbind(
  run_lm_clean("dLDL", "LDL1", merged_dt_clean, "dLDL"),
  run_lm_clean("dTAG", "TAG1", merged_dt_clean, "dTAG")
)
fwrite(lm_results, file.path(output_dir, "TableS1_Linear_Model_Coefficients.csv"))

# Figure S1: Population Structure
message("--- Generating Figure S1 ---")
plot_df <- copy(merged_dt_clean)
plot_df[, Sex := factor(Sex, levels = c("m", "f"), labels = c("Male", "Female"))]

# A: PC1 vs PC2
pA <- ggplot(plot_df, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 2.0, alpha = 0.7) +
  stat_ellipse(type = "t", linetype = 2, show.legend = FALSE) +
  labs(subtitle = "Genetic Ancestry (PC1 vs PC2)", x = "PC1", y = "PC2") +
  theme_classic() + theme(legend.position = "none", plot.subtitle = element_text(face="bold", size = 10))

# B: PC1 vs Age
label_B <- get_stats(plot_df$PC1, plot_df$Age, "corr")
pB <- ggplot(plot_df, aes(x = PC1, y = Age, color = Sex)) +
  geom_point(size = 2.0, alpha = 0.7) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", linetype = "dashed", se = FALSE, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = label_B, hjust = -0.1, vjust = 1.5, size = 2.8, fontface = "italic") +
  labs(subtitle = "Ancestry vs Age", x = "PC1", y = "Age (Years)") +
  theme_classic() + theme(legend.position = "none", plot.subtitle = element_text(face="bold", size = 10))

# C: PC1 vs dLDL
label_C <- get_stats(plot_df$PC1, plot_df$dLDL, "corr")
pC <- ggplot(plot_df, aes(x = PC1, y = dLDL, color = Sex)) +
  geom_point(size = 2.0, alpha = 0.7) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", linetype = "dashed", se = FALSE, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = label_C, hjust = -0.1, vjust = 1.5, size = 2.8, fontface = "italic") +
  labs(subtitle = "Ancestry vs ΔLDL-C", x = "PC1", y = "ΔLDL-C (mg/dL)") +
  theme_classic() + theme(legend.position = "none", plot.subtitle = element_text(face="bold", size = 10))

# D: PC1 vs dTAG
label_D <- get_stats(plot_df$PC1, plot_df$dTAG, "corr")
pD <- ggplot(plot_df, aes(x = PC1, y = dTAG, color = Sex)) +
  geom_point(size = 2.0, alpha = 0.7) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", linetype = "dashed", se = FALSE, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = label_D, hjust = -0.1, vjust = 1.5, size = 2.8, fontface = "italic") +
  labs(subtitle = "Ancestry vs ΔTAG", x = "PC1", y = "ΔTAG (mg/dL)") +
  theme_classic() + theme(legend.position = "none", plot.subtitle = element_text(face="bold", size = 10))

# E: Baseline vs dLDL
label_E <- get_stats(plot_df$LDL1, plot_df$dLDL, "lm")
pE <- ggplot(plot_df, aes(x = LDL1, y = dLDL, color = Sex)) +
  geom_point(size = 2.0, alpha = 0.7) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", linetype = "dashed", se = FALSE, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = label_E, hjust = -0.1, vjust = 1.5, size = 2.8, fontface = "italic") +
  labs(subtitle = "Baseline vs ΔLDL-C", x = "Baseline LDL-C (mg/dL)", y = "ΔLDL-C (mg/dL)") +
  theme_classic() + theme(legend.position = "none", plot.subtitle = element_text(face="bold", size = 10))

# F: Baseline vs dTAG
label_F <- get_stats(plot_df$TAG1, plot_df$dTAG, "lm")
pF <- ggplot(plot_df, aes(x = TAG1, y = dTAG, color = Sex)) +
  geom_point(size = 2.0, alpha = 0.7) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", linetype = "dashed", se = FALSE, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = label_F, hjust = -0.1, vjust = 1.5, size = 2.8, fontface = "italic") +
  labs(subtitle = "Baseline vs ΔTAG", x = "Baseline TAG (mg/dL)", y = "ΔTAG (mg/dL)") +
  theme_classic() + theme(plot.subtitle = element_text(face="bold", size = 10))

fig1_3x2 <- (pA | pB | pC) / (pD | pE | pF) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))

ggsave(file.path(output_dir, "Figure_S1_Population_Structure.png"), fig1_3x2, width = 14, height = 8, dpi = 600)

# Study Flow Table
n_raw_samp <- if(file.exists(paste0(out_raw, ".psam"))) nrow(fread(paste0(out_raw, ".psam"))) else NA
n_raw_var  <- if(file.exists(paste0(out_raw, ".pvar"))) length(readLines(paste0(out_raw, ".pvar"))) - length(grep("^#", readLines(paste0(out_raw, ".pvar")))) else NA
n_pheno_samp <- if(file.exists(paste0(out_fixed, ".psam"))) nrow(fread(paste0(out_fixed, ".psam"))) else NA
n_pheno_var  <- if(file.exists(paste0(out_clean_base, ".pvar"))) length(readLines(paste0(out_clean_base, ".pvar"))) - length(grep("^#", readLines(paste0(out_clean_base, ".pvar")))) else NA
n_pruned_var   <- if(file.exists(paste0(out_pruned, ".pvar"))) length(readLines(paste0(out_pruned, ".pvar"))) - length(grep("^#", readLines(paste0(out_pruned, ".pvar")))) else NA
n_final_model <- sum(complete.cases(merged_dt_clean[, .(Age, PC1, dLDL)]))

flow_table <- data.table(
  Stage = c("1. Initial Genotyped Cohort", "2. Phenotype Matching", "3. Relatedness Removed", "4. Final Analysis"),
  Samples = c(n_raw_samp, n_pheno_samp, nrow(merged_dt_clean), n_final_model),
  Variants_Used = c(n_raw_var, n_pheno_var, n_pruned_var, "N/A")
)
fwrite(flow_table, file.path(output_dir, "Study_Flow_Table.csv"))

message("✅ Part 1 Complete.")

# ==============================================================================
# PART 2: GWAS (Winsorized + FDR)
# ==============================================================================

# --- Setup & Dependencies ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(qqman)
  library(gridExtra)
  library(ggrepel)
})

# Winsorization Function
winsorize_3sd <- function(x) {
  if (all(is.na(x))) return(x)
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  lower <- mu - 3 * sigma
  upper <- mu + 3 * sigma
  x[x < lower] <- lower
  x[x > upper] <- upper
  return(x)
}

# --- GWAS Preparation ---
message("--- Step 1: Preparing Phenotypes & Running GWAS ---")

# Load analytical dataset from Part 1
merged_dt_clean <- fread(file.path(output_dir, "lipgen_merged_deltas_final.csv"))

# Genetic QC
out_clean_common <- file.path(output_dir, "lipgen_clean_common")
if (!file.exists(paste0(out_clean_common, ".pgen"))) {
  message("   > Running Genetic QC...")
  system2(plink_path, args = c(
    "--pfile", out_renamed,
    "--maf", "0.01",
    "--geno", "0.1",
    "--hwe", "1e-6",
    "--make-pgen",
    "--out", out_clean_common,
    "--threads", "4"
  ), stdout = FALSE, stderr = FALSE)
}

# --- LDL GWAS ---
gwas_pheno <- merged_dt_clean[, .(`#IID`, dLDL, Age, Sex, PC1)]
gwas_pheno[, Sex := as.numeric(factor(Sex, levels = c("m", "f")))]

message("   > Winsorizing dLDL at +/- 3SD...")
gwas_pheno[, dLDL := winsorize_3sd(dLDL)]

gwas_pheno[is.na(dLDL), dLDL := -99999]
gwas_pheno[is.na(Age), Age := -99999]
gwas_pheno[is.na(Sex), Sex := -99999]
gwas_pheno[is.na(PC1), PC1 := -99999] 

pheno_path_ldl <- file.path(output_dir, "lipgen_gwas_pheno.txt")
fwrite(gwas_pheno, pheno_path_ldl, sep = "\t", quote = FALSE)

message("   > Running GWAS for: dLDL")
out_gwas_ldl <- file.path(output_dir, "lipgen_gwas_results_LDL")
status <- system2(plink_path, args = c(
  "--pfile", out_clean_common,
  "--pheno", pheno_path_ldl,
  "--pheno-name", "dLDL",
  "--covar", pheno_path_ldl,
  "--covar-name", "Age,Sex",
  "--glm", "hide-covar", "allow-no-covars",
  "--input-missing-phenotype", "-99999",
  "--out", out_gwas_ldl,
  "--threads", "4"
), stdout = FALSE, stderr = FALSE)

if (status != 0) stop("GWAS failed for dLDL")

# --- TAG GWAS ---
message("--- Step 2: Running TAG GWAS ---")

psam_df <- read.table(paste0(out_clean_common, ".psam"), header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
colnames(psam_df)[1] <- "#IID"
psam_df$IID_norm <- normalize_id3(psam_df$`#IID`)

merged_df <- read.csv(file.path(output_dir, "lipgen_merged_deltas_final.csv"), stringsAsFactors = FALSE)
merged_df$IID_norm <- normalize_id3(merged_df$IID)

merged_ids <- merge(psam_df[, c("#IID", "IID_norm")], merged_df, by = "IID_norm")

final_tag_pheno <- data.frame(
  "#IID" = merged_ids$`#IID`,
  dTAG = merged_ids$dTAG,
  Age = merged_ids$Age,
  Sex = as.numeric(factor(merged_ids$Sex, levels = c("m", "f"))),
  PC1 = merged_ids$PC1,
  check.names = FALSE
)

message("   > Winsorizing dTAG at +/- 3SD...")
final_tag_pheno$dTAG <- winsorize_3sd(final_tag_pheno$dTAG)
final_tag_pheno[is.na(final_tag_pheno)] <- -99999

pheno_path_tag <- file.path(output_dir, "lipgen_gwas_pheno_tag.txt")
write.table(final_tag_pheno, pheno_path_tag, sep = "\t", quote = FALSE, row.names = FALSE)

message("   > Running GWAS for: dTAG")
out_gwas_tag <- file.path(output_dir, "lipgen_gwas_results_TAG")
status <- system2(plink_path, args = c(
  "--pfile", out_clean_common,
  "--pheno", pheno_path_tag,
  "--pheno-name", "dTAG",
  "--covar", pheno_path_tag,
  "--covar-name", "Age,Sex",
  "--glm", "hide-covar", "allow-no-covars",
  "--input-missing-phenotype", "-99999",
  "--out", out_gwas_tag,
  "--threads", "4"
), stdout = FALSE, stderr = FALSE)

if (status != 0) stop("GWAS failed for dTAG")

message("✅ GWAS Execution Complete.")

# --- Annotation & FDR ---
message("--- Step 3: Annotating Top Hits with FDR ---")

annotate_and_save <- function(plink_file, out_name) {
  if(!file.exists(plink_file)) return(NULL)
  
  res <- fread(plink_file)
  res <- res[TEST == "ADD" & !is.na(P)]
  
  res[, FDR := p.adjust(P, method = "BH")]
  top <- res[order(P)]
  
  fwrite(head(top, 100), paste0(out_name, ".csv"))
  
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    tryCatch({
      mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
      hits <- head(top, 100)
      hits[, query_region := paste(`#CHROM`, POS - 10, POS + 10, sep = ":")]
      
      ann <- biomaRt::getBM(attributes = c('chromosome_name','start_position','end_position','external_gene_name','description'),
                            filters = 'chromosomal_region', values = hits$query_region, mart = mart)
      setDT(ann)
      
      hits[, Gene := "Intergenic"]
      for (i in 1:nrow(hits)) {
        c <- as.character(hits$`#CHROM`[i]); p <- hits$POS[i]
        m <- ann[chromosome_name == c & start_position <= p & end_position >= p]
        if (nrow(m) > 0) hits$Gene[i] <- paste(unique(m$external_gene_name), collapse = "; ")
      }
      
      fwrite(hits[, .(Gene, `#CHROM`, POS, ID, P, FDR, BETA, REF, ALT)], paste0(out_name, "_annotated.csv"))
      message(paste("   > Saved annotated hits for:", out_name))
      
    }, error = function(e) message("   ! Annotation skipped (Network/Biomart error)"))
  }
}

annotate_and_save(paste0(out_gwas_ldl, ".dLDL.glm.linear"), file.path(output_dir, "lipgen_top_hits_dLDL"))
annotate_and_save(paste0(out_gwas_tag, ".dTAG.glm.linear"), file.path(output_dir, "lipgen_top_hits_dTAG"))

# ==============================================================================
# OFFLINE ANNOTATION UTILITIES
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_pkgs <- c("AnnotationDbi", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene", "GenomicRanges")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) BiocManager::install(new_pkgs, update=FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
})

# Robust Annotation Function
annotate_offline <- function(file_path) {
  if(!file.exists(file_path)) {
    message("   ! File not found: ", file_path)
    return(NULL)
  }
  
  message(paste("   > Annotating:", basename(file_path)))
  
  dt <- fread(file_path)
  if (!all(c("#CHROM", "POS") %in% names(dt))) return(NULL)
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  all_genes <- genes(txdb)
  
  chroms <- as.character(dt$`#CHROM`)
  if(!any(grepl("chr", chroms))) chroms <- paste0("chr", chroms)
  
  snp_ranges <- GRanges(seqnames = chroms, ranges = IRanges(start = dt$POS, end = dt$POS))
  
  # Find nearest genes
  hits_idx <- nearest(snp_ranges, all_genes, ignore.strand=TRUE)
  nearest_genes <- all_genes[hits_idx]
  dists <- distance(snp_ranges, nearest_genes, ignore.strand=TRUE)
  
  gene_symbols <- mapIds(org.Hs.eg.db, keys = nearest_genes$gene_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  
  dt$Gene <- ifelse(dists == 0, gene_symbols, paste0("Near: ", gene_symbols, " (", dists, "bp)"))
  dt[is.na(Gene), Gene := "Unknown"]
  
  fwrite(dt, file_path)
  message(paste("   ✅ Annotated and saved:", basename(file_path)))
}

annotate_offline(file.path(output_dir, "lipgen_top_hits_dLDL_annotated.csv"))
annotate_offline(file.path(output_dir, "lipgen_top_hits_dTAG_annotated.csv"))


# ==============================================================================
# PART 3: VISUALIZATION
# ==============================================================================
message("--- Step 4: Generating Global Plots ---")

# Manhattan Plot
plot_manhattan_clean <- function(file, title, filename) {
  if (!file.exists(file)) return(NULL)
  dt <- fread(file)[TEST == "ADD" & !is.na(P) & !is.na(POS)]
  dt <- dt[as.numeric(`#CHROM`) %in% 1:22]
  dt[, CHR := as.integer(`#CHROM`)]
  setorder(dt, CHR, POS)
  
  chr_max <- dt[, .(max_pos = max(POS)), by = CHR][order(CHR)]
  chr_max[, shift_pos := data.table::shift(cumsum(as.numeric(max_pos)), fill = 0)]
  dt <- merge(dt, chr_max[, .(CHR, shift_pos)], by = "CHR")
  dt[, BP_cum := POS + shift_pos]
  
  axis_pts <- dt[, .(center = min(BP_cum) + (max(BP_cum) - min(BP_cum))/2), by = CHR]
  
  dt[, neglog10P := -log10(pmax(P, 1e-300))]
  dt[, col_grp := fifelse(neglog10P > -log10(1e-5), "hit", as.character(CHR %% 2))]
  cols <- c("0" = "#2C3E50", "1" = "#ABB2B9", "hit" = "#E74C3C")
  
  p <- ggplot(dt, aes(x = BP_cum, y = neglog10P)) +
    geom_point(aes(color = col_grp), size = 0.6, alpha = 0.8) +
    scale_color_manual(values = cols, guide = "none") +
    geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_continuous(label = axis_pts$CHR, breaks = axis_pts$center, expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, subtitle = "Analysis of Autosomes 1-22", x = "Chromosome", y = expression(-log[10](italic(p)))) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", linewidth = 1))
  
  ggsave(filename, p, width = 10, height = 4.5, dpi = 600)
}

plot_manhattan_clean(paste0(out_gwas_ldl, ".dLDL.glm.linear"), "Exploratory GWAS of ΔLDL-C", file.path(output_dir, "Fig2_Manhattan_dLDL.png"))
plot_manhattan_clean(paste0(out_gwas_tag, ".dTAG.glm.linear"), "Exploratory GWAS of ΔTAG",   file.path(output_dir, "Fig3_Manhattan_dTAG.png"))

# QQ Plot
f_ldl <- paste0(out_gwas_ldl, ".dLDL.glm.linear")
f_tag <- paste0(out_gwas_tag, ".dTAG.glm.linear")

if (file.exists(f_ldl) && file.exists(f_tag)) {
  dt_ldl <- fread(f_ldl)[TEST=="ADD" & !is.na(P)]
  dt_tag <- fread(f_tag)[TEST=="ADD" & !is.na(P)]
  
  calc_lambda <- function(p) median(qchisq(1-p, df=1))/qchisq(0.5, df=1)
  qq_dat <- rbind(
    data.table(obs = -log10(sort(dt_ldl$P)), exp = -log10(ppoints(nrow(dt_ldl))), trait = "ΔLDL-C", lam = calc_lambda(dt_ldl$P)),
    data.table(obs = -log10(sort(dt_tag$P)), exp = -log10(ppoints(nrow(dt_tag))), trait = "ΔTAG",   lam = calc_lambda(dt_tag$P))
  )
  qq_dat[, label := sprintf("lambda[GC] == %.3f", lam)]
  
  p_qq <- ggplot(qq_dat, aes(x = exp, y = obs)) +
    geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed") +
    geom_point(size = 0.8, alpha = 0.6, color = "#2C3E50") +
    facet_wrap(~ trait, nrow = 1) +
    geom_text(data = unique(qq_dat[, .(trait, label)]), aes(x = 0, y = Inf, label = label), 
              hjust = 0, vjust = 1.5, parse = TRUE, inherit.aes = FALSE) +
    labs(x = expression(bold(Expected)~~-log[10](italic(p))), y = expression(bold(Observed)~~-log[10](italic(p)))) +
    theme_classic(base_size = 12) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, "Fig1_QQ_Combined.png"), p_qq, width = 8, height = 4, dpi = 600)
}

# ==============================================================================
# SECTION D: TARGETED PATHWAY ANALYSIS
# Filter: MAC >= 2
# ==============================================================================

rm(list=ls(pattern = "merged|gwas|clean")) # Clean up memory
suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
})

message("--- Section D: Targeted Pathway Analysis ---")

reactome_id <- "R-HSA-174824"
pheno <- fread(file.path(output_dir, "lipgen_merged_deltas_final.csv"))

if("IID_raw" %in% names(pheno)) {
  pheno[, IID_JOIN := as.character(IID_raw)] 
} else {
  pheno[, IID_JOIN := as.character(IID)]
}

# Winsorization (Redefined here for standalone safety)
winsorize_3sd <- function(x) {
  if (all(is.na(x))) return(x)
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  lower <- mu - (3 * sigma)
  upper <- mu + (3 * sigma)
  pmax(lower, pmin(x, upper))
}

if ("dLDL" %in% names(pheno)) pheno[, dLDL := winsorize_3sd(as.numeric(dLDL))]
if ("dTAG" %in% names(pheno)) pheno[, dTAG := winsorize_3sd(as.numeric(dTAG))]

# --- Query Ensembl ---
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
target_genes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                      filters = "reactome", values = reactome_id, mart = mart)
setDT(target_genes); setnames(target_genes, c("Gene", "Chr", "Start", "End"))
target_genes <- unique(target_genes[Chr %in% c(1:22, "X", "Y")], by="Gene")

# Simes Helper
calc_simes_p <- function(p_vals) {
  p_vals <- sort(p_vals[!is.na(p_vals)])
  n <- length(p_vals); if (n == 0) return(NA)
  min((n / 1:n) * p_vals)
}

# --- Analysis Loop ---
# Uses the 'renamed' file prefix defined in previous variables, constructed via output_dir
geno_prefix <- file.path(output_dir, paste0(genotype_prefix, "_renamed"))

for (target_pheno in c("dLDL", "dTAG")) {
  if (!target_pheno %in% names(pheno)) next
  
  message(sprintf("   > Analyzing: %s", target_pheno))
  results_table <- data.table()
  
  for (i in 1:nrow(target_genes)) {
    g <- target_genes[i]
    out_prefix <- file.path(output_dir, paste0("temp_", g$Gene))
    
    cmd <- paste(plink_path, "--pfile", geno_prefix, "--chr", g$Chr, "--from-bp", g$Start, "--to-bp", g$End, "--export A", "--out", out_prefix)
    system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)
    
    raw_file <- paste0(out_prefix, ".raw")
    if (!file.exists(raw_file)) next
    
    geno_data <- fread(raw_file)
    
    if ("IID" %in% names(geno_data)) {
      geno_data[, IID_JOIN := as.character(IID)] 
    } else {
      geno_data[, IID_JOIN := as.character(geno_data[[2]])]
    }
    
    analysis_set <- merge(pheno, geno_data, by="IID_JOIN")
    file.remove(Sys.glob(paste0(out_prefix, "*")))
    
    if(nrow(analysis_set) < 10) next 
    
    vars <- setdiff(names(geno_data), c("FID","IID","PAT","MAT","SEX","PHENOTYPE","IID_JOIN"))
    p_values <- numeric()
    
    for (v in vars) {
      v_vec <- analysis_set[[v]]
      
      if(sum(!is.na(v_vec)) < 10) next
      
      counts <- table(v_vec)
      if(length(counts) < 2) next
      
      minor_group_count <- min(counts)
      if(minor_group_count < 2) next
      
      form <- as.formula(paste(target_pheno, "~ `", v, "` + Age + Sex", sep=""))
      fit <- try(lm(form, data=analysis_set), silent=TRUE)
      
      if (!inherits(fit, "try-error")) {
        cf <- summary(fit)$coefficients
        idx <- grep(v, rownames(cf), fixed=TRUE)
        if(length(idx) > 0) p_values <- c(p_values, cf[idx, "Pr(>|t|)"])
      }
    }
    
    if (length(p_values) > 0) {
      results_table <- rbind(results_table, data.table(Gene = g$Gene, N_Vars = length(p_values), Simes_P = calc_simes_p(p_values)))
    }
  }
  
  if (nrow(results_table) > 0) {
    results_table[, FDR := p.adjust(Simes_P, method="BH")]
    fname <- file.path(output_dir, paste0("Reactome_", target_pheno, "_Winsor_Strict.csv"))
    fwrite(results_table[order(Simes_P)], fname)
    message(sprintf("   ✅ Saved %s (Top: %s, P=%.3e)", basename(fname), results_table$Gene[1], results_table$Simes_P[1]))
  } else {
    message(sprintf("   ⚠️ No hits for %s.", target_pheno))
  }
}

# ==============================================================================
# SECTION E: REGIONAL VISUALIZATION (LOCUSZOOM + BOXPLOTS)
# ==============================================================================

rm(list=ls(pattern = "geno_data|analysis_set"))
suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(grid)
  library(scales)
})

# --- Configuration ---
targets <- data.table(
  Gene  = c("APOC3", "APOC2", "AP2A1", "AP2A2", "A2M", "APOB"),
  Pheno = c("dLDL",  "dLDL",  "dLDL",  "dLDL",  "dLDL", "dTAG")
)

window_size <- 200000

pheno <- fread(file.path(output_dir, "lipgen_merged_deltas_final.csv"))
if ("IID_raw" %in% names(pheno)) {
  pheno[, IID_JOIN := as.character(IID_raw)]
} else {
  pheno[, IID_JOIN := as.character(IID)]
}

pheno[, dLDL_Plot := winsorize_3sd(as.numeric(dLDL))]
pheno[, dTAG_Plot := winsorize_3sd(as.numeric(dTAG))]

# Ensembl Mapping
mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "grch37.ensembl.org") 

gene_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = targets$Gene,
  mart       = mart
)
setDT(gene_coords)
setnames(gene_coords, c("Gene", "Chr", "Start", "End"))

gene_coords <- gene_coords[Chr %in% c(as.character(1:22), "X", "Y")]
gene_coords[, Span := End - Start]
setorder(gene_coords, Gene, -Span, -End)
gene_coords <- gene_coords[, .SD[1], by = Gene][, Span := NULL]

# --- Styling ---
ref_style <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face="bold", hjust=0, size=12),
    axis.title = element_text(face="bold", size=10),
    axis.text  = element_text(color="black", size=9),
    axis.line  = element_line(linewidth=0.5, color="black"),
    plot.margin = margin(2, 6, 2, 6)
  )

plot_list <- list()

# --- Plot Generation Loop ---
geno_prefix <- file.path(output_dir, paste0(genotype_prefix, "_renamed"))

for (i in 1:nrow(targets)) {
  
  t_gene <- targets$Gene[i]
  t_pheno_name <- targets$Pheno[i]
  t_pheno_col <- if (t_pheno_name == "dLDL") "dLDL_Plot" else "dTAG_Plot"
  
  y_label_box <- if (t_pheno_name == "dLDL") expression(paste(Delta, "LDL")) else expression(paste(Delta, "TAG"))
  box_limits <- if (t_pheno_name == "dLDL") c(-120, 0) else c(-150, 50)
  lz_legend_pos <- if (t_gene == "APOC2") c(0.85, 0.75) else "none"
  
  g_info <- gene_coords[Gene == t_gene]
  if (nrow(g_info) == 0) next
  
  target_chr <- as.character(g_info$Chr[1])
  center_bp  <- floor((g_info$Start[1] + g_info$End[1]) / 2)
  start_bp   <- max(1, center_bp - window_size)
  end_bp     <- center_bp + window_size
  
  out_prefix <- file.path(output_dir, paste0("temp_reg_", t_gene))
  
  cmd <- paste(
    shQuote(plink_path), "--pfile", geno_prefix, "--chr", target_chr,
    "--from-bp", start_bp, "--to-bp", end_bp, "--export A", "--out", out_prefix
  )
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  raw_file <- paste0(out_prefix, ".raw")
  if (!file.exists(raw_file)) { file.remove(Sys.glob(paste0(out_prefix, "*"))); next }
  
  geno_data <- fread(raw_file)
  if ("IID" %in% names(geno_data)) {
    geno_data[, IID_JOIN := as.character(IID)]
  } else {
    geno_data[, IID_JOIN := as.character(geno_data[[2]])]
  }
  
  analysis_set <- merge(pheno, geno_data, by = "IID_JOIN")
  if (nrow(analysis_set) == 0) { file.remove(Sys.glob(paste0(out_prefix, "*"))); next }
  
  variant_cols <- setdiff(names(geno_data), c("FID","IID","PAT","MAT","SEX","PHENOTYPE","IID_JOIN"))
  if (length(variant_cols) == 0) { file.remove(Sys.glob(paste0(out_prefix, "*"))); next }
  
  results <- data.table()
  dosages <- matrix(NA_real_, nrow = nrow(analysis_set), ncol = length(variant_cols))
  colnames(dosages) <- variant_cols
  
  for (k in seq_along(variant_cols)) {
    v <- variant_cols[k]
    v_vec <- analysis_set[[v]]
    if (sum(!is.na(v_vec)) < 10) next
    counts <- table(v_vec)
    if (length(counts) < 2 || min(counts) < 2) { dosages[, k] <- NA_real_; next }
    dosages[, k] <- ifelse(is.na(v_vec), mean(v_vec, na.rm = TRUE), v_vec)
    
    if (!all(c("Age","Sex", t_pheno_col) %in% names(analysis_set))) next
    form <- as.formula(paste0(t_pheno_col, " ~ `", v, "` + Age + Sex"))
    fit <- try(lm(form, data = analysis_set), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cf <- summary(fit)$coefficients
      idx <- grep(v, rownames(cf), fixed = TRUE)
      if (length(idx) > 0) results <- rbind(results, data.table(SNP = v, P = cf[idx, "Pr(>|t|)"]))
    }
  }
  
  if (nrow(results) == 0) { file.remove(Sys.glob(paste0(out_prefix, "*"))); next }
  
  results[, NegLogP := -log10(P)]
  results[, Pos := as.integer(tstrsplit(SNP, ":", fixed = TRUE)[[2]])]
  results <- results[!is.na(Pos)]
  if (nrow(results) == 0) { file.remove(Sys.glob(paste0(out_prefix, "*"))); next }
  
  top_hit <- results[which.min(P)]
  target_snp_id <- top_hit$SNP
  
  top_idx <- which(variant_cols == target_snp_id)
  r2_vals <- rep(0, length(variant_cols))
  if (length(top_idx) > 0 && !all(is.na(dosages[, top_idx]))) {
    r2_vals <- apply(dosages, 2, function(col) {
      if (all(is.na(col))) return(0)
      suppressWarnings(cor(dosages[, top_idx], col, use = "pairwise.complete.obs")^2)
    })
    r2_vals[is.na(r2_vals)] <- 0
  }
  
  plot_data <- merge(results, data.table(SNP = variant_cols, R2 = r2_vals), by = "SNP", all.x = TRUE)
  plot_data[is.na(R2), R2 := 0]
  
  # Regional Scatter Plot
  pA <- ggplot(plot_data, aes(x = Pos, y = NegLogP, fill = R2)) +
    geom_point(shape = 21, size = 2.5, alpha = 0.8, color = "black", stroke = 0.2) +
    scale_fill_gradientn(
      colors = c("#1e1e2c", "#8edbc3", "#eaa239", "#832b29"),
      limits = c(0, 1),
      name   = expression(italic(R)^2)
    ) +
    geom_point(data = plot_data[SNP == target_snp_id], shape = 23, size = 3.5, fill = "white", color = "black", stroke = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    coord_cartesian(ylim = c(0, 4.1), xlim = c(start_bp, end_bp)) + 
    labs(y = expression(-log[10](italic(p))), x = NULL) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = " Mb")) +
    ref_style +
    theme(
      legend.position = lz_legend_pos,
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(size=7),
      legend.text  = element_text(size=10),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  # Boxplot
  box_df <- analysis_set[!is.na(get(target_snp_id))]
  box_df[, Genotype := factor(round(get(target_snp_id)), levels = c(0,1,2), labels = c("0","1","2"))]
  
  pB <- ggplot(box_df, aes(x = Genotype, y = .data[[t_pheno_col]])) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5, color = "black", fill = "grey88") +
    geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.45, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white", color = "black") +
    coord_cartesian(ylim = box_limits) +
    labs(x = "Genotype", y = y_label_box, title = t_gene) +
    ref_style +
    theme(legend.position = "none")
  
  combined <- (pB | pA) + plot_layout(widths = c(1, 2))
  plot_list[[t_gene]] <- combined
  
  file.remove(Sys.glob(paste0(out_prefix, "*")))
  message(paste("   Processed:", t_gene))
}

# --- Final Assembly ---
if (length(plot_list) > 0) {
  final_grid <- wrap_plots(plot_list, ncol = 2)
  ggsave(file.path(output_dir, "Final_Grid_SpecificScales.png"), final_grid, width = 14, height = 10, dpi = 300)
  message("✅ Final_Grid_SpecificScales.png saved successfully.")
} else {
  message("⚠️ No plots generated.")
}