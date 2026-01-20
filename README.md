Lipid Genomics Analysis Pipeline
This repository contains the analytical pipeline used for the study:

"Genetic Modulation of LDL-C and Triglyceride Response to Mediterranean-Style Diet: Whole-Exome Evidence with a Focus on Plasma Lipoprotein Pathways". Iordanishvili S., et al.

Description
This pipeline performs end-to-end analysis for lipid genomic traits (LDL-C and Triglycerides), including:

Data Preparation: VCF to PLINK conversion, phenotype matching, and QC.

Population Structure: PCA generation and kinship analysis.

GWAS: Linear regression (winsorized phenotypes) with Benjamini-Hochberg FDR correction.

Visualization: Generation of Manhattan plots, QQ plots, and regional association plots.

Files
LIPID_GENOMICS_PIPELINE.R : Complete analysis script containing all processing steps.

Requirements
R >= 4.1.0

PLINK v2.0 (Must be installed and accessible)

R Packages:

CRAN: data.table, ggplot2, patchwork, qqman, gridExtra, ggrepel, scales, BiocManager

Bioconductor: biomaRt, AnnotationDbi, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg38.knownGene, GenomicRanges

Usage
Open LIPID_DIET_GWAS_ANALYSIS.R.

Edit the "USER-DEFINED PATHS" section at the top of the script to point to your local PLINK executable and input data.

Run the script source in R.

Data
No individual-level genomic or clinical data are included in this repository to protect participant privacy. Users must supply their own VCF and phenotype files formatted as described in the script comments.
