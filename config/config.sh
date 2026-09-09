#!/usr/bin/env bash
set -euo pipefail

# Generic project configuration. Edit these values before running.
PROJECT_DIR="/path/to/project"
REF="${PROJECT_DIR}/reference/masked_palm_civet.fa"
GFF="${PROJECT_DIR}/reference/masked_palm_civet.gff3"
AUTOSOME_RENAME="${PROJECT_DIR}/metadata/autosomes.rename.tsv"
SAMPLES_WGS="${PROJECT_DIR}/metadata/wgs_samples.tsv"
SAMPLES_RNA="${PROJECT_DIR}/metadata/rna_samples.tsv"
TRANSCRIPT_TO_GENE="${PROJECT_DIR}/metadata/transcript_to_gene.tsv"
GENE2GO="${PROJECT_DIR}/metadata/gene2go.tsv"
GENE2KO="${PROJECT_DIR}/metadata/gene2ko.tsv"
KO2PATHWAY="${PROJECT_DIR}/metadata/ko2pathway.tsv"
IMMUNE_PANEL="${PROJECT_DIR}/metadata/immune_gene_panel.tsv"

THREADS=24
JAVA_MEM="32g"

WGS_DIR="${PROJECT_DIR}/wgs"
RNA_DIR="${PROJECT_DIR}/rnaseq"
RESULTS="${PROJECT_DIR}/results"

STRICT_VCF="${RESULTS}/variants/strict.autosomes.renamed.vcf.gz"
RELAXED_VCF="${RESULTS}/variants/relaxed.autosomes.renamed.vcf.gz"

# Population-structure analyses use STRICT_VCF. Observed heterozygosity,
# nucleotide diversity, LD, ROH and differentiation analyses use RELAXED_VCF
# unless stated otherwise.
