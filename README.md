# Masked palm civet population-genomic and transcriptomic workflow

Reproducible analysis scripts corresponding to the revised manuscript on genomic differentiation between farmed and wild masked palm civets (*Paguma larvata*).

## Important setup

1. Edit `config/config.sh` and set `PROJECT_DIR`, the reference FASTA and GFF3 files, thread count, and Java memory.
2. Populate the template TSV files in `metadata/` with the study-specific sample information and annotation mappings. Local filesystem paths are intentionally not included in this public repository.
3. Verify `metadata/autosomes.rename.tsv` against the exact GCA_030068075.1 reference assembly used for the analysis. The supplied mapping excludes `HIC_ASM_7` from the autosomal set and retains 21 autosomes.
4. Run analysis scripts individually or use:

```bash
bash run_pipeline.sh
```

Figure-generation scripts for Figures 1–4 can subsequently be run using:

```bash
bash run_figures.sh
```

Figure 5 and Figure S2 are generated as part of the transcriptomic workflow.

## Metadata and annotation tables

The files in the `metadata/` directory are provided as templates for the sample information and annotation mappings required by the workflow.

Before running the pipeline, users should populate these files with the corresponding project-specific information and resources.

The required metadata include WGS and RNA-seq sample information, autosomal chromosome-renaming information, GO and KEGG annotation mappings, and the predefined immune- and viral-entry-associated gene panel.

Local filesystem paths and project-specific annotation resources are intentionally not distributed in the public repository.

## Software

The workflow uses fastp 0.22.0; BWA-MEM2; SAMtools; Picard; GATK 4.2.0.0; BCFtools; VCFtools 0.1.13; PLINK 1.9; ADMIXTURE 1.3; IQ-TREE 3; KING; PopLDdecay; bedtools; HISAT2 2.2.1; featureCounts/Subread; Python 3; and R with `clusterProfiler`, `ggplot2`, `patchwork`, and `pheatmap`.

## Variant datasets

`04_filter_variants.sh` generates two autosomal biallelic SNP datasets.

### Strict/common-variant dataset

```text
--remove-indels
--min-alleles 2
--max-alleles 2
--min-meanDP 5
--max-meanDP 200
--max-missing 0.8
--maf 0.05
--hwe 0.001
```

This dataset is used for LD-pruned population-structure analyses.

### Less-restrictive dataset

```text
--remove-indels
--min-alleles 2
--max-alleles 2
--min-meanDP 5
--max-meanDP 200
--max-missing 0.8
```

No MAF or Hardy-Weinberg equilibrium filter is applied at this stage.

This dataset is used for observed heterozygosity, nucleotide diversity, LD, ROH, and genomic differentiation.

For LD analysis only, allele frequencies and genotype missingness are recalculated independently within each population before PopLDdecay analysis.

## Population structure

`05_population_structure.sh` uses the strict dataset and performs:

- PLINK LD pruning using `--indep-pairwise 50 5 0.2`.
- PCA of all individuals.
- An additional wild-only PCA using the LD-pruned autosomal marker set.
- ADMIXTURE analyses for K = 2–5 with 10 independent runs per K using seeds 1–10.
- Cross-validation using the `--cv` option.
- Selection of the replicate with the lowest cross-validation error for visualization.
- Calculation of the mean cross-validation error across replicates.

The script also reports marker counts before and after LD pruning.

## Phylogenetic analysis

`07_phylogeny_iqtree.sh` converts the LD-pruned SNP dataset to an IUPAC nucleotide alignment and retains variable sites.

IQ-TREE 3 is run using `MFP+ASC`, 1,000 ultrafast bootstrap replicates, and 1,000 SH-aLRT replicates.

No external outgroup is specified. The resulting phylogeny is therefore interpreted as unrooted.

## Relatedness

`08_relatedness_king.sh` estimates genome-wide pairwise relatedness using KING based on the LD-pruned autosomal SNP dataset.

Because the sampled animals show pronounced population structure, KING coefficients are interpreted primarily as relative measures of genomic similarity rather than definitive pedigree assignments.

## Diversity, LD, and ROH

### Observed heterozygosity

`09_diversity_heterozygosity.sh` calculates individual observed heterozygosity from the less-restrictive autosomal SNP dataset.

Farmed and wild individuals are compared using a two-sided Wilcoxon rank-sum test.

### Nucleotide diversity

`10_nucleotide_diversity.sh` calculates descriptive nucleotide diversity in 50-kb windows with a 10-kb step and nucleotide diversity in non-overlapping 5-Mb genomic blocks.

The paired Wild-minus-Farm nucleotide-diversity difference is evaluated using 10,000 block-bootstrap replicates.

No additional population-specific missingness filter is applied for nucleotide-diversity analysis.

### LD decay

`11_ld_decay.sh` independently filters each population to MAF ≥ 0.05 and genotype call rate ≥ 80%, followed by PopLDdecay using `-MaxDist 300`.

`summarize_ld.R` calculates the mean r² summary across distances of 0–30 kb used in the manuscript.

### Runs of homozygosity

`12_roh.sh` identifies ROH using the less-restrictive autosomal SNP dataset without additional MAF filtering, Hardy-Weinberg equilibrium filtering, or LD pruning.

ROH parameters match those reported in the revised Methods.

FROH is calculated as total autosomal ROH length divided by the total physical length of the 21 autosomes.

ROH are also classified according to their physical lengths.

## Highly differentiated genomic regions

`13_genomic_differentiation.sh` calculates wild and farmed nucleotide diversity and weighted Weir-Cockerham FST in 50-kb windows with a 10-kb step directly from the less-restrictive SNP dataset.

Windows from the wild π, farmed π, and FST outputs are matched by genomic coordinates.

Only windows containing at least 10 informative SNPs in each corresponding estimate are retained.

For retained windows:

```text
pi ratio = piWild / piFarm
```

Highly differentiated windows are defined empirically as those simultaneously falling within the upper 5% tails of weighted FST and piWild/piFarm.

Physically overlapping highly differentiated windows are subsequently merged into non-overlapping genomic regions.

`14_overlap_annotation.sh` retains transcript models showing any physical overlap with a merged highly differentiated genomic region.

No minimum percentage-overlap criterion or additional flanking distance is applied.

## GO and KEGG enrichment

`15_run_enrichment.sh` performs GO and KEGG enrichment analyses using `clusterProfiler::enricher`.

For GO enrichment, genes with GO annotations are used as the background universe.

For KEGG enrichment, genes that can be mapped to pathways through the supplied KO-to-pathway mapping are used as the background universe.

The analyses use:

```text
pvalueCutoff = 1
qvalueCutoff = 1
minGSSize = 3
```

These permissive output thresholds allow both nominal and multiple-testing-adjusted results to be retained even when no category passes FDR correction.

Benjamini-Hochberg adjustment is applied for multiple testing.

## RNA-seq and gene-level TPM

The transcriptomic workflow is coordinated by:

```text
19_transcriptome_expression.sh
```

The workflow runs in the following order:

```text
17_build_gene_annotation.py
        ↓
16_rnaseq_mapping_counts.sh
        ↓
18_calculate_tpm.py
        ↓
20_transcriptome_descriptive_plots.R
```

### Gene-level exon annotation

`17_build_gene_annotation.py` parses the reference GFF3 annotation and uses mRNA-to-gene relationships to assign exon features to genes.

Exonic intervals belonging to alternative transcript isoforms of the same gene are combined, and overlapping or directly adjacent intervals are merged to generate non-redundant gene-level exon intervals.

The script generates:

```text
gene_exons.saf
gene_lengths.tsv
gene_map.tsv
```

`gene_exons.saf` contains the non-redundant exon intervals used for gene-level featureCounts analysis.

For each gene, effective gene length is defined as the total physical length of the merged non-redundant exon intervals.

Genes with an effective exon length ≤50 bp are excluded.

`gene_map.tsv` records gene identifiers together with available gene names or symbols for downstream visualization.

### RNA-seq alignment and gene-level counting

`16_rnaseq_mapping_counts.sh` aligns clean paired-end RNA-seq reads to the masked palm civet reference genome using HISAT2 with the `--dta` option.

Alignment files are coordinate-sorted and indexed using SAMtools.

Gene-level read counts are generated using featureCounts in paired-end mode using the SAF annotation produced by `17_build_gene_annotation.py`.

Counting is therefore performed directly at the gene level rather than by averaging transcript-level expression estimates.

### TPM calculation

`18_calculate_tpm.py` calculates gene-level TPM values from gene-level read counts and effective non-redundant exon lengths.

For each gene:

```text
RPK = read count / effective exon length in kilobases
```

For each tissue:

```text
TPM = RPK / sum(RPK) × 1,000,000
```

TPM normalization is performed using all genes with valid effective exon lengths before downstream low-expression filtering.

After TPM normalization, genes are retained for exploratory transcriptomic analyses when:

```text
total read count across the seven tissues > 10
```

No additional mean-read-count threshold is applied.

The main outputs are:

```text
gene_TPM_all.tsv
gene_TPM.tsv
gene_counts_filtered.tsv
```

`gene_TPM_all.tsv` contains TPM values before downstream low-expression filtering.

`gene_TPM.tsv` contains genes retained after applying the total-read-count criterion and is used for downstream descriptive analyses.

## Transcriptomic descriptive analyses

`20_transcriptome_descriptive_plots.R` performs descriptive analysis of the seven-tissue transcriptomic dataset.

Expression values are transformed as:

```text
log2(TPM + 1)
```

PCA is performed using all retained genes showing non-zero variance across tissues.

The 100 genes showing the greatest variance in log2(TPM + 1) expression are used to visualize broader tissue-specific transcriptional variation.

Heatmap expression values are standardized independently by gene across tissues.

The predefined immune- and viral-entry-associated gene panel is matched to detected gene loci using gene identifiers and gene-name mappings.

Where more than one annotated gene locus is assigned the same gene symbol, the individual loci are retained separately rather than averaged.

For visualization, duplicated symbols are distinguished using numerical suffixes such as `CTSL-1` and `CTSL-2`.

The script generates:

```text
Figure5A_PCA.pdf
Figure5B_immune_heatmap.pdf
Figure5C_viral_entry_heatmap.pdf
FigureS2_top100_heatmap.pdf
```

The transcriptomic workflow is descriptive only.

No differential-expression testing is performed.

The transcriptomic data are derived from seven tissues collected from a single adult farmed female masked palm civet and are therefore not interpreted as evidence of farmed-versus-wild expression differences, statistical evidence of tissue-level differential expression, tissue susceptibility to infection, or functional validation of the population-genomic differentiation results.

## Main scripts

| Script | Purpose |
|---|---|
| `01_prepare_reference.sh` | Build reference indices |
| `02_wgs_qc_mapping.sh` | WGS read QC, alignment, read-group assignment, and duplicate removal |
| `02b_mapping_qc.sh` | Per-sample depth, mapping, pairing, and genome-coverage summaries |
| `03_gatk_joint_calling.sh` | Per-sample GVCF calling and chromosome-wise joint genotyping |
| `04_filter_variants.sh` | Generate strict and less-restrictive autosomal SNP datasets |
| `05_population_structure.sh` | LD pruning, PCA, wild-only PCA, and ADMIXTURE |
| `06_vcf_to_iupac_fasta.py` | Convert biallelic genotypes to a variable-site IUPAC alignment |
| `07_phylogeny_iqtree.sh` | IQ-TREE 3 phylogenetic analysis with ascertainment-bias correction |
| `08_relatedness_king.sh` | KING genome-wide relatedness analysis |
| `09_diversity_heterozygosity.sh` | Individual observed heterozygosity and Wilcoxon test |
| `10_nucleotide_diversity.sh` | Sliding-window π, paired 5-Mb blocks, and block bootstrap |
| `11_ld_decay.sh` | Population-specific filtering and PopLDdecay analysis |
| `12_roh.sh` | ROH classification and FROH calculation |
| `13_genomic_differentiation.sh` | Empirical FST/π-ratio differentiated-region scan |
| `14_overlap_annotation.sh` | Transcript overlap with merged differentiated regions |
| `15_run_enrichment.sh` | GO and KEGG enrichment using analysis-specific background universes |
| `16_rnaseq_mapping_counts.sh` | HISAT2 alignment and direct gene-level featureCounts quantification |
| `17_build_gene_annotation.py` | Generate non-redundant gene-level exon SAF, effective gene lengths, and gene-name mappings |
| `18_calculate_tpm.py` | Gene-level TPM calculation followed by low-expression filtering |
| `19_transcriptome_expression.sh` | Coordinate the complete transcriptomic workflow |
| `20_transcriptome_descriptive_plots.R` | Transcriptomic PCA and immune/viral-entry heatmaps |
| `21–24_plot_*.R` | Reproducible plotting scripts for manuscript Figures 1–4 |

## Interpretation

The genome-wide outlier analysis identifies empirically highly differentiated genomic regions.

These regions should not be interpreted as a demographic-null test of natural or artificial selection because extreme differentiation may also arise from founder effects, genetic drift, demographic history, geographic population structure, and relatedness.

The seven-tissue RNA-seq dataset was generated from a single adult farmed female masked palm civet.

The transcriptomic analysis is therefore descriptive and is not intended as statistical evidence of population-level or tissue-level differential expression, nor as functional validation of the population-genomic results.
