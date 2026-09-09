# Masked palm civet population-genomic and transcriptomic workflow

Reproducible analysis scripts corresponding to the revised manuscript on genomic differentiation between farmed and wild masked palm civets (*Paguma larvata*).

## Important setup

1. Edit `config/config.sh` and set `PROJECT_DIR`, the reference FASTA and GFF3 files, thread count, and Java memory.
2. Populate the template TSV files in `metadata/` with the study-specific sample information and annotation mappings. Local filesystem paths are intentionally not included in this public repository.
3. Verify `metadata/autosomes.rename.tsv` against the exact GCA_030068075.1 reference assembly used for the analysis. The supplied mapping excludes `HIC_ASM_7` from the autosomal set and retains 21 autosomes.
4. Run analysis scripts individually or use `bash run_pipeline.sh`. Figure-generation scripts can be run with `bash run_figures.sh` after the required analysis outputs have been generated.

## Metadata and annotation tables

The files in the `metadata/` directory are provided as templates for the sample information and annotation mappings required by the workflow. Before running the pipeline, users should populate these files with the corresponding project-specific information and resources.

The required metadata include WGS and RNA-seq sample information, transcript-to-gene mappings, GO and KEGG annotation mappings, and the predefined immune- and viral-entry-associated gene panel. Local filesystem paths and project-specific annotation resources are intentionally not distributed in the public repository.

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

This dataset is used for LD-pruned population-structure analyses.

Less-restrictive dataset
--remove-indels
--min-alleles 2
--max-alleles 2
--min-meanDP 5
--max-meanDP 200
--max-missing 0.8

No MAF or Hardy-Weinberg equilibrium filter is applied at this stage. This dataset is used for observed heterozygosity, nucleotide diversity, LD, ROH, and genomic differentiation. For LD analysis only, allele frequencies and genotype missingness are recalculated independently within each population before PopLDdecay analysis.

Population structure

05_population_structure.sh uses the strict dataset and performs:

PLINK LD pruning using --indep-pairwise 50 5 0.2.
PCA of all individuals and an additional wild-only PCA using the same LD-pruned marker set.
ADMIXTURE analyses for K = 2–5 with 10 independent runs per K and seeds 1–10, using the --cv option without explicitly specifying the number of cross-validation folds.
Selection of the replicate with the lowest cross-validation error for visualization, together with calculation of the mean cross-validation error across replicates.

The script reports the marker counts before and after LD pruning for direct reproducibility checks.

07_phylogeny_iqtree.sh converts the same LD-pruned marker set to an IUPAC nucleotide alignment, retains variable sites, and runs IQ-TREE 3 using MFP+ASC, 1,000 ultrafast bootstrap replicates, and 1,000 SH-aLRT replicates. No external outgroup is specified, and the inferred phylogeny should therefore be interpreted as unrooted.

Diversity, LD, and ROH
09_diversity_heterozygosity.sh: calculates individual observed heterozygosity from the less-restrictive autosomal dataset and performs a two-sided Wilcoxon rank-sum test between farmed and wild individuals.
10_nucleotide_diversity.sh: calculates descriptive 50-kb windows with a 10-kb step and non-overlapping 5-Mb genomic blocks. The paired Wild-minus-Farm difference in nucleotide diversity is evaluated using 10,000 block-bootstrap replicates. No additional population-specific missingness filter is applied for nucleotide-diversity analysis.
11_ld_decay.sh: independently filters each population to MAF ≥0.05 and genotype call rate ≥80%, followed by PopLDdecay using -MaxDist 300. summarize_ld.R reproduces the mean r² summary across distances of 0–30 kb used in the manuscript.
12_roh.sh: identifies ROH using the less-restrictive dataset without additional MAF filtering, Hardy-Weinberg equilibrium filtering, or LD pruning. Parameters match the revised Methods, and FROH is calculated as total autosomal ROH length divided by the physical length of the 21 autosomes.
Highly differentiated genomic regions

13_genomic_differentiation.sh calculates wild and farmed nucleotide diversity and weighted Weir-Cockerham FST in 50-kb windows with 10-kb steps directly from the less-restrictive dataset. Windows are matched among the wild π, farmed π, and FST outputs and are retained only when each corresponding estimate contains at least 10 informative SNPs.

For retained windows:

pi ratio = piWild / piFarm

Highly differentiated windows are defined as those simultaneously at or above the genome-wide 95th percentile of weighted FST and the 95th percentile of piWild/piFarm. Physically overlapping highly differentiated windows are subsequently merged into non-overlapping genomic regions.

14_overlap_annotation.sh retains transcript models showing any physical overlap with a merged highly differentiated region. No minimum percentage-overlap criterion or flanking distance is applied.

15_run_enrichment.sh maps transcript identifiers to genes through transcript_to_gene.tsv rather than inferring gene names directly from GFF3 attribute strings.

GO and KEGG enrichment analyses are performed using clusterProfiler::enricher. For GO analysis, genes with GO annotations are used as the background universe. For KEGG analysis, genes that can be mapped to pathways through the supplied KO-to-pathway mapping are used as the background universe.

pvalueCutoff=1 and qvalueCutoff=1 are used so that both nominal and multiple-testing-adjusted results are retained in the output even when no term passes FDR correction. Benjamini-Hochberg adjustment is applied, with a minimum gene-set size of three.

RNA-seq and TPM

16_rnaseq_mapping_counts.sh aligns clean paired-end reads using HISAT2 with --dta, sorts and indexes alignment files using SAMtools, and performs transcript-level read counting with featureCounts using -p -t exon -g Parent.

17_build_transcript_lengths.py calculates transcript lengths from CDS annotations. CDS records explicitly annotated as partial=true are excluded, remaining CDS lengths are summed by transcript, and transcripts with a summed CDS length ≤50 bp are removed.

18_calculate_tpm.py retains transcript features with a total read count >10 and a mean read count >1 across the seven tissues, calculates reads per kilobase using summed CDS length, and then calculates TPM. Where multiple transcript models are associated with the same gene name, transcript-level TPM values are averaged to obtain gene-level expression estimates. RNA metadata are used to replace featureCounts BAM-path column names with tissue labels.

20_transcriptome_descriptive_plots.R performs PCA using all retained genes with non-zero variance, generates a row-scaled heatmap of the 100 genes showing the greatest expression variance across tissues, matches the complete predefined immune- and viral-entry-associated gene panel to the detected transcriptomic genes, and generates separate immune and viral-entry heatmaps.

The transcriptomic workflow is descriptive only. It does not perform differential-expression testing and does not calculate the removed susceptibility index.

Main scripts
Script	Purpose
01_prepare_reference.sh	Build reference indices
02_wgs_qc_mapping.sh	WGS read QC, alignment, read-group assignment, and duplicate removal
02b_mapping_qc.sh	Per-sample depth, mapping, pairing, and coverage-threshold summaries
03_gatk_joint_calling.sh	Per-sample GVCF calling and chromosome-wise joint genotyping
04_filter_variants.sh	Generate strict and less-restrictive SNP datasets
05_population_structure.sh	LD pruning, PCA, wild-only PCA, and ADMIXTURE
06_vcf_to_iupac_fasta.py	Convert biallelic genotypes to a variable-site IUPAC alignment
07_phylogeny_iqtree.sh	IQ-TREE 3 phylogenetic analysis with ascertainment-bias correction
08_relatedness_king.sh	KING genome-wide relatedness analysis
09_diversity_heterozygosity.sh	Individual observed heterozygosity and Wilcoxon test
10_nucleotide_diversity.sh	Sliding-window π, paired 5-Mb blocks, and block bootstrap
11_ld_decay.sh	Population-specific filtering and PopLDdecay analysis
12_roh.sh	ROH classification and FROH calculation
13_genomic_differentiation.sh	Empirical FST/π-ratio outlier intersection and region merging
14_overlap_annotation.sh	Transcript overlap with merged differentiated regions
15_run_enrichment.sh	GO and KEGG enrichment using analysis-specific background universes
16_rnaseq_mapping_counts.sh	HISAT2 alignment and featureCounts quantification
17_build_transcript_lengths.py	Summed non-partial CDS lengths by transcript
18_calculate_tpm.py	Transcript- and gene-level TPM calculation
20_transcriptome_descriptive_plots.R	Transcriptomic PCA and heatmaps
21–24_plot_*.R	Reproducible plotting scripts for manuscript figures
Interpretation

The genomic outlier analysis identifies empirically highly differentiated regions and should not be interpreted as a demographic-null test of selection.

The seven-tissue RNA-seq dataset was generated from a single adult farmed female masked palm civet and is therefore analyzed descriptively. The transcriptomic results are not intended as statistical evidence of tissue-level differential expression or as functional validation of the population-genomic results.