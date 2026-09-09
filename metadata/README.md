# Metadata files

The TSV files in this directory are deliberately distributed as templates rather than with local filesystem paths or private sample metadata.

Before running the workflow, populate:
- `wgs_samples.tsv`: all WGS sample IDs, `Farm`/`Wild` ecotype, and paired FASTQ paths.
- `rna_samples.tsv`: the seven RNA-seq samples, tissue labels, and paired FASTQ paths.
- `transcript_to_gene.tsv`: transcript ID to the gene symbol/identifier used by the functional annotation.
- `gene2go.tsv`, `gene2ko.tsv`, `ko2pathway.tsv`: annotation mappings used for enrichment.
- `immune_gene_panel.tsv`: the complete predefined panel reported in Supplementary Table S4, with one of the seven categories used in the manuscript. Use `Viral-entry` for the viral-entry category.

`autosomes.rename.tsv` is study-specific and excludes `HIC_ASM_7`; verify it against the exact reference assembly before analysis.
