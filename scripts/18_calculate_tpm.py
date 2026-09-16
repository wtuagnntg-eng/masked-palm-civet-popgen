#!/usr/bin/env python3

import argparse
import csv
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate gene-level TPM from featureCounts gene counts "
            "and effective non-redundant exon lengths."
        )
    )

    parser.add_argument(
        "--counts",
        required=True,
        help="featureCounts gene-level count file"
    )

    parser.add_argument(
        "--lengths",
        required=True,
        help="Gene effective exon length table"
    )

    parser.add_argument(
        "--sample-metadata",
        required=True,
        help="RNA-seq sample metadata TSV"
    )

    parser.add_argument(
        "--output-all",
        required=True,
        help="Output TPM table for all genes used in TPM normalization"
    )

    parser.add_argument(
        "--output-filtered",
        required=True,
        help="Output TPM table after total read count >10 filtering"
    )

    parser.add_argument(
        "--output-counts",
        required=True,
        help="Output filtered raw count table"
    )

    return parser.parse_args()


def load_lengths(path):
    lengths = {}

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        required = {
            "gene_id",
            "effective_exon_length_bp"
        }

        if not required.issubset(reader.fieldnames):
            raise SystemExit(
                "ERROR: gene length table must contain columns: "
                "gene_id and effective_exon_length_bp"
            )

        for row in reader:
            gene_id = row["gene_id"]
            length = float(row["effective_exon_length_bp"])

            if length <= 50:
                continue

            lengths[gene_id] = length

    return lengths


def load_sample_labels(path):
    """
    Map BAM basename to tissue/sample label.

    Expected metadata columns:
        sample_id
        tissue
        r1
        r2
    """
    label_by_bam = {}

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        if "sample_id" not in reader.fieldnames:
            raise SystemExit(
                "ERROR: sample metadata must contain a sample_id column"
            )

        for row in reader:
            sample_id = row["sample_id"].strip()

            if not sample_id:
                continue

            tissue = row.get("tissue", "").strip()
            label = tissue if tissue else sample_id

            bam_name = f"{sample_id}.sorted.bam"
            label_by_bam[bam_name] = label

    return label_by_bam


def read_featurecounts(path, lengths, label_by_bam):

    with open(path) as f:

        lines = [
            line
            for line in f
            if not line.startswith("#")
        ]

    reader = csv.reader(
        lines,
        delimiter="\t"
    )

    header = next(reader)

    if len(header) < 7:
        raise SystemExit(
            "ERROR: unexpected featureCounts format"
        )

    raw_sample_columns = header[6:]

    sample_labels = []

    for column in raw_sample_columns:

        base = os.path.basename(column)

        label = label_by_bam.get(
            base,
            os.path.splitext(
                os.path.splitext(base)[0]
            )[0]
        )

        sample_labels.append(label)

    if len(set(sample_labels)) != len(sample_labels):
        raise SystemExit(
            "ERROR: sample/tissue labels are not unique "
            "after metadata mapping"
        )

    genes = []

    for row in reader:

        if len(row) < 6 + len(sample_labels):
            continue

        gene_id = row[0]

        if gene_id not in lengths:
            continue

        counts = [
            float(x)
            for x in row[6:]
        ]

        genes.append(
            (
                gene_id,
                counts
            )
        )

    return sample_labels, genes


def calculate_tpm(genes, lengths, n_samples):

    rpk_records = []

    for gene_id, counts in genes:

        kb = lengths[gene_id] / 1000.0

        rpk = [
            count / kb
            for count in counts
        ]

        rpk_records.append(
            (
                gene_id,
                counts,
                rpk
            )
        )

    scale_factors = []

    for j in range(n_samples):

        total_rpk = sum(
            rpk[j]
            for _, _, rpk in rpk_records
        )

        scale = total_rpk / 1e6
        scale_factors.append(scale)

    tpm_records = []

    for gene_id, counts, rpk in rpk_records:

        tpm = []

        for j in range(n_samples):

            if scale_factors[j] > 0:
                value = rpk[j] / scale_factors[j]
            else:
                value = 0.0

            tpm.append(value)

        tpm_records.append(
            (
                gene_id,
                counts,
                tpm
            )
        )

    return tpm_records


def write_tpm(path, sample_labels, records):

    with open(path, "w") as out:

        out.write(
            "gene_id\t"
            + "\t".join(sample_labels)
            + "\n"
        )

        for gene_id, counts, tpm in records:

            out.write(
                gene_id
                + "\t"
                + "\t".join(
                    f"{x:.10g}"
                    for x in tpm
                )
                + "\n"
            )


def write_counts(path, sample_labels, records):

    with open(path, "w") as out:

        out.write(
            "gene_id\t"
            + "\t".join(sample_labels)
            + "\n"
        )

        for gene_id, counts, tpm in records:

            out.write(
                gene_id
                + "\t"
                + "\t".join(
                    f"{x:g}"
                    for x in counts
                )
                + "\n"
            )


def main():

    args = parse_args()

    ########################################
    # Load effective exon lengths
    ########################################

    lengths = load_lengths(
        args.lengths
    )

    if not lengths:
        raise SystemExit(
            "ERROR: no valid genes found in length table"
        )

    ########################################
    # Load sample metadata
    ########################################

    label_by_bam = load_sample_labels(
        args.sample_metadata
    )

    ########################################
    # Read gene-level featureCounts output
    ########################################

    sample_labels, genes = read_featurecounts(
        args.counts,
        lengths,
        label_by_bam
    )

    if not genes:
        raise SystemExit(
            "ERROR: no gene counts matched the gene-length table"
        )

    ########################################
    # Calculate TPM BEFORE low-expression filtering
    ########################################

    tpm_records = calculate_tpm(
        genes,
        lengths,
        len(sample_labels)
    )

    ########################################
    # Write TPM for all genes
    ########################################

    write_tpm(
        args.output_all,
        sample_labels,
        tpm_records
    )

    ########################################
    # Downstream low-expression filtering
    #
    # Manuscript criterion:
    # total read count >10 across seven tissues
    ########################################

    filtered_records = [
        record
        for record in tpm_records
        if sum(record[1]) > 10
    ]

    ########################################
    # Write filtered TPM and counts
    ########################################

    write_tpm(
        args.output_filtered,
        sample_labels,
        filtered_records
    )

    write_counts(
        args.output_counts,
        sample_labels,
        filtered_records
    )

    ########################################
    # Summary
    ########################################

    print(
        f"Genes with valid effective exon lengths: "
        f"{len(lengths):,}"
    )

    print(
        f"Genes present in featureCounts output: "
        f"{len(genes):,}"
    )

    print(
        f"Genes retained after total read count >10: "
        f"{len(filtered_records):,}"
    )

    print(
        "TPM normalization was performed before "
        "low-expression filtering."
    )


if __name__ == "__main__":
    main()
