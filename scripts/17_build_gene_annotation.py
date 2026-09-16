#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
from urllib.parse import unquote


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build non-redundant gene-level exon intervals from a GFF3 file. "
            "Outputs a SAF file for featureCounts, effective gene exon lengths, "
            "and a gene ID-to-name mapping table."
        )
    )
    parser.add_argument(
        "--gff",
        required=True,
        help="Input GFF3 annotation file"
    )
    parser.add_argument(
        "--saf",
        required=True,
        help="Output SAF file containing merged non-redundant exon intervals"
    )
    parser.add_argument(
        "--lengths",
        required=True,
        help="Output gene effective exon length table"
    )
    parser.add_argument(
        "--gene-map",
        required=True,
        help="Output gene ID to gene name/symbol mapping table"
    )
    return parser.parse_args()


def parse_attributes(attr_string):
    """
    Parse a GFF3 attribute string into a dictionary.
    Example:
        ID=gene1;Name=ABC1;Parent=xxx
    """
    attrs = {}

    for item in attr_string.strip().split(";"):
        item = item.strip()

        if not item or "=" not in item:
            continue

        key, value = item.split("=", 1)
        attrs[key] = unquote(value)

    return attrs


def clean_id(x):
    """
    Remove common GFF3 prefixes such as:
        gene:
        transcript:
        mRNA:
    only when present at the beginning.
    """
    x = x.strip()

    prefixes = (
        "gene:",
        "transcript:",
        "mRNA:",
        "mrna:"
    )

    for prefix in prefixes:
        if x.startswith(prefix):
            return x[len(prefix):]

    return x


def merge_intervals(intervals):
    """
    Merge overlapping or directly adjacent genomic intervals.

    Input:
        [(start, end), ...]

    Output:
        [(merged_start, merged_end), ...]
    """
    if not intervals:
        return []

    intervals = sorted(intervals)

    merged = []
    cur_start, cur_end = intervals[0]

    for start, end in intervals[1:]:

        if start <= cur_end + 1:
            if end > cur_end:
                cur_end = end
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end

    merged.append((cur_start, cur_end))

    return merged


def main():

    args = parse_args()

    # ---------------------------------------------------------
    # Containers
    # ---------------------------------------------------------

    # transcript_id -> gene_id
    transcript_to_gene = {}

    # gene_id -> gene name/symbol
    gene_names = {}

    # Temporary exon records:
    # [(chrom, start, end, strand, [parent transcripts]), ...]
    exon_records = []

    # ---------------------------------------------------------
    # Parse GFF3
    # ---------------------------------------------------------

    with open(args.gff) as fh:

        for line in fh:

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])

            # -------------------------------------------------
            # Gene records
            # -------------------------------------------------

            if feature_type == "gene":

                raw_gene_id = attrs.get("ID")

                if not raw_gene_id:
                    continue

                gene_id = clean_id(raw_gene_id)

                gene_name = (
                    attrs.get("gene")
                    or attrs.get("gene_name")
                    or attrs.get("Name")
                    or attrs.get("gene_symbol")
                    or gene_id
                )

                gene_names[gene_id] = gene_name

            # -------------------------------------------------
            # Transcript / mRNA records
            # -------------------------------------------------

            elif feature_type in ("mRNA", "transcript"):

                raw_transcript_id = attrs.get("ID")
                raw_parent = attrs.get("Parent")

                if not raw_transcript_id or not raw_parent:
                    continue

                transcript_id = clean_id(raw_transcript_id)

                # Normally a transcript has one parent gene.
                # If more than one occurs, retain the first one.
                parent_gene = raw_parent.split(",")[0]
                gene_id = clean_id(parent_gene)

                transcript_to_gene[transcript_id] = gene_id

                # Some GFF files store the useful gene symbol
                # on the transcript record rather than gene record.
                if gene_id not in gene_names:

                    gene_name = (
                        attrs.get("gene")
                        or attrs.get("gene_name")
                        or attrs.get("Name")
                        or attrs.get("gene_symbol")
                    )

                    if gene_name:
                        gene_names[gene_id] = gene_name

            # -------------------------------------------------
            # Exon records
            # -------------------------------------------------

            elif feature_type == "exon":

                raw_parent = attrs.get("Parent")

                if not raw_parent:
                    continue

                parents = [
                    clean_id(x)
                    for x in raw_parent.split(",")
                    if x.strip()
                ]

                exon_records.append(
                    (
                        chrom,
                        start,
                        end,
                        strand,
                        parents
                    )
                )

    # ---------------------------------------------------------
    # Assign exons to genes through transcript -> gene mapping
    # ---------------------------------------------------------

    # key:
    # (gene_id, chrom, strand) -> [(start, end), ...]
    gene_intervals = defaultdict(list)

    missing_transcripts = set()

    for chrom, start, end, strand, parents in exon_records:

        gene_ids = set()

        for transcript_id in parents:

            gene_id = transcript_to_gene.get(transcript_id)

            if gene_id is None:
                missing_transcripts.add(transcript_id)
                continue

            gene_ids.add(gene_id)

        # A shared exon may belong to multiple transcript isoforms
        # of the same gene. Using a set prevents duplicate insertion.
        for gene_id in gene_ids:
            gene_intervals[(gene_id, chrom, strand)].append(
                (start, end)
            )

    if not gene_intervals:
        raise SystemExit(
            "ERROR: no exon intervals could be assigned to genes. "
            "Check the GFF3 ID/Parent relationships."
        )

    # ---------------------------------------------------------
    # Merge exon intervals within each gene
    # ---------------------------------------------------------

    merged_by_gene = {}

    for key, intervals in gene_intervals.items():
        merged_by_gene[key] = merge_intervals(intervals)

    # ---------------------------------------------------------
    # Write SAF
    # ---------------------------------------------------------

    with open(args.saf, "w") as out:

        out.write(
            "GeneID\tChr\tStart\tEnd\tStrand\n"
        )

        for gene_id, chrom, strand in sorted(
            merged_by_gene,
            key=lambda x: (x[1], x[0], x[2])
        ):

            for start, end in merged_by_gene[
                (gene_id, chrom, strand)
            ]:

                out.write(
                    f"{gene_id}\t"
                    f"{chrom}\t"
                    f"{start}\t"
                    f"{end}\t"
                    f"{strand}\n"
                )

    # ---------------------------------------------------------
    # Calculate effective exon length
    # ---------------------------------------------------------

    gene_lengths = defaultdict(int)

    for (gene_id, chrom, strand), intervals in merged_by_gene.items():

        for start, end in intervals:
            gene_lengths[gene_id] += end - start + 1

    # ---------------------------------------------------------
    # Write gene length table
    # ---------------------------------------------------------

    with open(args.lengths, "w") as out:

        out.write(
            "gene_id\teffective_exon_length_bp\n"
        )

        for gene_id in sorted(gene_lengths):

            length = gene_lengths[gene_id]

            # Match manuscript:
            # genes with effective exon length <= 50 bp excluded
            if length <= 50:
                continue

            out.write(
                f"{gene_id}\t{length}\n"
            )

    # ---------------------------------------------------------
    # Write gene ID -> gene name mapping
    # ---------------------------------------------------------

    with open(args.gene_map, "w") as out:

        out.write(
            "gene_id\tgene_name\n"
        )

        for gene_id in sorted(gene_lengths):

            if gene_lengths[gene_id] <= 50:
                continue

            gene_name = gene_names.get(
                gene_id,
                gene_id
            )

            out.write(
                f"{gene_id}\t{gene_name}\n"
            )

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------

    retained = sum(
        length > 50
        for length in gene_lengths.values()
    )

    n_intervals = sum(
        len(x)
        for x in merged_by_gene.values()
    )

    print(
        f"Transcripts mapped to genes: "
        f"{len(transcript_to_gene):,}"
    )

    print(
        f"Genes with exon annotation: "
        f"{len(gene_lengths):,}"
    )

    print(
        f"Genes retained (>50 bp effective exon length): "
        f"{retained:,}"
    )

    print(
        f"Merged non-redundant exon intervals: "
        f"{n_intervals:,}"
    )

    if missing_transcripts:

        print(
            f"WARNING: {len(missing_transcripts):,} "
            f"exon Parent transcript IDs were not found "
            f"in mRNA/transcript records."
        )


if __name__ == "__main__":
    main()
