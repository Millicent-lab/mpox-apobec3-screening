import pandas as pd
from Bio import SeqIO
import re

# Load reference sequence from FASTA file
def load_reference_sequence(fasta_path):
    ref_seq = ""
    for record in SeqIO.parse(fasta_path, "fasta"):
        ref_seq += str(record.seq)
    return ref_seq

# Extract context ±10 nucleotides around a position
def extract_context(ref_seq, position, window=10):
    start = max(0, position - window - 1)
    end = min(len(ref_seq), position + window)
    return ref_seq[start:position - 1], ref_seq[position:end]

# Parse GFF file to extract gene annotations
def parse_gff(gff_path):
    gene_annotations = []
    with open(gff_path, "r") as gff:
        for line in gff:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                if len(columns) >= 9:
                    start = int(columns[3])
                    end = int(columns[4])
                    strand = columns[6]
                    gene_info = re.search(r"Name=([^;]+)", columns[8])
                    gene_name = gene_info.group(1) if gene_info else "Unknown"
                    gene_annotations.append({
                        "Start": start,
                        "End": end,
                        "Strand": strand,
                        "Gene_Name": gene_name
                    })
    return pd.DataFrame(gene_annotations)

# Map positions from the VCF file to genes and strands using the GFF annotations
def map_positions_to_genes(variant_data, gene_data):
    gene_list = []
    for _, row in variant_data.iterrows():
        position = row["POS"]
        gene_match = gene_data[(gene_data["Start"] <= position) & (gene_data["End"] >= position)]
        if not gene_match.empty:
            gene_list.append({
                "Gene_Name": gene_match.iloc[0]["Gene_Name"],
                "Strand": gene_match.iloc[0]["Strand"]
            })
        else:
            gene_list.append({"Gene_Name": "Intergenic", "Strand": "Unknown"})
    return pd.DataFrame(gene_list)

# Analyze mutations and classify as APOBEC or non-APOBEC
def analyze_variants(variant_data, ref_seq, gene_data):
    results = []
    gene_list = map_positions_to_genes(variant_data, gene_data)
    for i, row in variant_data.iterrows():
        position = row["POS"]
        ref_allele = row["REF"]
        alt_allele = row["ALT"]
        variant_type = row["Type"]
        aa_change = row["AA_Change"]
        gene_name = gene_list.iloc[i]["Gene_Name"]
        strand = gene_list.iloc[i]["Strand"]

        upstream, downstream = extract_context(ref_seq, position)
        context = upstream + ref_allele + downstream
        mutant_context = upstream + alt_allele + downstream

        # Check for known APOBEC motifs
        matched_motif = "None"
        motif_match = False

        if ref_allele == "C" and alt_allele == "T":
            if re.search(r"TC", context):
                motif_match = True
                matched_motif = "TC→TT"
        elif ref_allele == "G" and alt_allele == "A":
            if re.search(r"GA", context):
                motif_match = True
                matched_motif = "GA→AA"
            elif re.search(r"GG", context):
                motif_match = True
                matched_motif = "GG→AG"

        apobec_nucleotide_match = (
            (ref_allele == "C" and alt_allele == "T") or
            (ref_allele == "G" and alt_allele == "A")
        )

        is_apobec = motif_match and apobec_nucleotide_match and strand in ["+", "-"]

        results.append({
            "Position": position,
            "Mutation": f"{ref_allele}->{alt_allele}",
            "Context": context,
            "Mutant_Context": mutant_context,
            "Strand": strand,
            "Gene_Name": gene_name,
            "Motif_Match": motif_match,
            "Matched_Motif": matched_motif,
            "APOBEC_Nucleotide_Match": apobec_nucleotide_match,
            "APOBEC": "Yes" if is_apobec else "No",
            "Expected_Motif": "TC, GA, GG",
            "Expected_AA_Change": "Depends on mutation type",
            "Type": variant_type,
            "AA_Change": aa_change
        })

    return pd.DataFrame(results)

# Process VCF file into a DataFrame
def process_vcf(vcf_path):
    variants = []
    with open(vcf_path, "r") as vcf:
        for line in vcf:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                pos = int(columns[1])
                ref = columns[3]
                alt = columns[4]
                info = columns[7]

                type_match = re.search(r"ANN=.*?\|([^\|]*)\|", info)
                aa_match = re.search(r"ANN=.*?\|.*?\|.*?\|.*?\|.*?\|([^\|]*)\|", info)
                variant_type = type_match.group(1) if type_match else "Unknown"
                aa_change = aa_match.group(1) if aa_match else "None"

                variants.append({
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "Type": variant_type,
                    "AA_Change": aa_change
                })
    return pd.DataFrame(variants)

# Main workflow
fasta_path = "C:\\Users\\millicent\\Downloads\\SnEff\\SnEff_Results\\sequences.fa"
vcf_path = "C:\\Users\\millicent\\Downloads\\Mutation_Analysis\\SnEff\\SnEff_Results\\representativeseqfiles\\Ia\\annotated_variants_hgvs.vcf"
gff_path = "C:\\Users\\millicent\\Downloads\\Mutation_Analysis\\SnEff\\data\\trimmed_refseq\\genes.gff"
output_path = "C:\\Users\\millicent\\Downloads\\apobec_analysis_results.csv"

print("Loading reference sequence...")
ref_sequence = load_reference_sequence(fasta_path)

print("Parsing GFF file...")
gene_data = parse_gff(gff_path)

print("Processing VCF file...")
variant_data = process_vcf(vcf_path)

print("Analyzing variants for APOBEC patterns and annotating...")
results = analyze_variants(variant_data, ref_sequence, gene_data)

# Summary statistics
synonymous = results[results["Type"] == "synonymous_variant"]
nonsynonymous = results[results["Type"] == "missense_variant"]
apobec = results[results["APOBEC"] == "Yes"]
non_apobec = results[results["APOBEC"] == "No"]

print("=== Summary ===")
print(f"Total variants: {len(results)}")
print(f"Synonymous: {len(synonymous)}, Non-synonymous: {len(nonsynonymous)}")
print(f"APOBEC mutations: {len(apobec)}, Non-APOBEC mutations: {len(non_apobec)}")

# Detailed breakdown
syn_apobec = synonymous[synonymous["APOBEC"] == "Yes"]
syn_non_apobec = synonymous[synonymous["APOBEC"] == "No"]
nonsyn_apobec = nonsynonymous[nonsynonymous["APOBEC"] == "Yes"]
nonsyn_non_apobec = nonsynonymous[nonsynonymous["APOBEC"] == "No"]

print("\n=== Detailed Breakdown ===")
print(f"Synonymous - APOBEC: {len(syn_apobec)}, Non-APOBEC: {len(syn_non_apobec)}")
print(f"Non-synonymous - APOBEC: {len(nonsyn_apobec)}, Non-APOBEC: {len(nonsyn_non_apobec)}")

# New: Classification by variant type (all types)
print("\n=== APOBEC Classification by Variant Type ===")
variant_type_summary = results.groupby(["Type", "APOBEC"]).size().unstack(fill_value=0)
print(variant_type_summary)

# Save results
print("Saving results to Downloads folder...")
results.to_csv(output_path, index=False)

# Save the new summary
summary_output_path = output_path.replace(".csv", "_summary_by_type.csv")
variant_type_summary.to_csv(summary_output_path)

print(f"Results saved to '{output_path}'")
print(f"Summary of APOBEC classification by type saved to '{summary_output_path}'")
