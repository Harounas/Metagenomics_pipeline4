import os
import glob
import argparse
import pandas as pd
import sys
import subprocess
import logging
from collections import defaultdict
from Metagenomics_pipeline.kraken_abundance_pipeline import (
    process_sample,
    generate_abundance_plots,
    process_all_ranks,
    generate_unfiltered_merged_tsv,
    run_multiqc,
    aggregate_kraken_results,
    process_kraken_reports,
    extract_domains_from_kraken_report,
    run_fastqc
)
from Metagenomics_pipeline.ref_based_assembly import ref_based
from Metagenomics_pipeline.deno_ref_assembly2 import deno_ref_based

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def create_sample_id_df(input_dir):
    """
    Create a DataFrame with sample IDs based on the input FASTQ file names.
    """
    sample_ids = []
    for f in glob.glob(os.path.join(input_dir, "*_R1*.fastq*")):
        sample_id = os.path.basename(f).split('_R1')[0]
        sample_ids.append(sample_id)
    return pd.DataFrame(sample_ids, columns=["Sample_IDs"])

def read_contig_files(contig_file):
    """
    Reads a file containing paths to contig fasta files and returns a list of file paths.
    """
    try:
        with open(contig_file, 'r') as f:
            return [line.strip() for line in f.readlines() if line.strip()]
    except FileNotFoundError:
        logging.error(f"Contig file '{contig_file}' not found.")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Pipeline for FastQC, trimming, host depletion, and Kraken2 taxonomic classification."
    )
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", help="Path to Bowtie2 index (optional).")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use.")
    parser.add_argument("--metadata_file", help="Path to the metadata CSV file (optional).")
    parser.add_argument("--no_metadata", action='store_true', help="Use sample IDs as metadata instead of a metadata file.")
    parser.add_argument("--read_count", type=int, default=1, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=10**13, help="Select the top N most common taxa.")
    parser.add_argument("--no_bowtie2", action='store_true', help="Skip Bowtie2 host depletion.")
    parser.add_argument("--use_precomputed_reports", action='store_true', help="Use precomputed Kraken reports instead of running Kraken2.")
    parser.add_argument("--bacteria", action='store_true', help="Process Bacteria domain.")
    parser.add_argument("--virus", action='store_true', help="Process Viruses domain.")
    parser.add_argument("--archaea", action='store_true', help="Process Archaea domain.")
    parser.add_argument("--eukaryota", action='store_true', help="Process Eukaryota domain.")
    parser.add_argument("--run_ref_base", action="store_true", help="Run reference-based assembly pipeline.")
    parser.add_argument("--run_deno_ref", action="store_true", help="Run de novo reference assembly pipeline.")
    parser.add_argument("--process_all_ranks", action='store_true', help="Process all taxonomic ranks.")
    parser.add_argument("--skip_fastqc", action='store_true', help="Skip FastQC quality control.")
    parser.add_argument("--skip_multiqc", action='store_true', help="Skip MultiQC report generation.")
    parser.add_argument("--skip_reports", action='store_true', help="Skip processing Kraken reports and generating plots; directly run genome assembly.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)

    if not args.skip_fastqc:
        run_fastqc(args.output_dir, args.threads)
    else:
        logging.info("Skipping FastQC per user request.")

    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward).split('_R1')[0]
        reverse = next((f for f in [
            os.path.join(args.input_dir, f"{base_name}_R2_001.fastq.gz"),
            os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz"),
            os.path.join(args.input_dir, f"{base_name}_R2.fastq")
        ] if os.path.isfile(f)), None)

        if reverse:
            logging.info(f"Processing sample {base_name} with paired files.")
            process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db,
                           args.output_dir, args.threads, run_bowtie, args.use_precomputed_reports)
        else:
            logging.warning(f"No matching R2 file found for {base_name}. Skipping.")

    if not args.skip_multiqc:
        run_multiqc(args.output_dir)
    else:
        logging.info("Skipping MultiQC per user request.")

    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)

    if not args.skip_reports:
        process_kraken_reports(args.output_dir)

if __name__ == "__main__":
    main()
