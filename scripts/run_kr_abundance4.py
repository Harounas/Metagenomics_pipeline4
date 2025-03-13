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
        sample_id = os.path.basename(f)
        sample_id = sample_id.replace("_R1_001.fastq.gz", "").replace("_R1_001.fastq", "")
        sample_id = sample_id.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        sample_id = sample_id.replace("R1.fastq.gz", "").replace("R1.fastq", "")
        sample_id = sample_id.replace("_R1_001", "").replace("_R1", "")
        sample_ids.append(sample_id)
    sample_id_df = pd.DataFrame(sample_ids, columns=["Sample_IDs"])
    return sample_id_df

def read_contig_files(contig_file):
    """
    Reads a file containing paths to contig fasta files and returns a list of file paths.
    """
    contig_paths = []
    try:
        with open(contig_file, 'r') as f:
            contig_paths = [line.strip() for line in f.readlines() if line.strip()]
    except FileNotFoundError:
        logging.error(f"Contig file '{contig_file}' not found.")
        sys.exit(1)
    return contig_paths

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
    # Domain-specific flags for abundance plotting and assembly:
    parser.add_argument("--bacteria", action='store_true', help="Process Bacteria domain.")
    parser.add_argument("--virus", action='store_true', help="Process Viruses domain.")
    parser.add_argument("--archaea", action='store_true', help="Process Archaea domain.")
    parser.add_argument("--eukaryota", action='store_true', help="Process Eukaryota domain.")
    parser.add_argument("--run_ref_base", action="store_true", help="Run reference-based assembly pipeline.")
    parser.add_argument("--run_deno_ref", action="store_true", help="Run de novo reference assembly pipeline.")
    parser.add_argument("--process_all_ranks", action='store_true', help="Process all taxonomic ranks.")
    parser.add_argument("--filtered_tsv", help="Path to the filtered merged Kraken output TSV for assembly (optional).")
    parser.add_argument("--skip_fastqc", action='store_true', help="Skip FastQC quality control.")
    parser.add_argument("--skip_multiqc", action='store_true', help="Skip MultiQC report generation.")
    parser.add_argument("--max_read_count", type=int, default=10**13, help="Maximum number of read counts")
    parser.add_argument("--min_read_bacteria", type=int, default=1, help="Minimum read count for Bacteria.")
    parser.add_argument("--max_read_bacteria", type=int, default=10**13, help="Maximum read count for Bacteria.")
    parser.add_argument("--min_read_virus", type=int, default=1, help="Minimum read count for Viruses.")
    parser.add_argument("--max_read_virus", type=int, default=10**13, help="Maximum read count for Viruses.")
    parser.add_argument("--min_read_archaea", type=int, default=1, help="Minimum read count for Archaea.")
    parser.add_argument("--max_read_archaea", type=int, default=10**13, help="Maximum read count for Archaea.")
    parser.add_argument("--min_read_eukaryota", type=int, default=1, help="Minimum read count for Eukaryota.")
    parser.add_argument("--max_read_eukaryota", type=int, default=10**13, help="Maximum read count for Eukaryota.")
    # New arguments for filtering:
    parser.add_argument("--col_filter", nargs='+', help="List of taxa to filter out.")
    parser.add_argument("--pat_to_keep", nargs='+', help="List of taxa to retain.")
    # New flag: if set, skip Kraken report processing and abundance plot generation,
    # and directly generate the unfiltered TSV to pass to genome assembly.
    parser.add_argument("--skip_reports", action='store_true', help="Skip processing Kraken reports and generating plots; directly run genome assembly.")
    args = parser.parse_args()

    # Define per-domain read count thresholds
    min_read_counts = {
        "Bacteria": args.min_read_bacteria,
        "Viruses": args.min_read_virus,
        "Archaea": args.min_read_archaea,
        "Eukaryota": args.min_read_eukaryota
    }
    max_read_counts = {
        "Bacteria": args.max_read_bacteria,
        "Viruses": args.max_read_virus,
        "Archaea": args.max_read_archaea,
        "Eukaryota": args.max_read_eukaryota
    }

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)

    # Run FastQC conditionally
    if not args.skip_fastqc:
        run_fastqc(args.output_dir, args.threads)
    else:
        logging.info("Skipping FastQC per user request.")
    
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward)
        base_name = base_name.replace("_R1_001.fastq.gz", "").replace("_R1_001.fastq", "")
        base_name = base_name.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        base_name = base_name.replace("R1.fastq.gz", "").replace("R1.fastq", "")
        base_name = base_name.replace("_R1_001", "").replace("_R1", "")
        reverse_candidates = [
            os.path.join(args.input_dir, f"{base_name}_R2_001.fastq.gz"),
            os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz"),
            os.path.join(args.input_dir, f"{base_name}_R2.fastq"),
        ]
        reverse = next((f for f in reverse_candidates if os.path.isfile(f)), None)
        if reverse:
            logging.info(f"Processing sample {base_name} with paired files.")
            process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db,
                           args.output_dir, args.threads, run_bowtie, args.use_precomputed_reports)
        else:
            logging.warning(f"No matching R2 file found for {base_name}. Skipping.")
    
    # Run MultiQC conditionally
    if not args.skip_multiqc:
        run_multiqc(args.output_dir)
    else:
        logging.info("Skipping MultiQC per user request.")
    
    # Generate sample IDs CSV (if needed)
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
    else:
        sample_id_df = None

    if not args.skip_reports:
        # Process domain-specific Kraken reports
        process_kraken_reports(args.output_dir)
        
        # For each domain flag provided, aggregate results and run abundance plots and assembly pipelines.
        domains_to_process = []
        if args.bacteria:
            domains_to_process.append("Bacteria")
        if args.virus:
            domains_to_process.append("Viruses")
        if args.archaea:
            domains_to_process.append("Archaea")
        if args.eukaryota:
            domains_to_process.append("Eukaryota")
        
        # Define the rank codes to generate for each domain.
        # For Bacteria, Archaea, and Eukaryota: species, family, domain, and kingdom.
        # For Viruses: species, family, and domain.
        domain_rank_codes = {
            "Bacteria": ['S', 'F', 'D', 'K'],
            "Viruses": ['S', 'F', 'D'],
            "Archaea": ['S', 'F', 'D', 'K'],
            "Eukaryota": ['S', 'F', 'D', 'K']
        }
        
        for domain in domains_to_process:
            logging.info(f"Aggregating results for domain: {domain}")
            for rank in domain_rank_codes.get(domain, ['S']):
                merged_tsv = aggregate_kraken_results(
                    args.output_dir,
                    args.metadata_file,
                    sample_id_df,
                    min_read_counts=min_read_counts,
                    max_read_counts=max_read_counts,
                    rank_code=rank,
                    domain_filter=domain
                )
                if merged_tsv and os.path.isfile(merged_tsv):
                    logging.info(f"Generating abundance plots for {domain} at rank {rank}.")
                    generate_abundance_plots(merged_tsv, args.top_N, args.col_filter, args.pat_to_keep, rank)
                    if args.run_ref_base:
                        df = pd.read_csv(merged_tsv, sep="\t")
                        df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                        df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                        logging.info(f"Starting reference-based assembly for {domain} at rank {rank}.")
                        ref_based(df, run_bowtie, args.output_dir)
                    if args.run_deno_ref:
                        df = pd.read_csv(merged_tsv, sep="\t")
                        df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                        df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                        logging.info(f"Starting de novo reference assembly for {domain} at rank {rank}.")
                        deno_ref_based(df, args.output_dir, args.output_dir, run_bowtie)
    else:
        logging.info("Skipping domain-specific report processing and plot generation. Running genome assembly directly using unfiltered merged TSV.")
        # When skipping report processing, generate an unfiltered TSV and use it for assembly.
        merged_tsv = generate_unfiltered_merged_tsv(args.output_dir, args.metadata_file, sample_id_df)
        if merged_tsv and os.path.isfile(merged_tsv):
            if args.run_ref_base:
                df = pd.read_csv(merged_tsv, sep="\t")
                df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                logging.info("Starting reference-based assembly for all genomes using unfiltered TSV.")
                ref_based(df, run_bowtie, args.output_dir)
            if args.run_deno_ref:
                df = pd.read_csv(merged_tsv, sep="\t")
                df = df[~df['Scientific_name'].str.contains('Homo sapiens', case=False, na=False)]
                df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                logging.info("Starting de novo reference assembly for all genomes using unfiltered TSV.")
                deno_ref_based(df, args.output_dir, args.output_dir, run_bowtie)
    
    # Optionally process all ranks (if --process_all_ranks is set)
    if args.process_all_ranks:
        if args.no_metadata:
            sample_id_df = create_sample_id_df(args.input_dir)
            sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
            process_all_ranks(args.output_dir, sample_id_df=sample_id_df,
                              read_count=args.read_count, max_read_count=args.top_N,
                              top_N=args.top_N, col_filter=args.col_filter, pat_to_keep=args.pat_to_keep)
        else:
            if not args.metadata_file or not os.path.isfile(args.metadata_file):
                logging.error(f"Metadata file '{args.metadata_file}' not found.")
                sys.exit(1)
            process_all_ranks(args.output_dir, metadata_file=args.metadata_file,
                              read_count=args.read_count, max_read_count=args.top_N,
                              top_N=args.top_N, col_filter=args.col_filter, pat_to_keep=args.pat_to_keep)

if __name__ == "__main__":
    main()
