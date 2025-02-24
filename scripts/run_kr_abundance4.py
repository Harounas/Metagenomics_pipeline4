import os
import glob
import argparse
import pandas as pd
import sys
from Metagenomics_pipeline.kraken_abundance_pipeline import process_sample, aggregate_kraken_results, generate_abundance_plots, process_all_ranks
from Metagenomics_pipeline.ref_based_assembly import ref_based
from Metagenomics_pipeline.deno_ref_assembly2 import deno_ref_based
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def create_sample_id_df(input_dir):
    """
    Create a DataFrame with sample IDs based on the input FASTQ file names.
    
    Args:
        input_dir (str): Directory containing input FASTQ files.
    
    Returns:
        pd.DataFrame: DataFrame with a single column "Sample_IDs" containing sample identifiers.
    """
    sample_ids = []
    for f in glob.glob(os.path.join(input_dir, "*_R1*.fastq*")):
        sample_id = os.path.basename(f)
        sample_id = sample_id.replace("_R1_001.fastq.gz", "").replace("_R1_001.fastq", "")
        sample_id = sample_id.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        sample_id = sample_id.replace("R1.fastq.gz", "").replace("R1.fastq", "")
        sample_id = sample_id.replace("_R1_001", "").replace("_R1", "")  # For cases without ".fastq"
        sample_ids.append(sample_id)

    sample_id_df = pd.DataFrame(sample_ids, columns=["Sample_IDs"])
    return sample_id_df

def read_contig_files(contig_file):
    """
    Reads a file containing paths to contig fasta files and returns a list of file paths.
    
    Args:
        contig_file (str): Path to the file listing contig fasta file paths.
    
    Returns:
        list: List of contig file paths.
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
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Pipeline for Trimmomatic trimming, Bowtie2 host depletion (optional), and Kraken2 taxonomic classification.")
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", help="Path to Bowtie2 index (optional).")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use.")
    parser.add_argument("--metadata_file", help="Path to the metadata CSV file (optional).")
    parser.add_argument("--no_metadata", action='store_true', help="Use sample IDs as metadata instead of a metadata file.")
    parser.add_argument("--read_count", type=int, default=1, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=10000, help="Select the top N most common viruses or bacteria.")
    parser.add_argument("--no_bowtie2", action='store_true', help="Skip Bowtie2 host depletion.")
    parser.add_argument("--bacteria", action='store_true', help="Generate bacterial abundance plots.")
    parser.add_argument("--virus", action='store_true', help="Generate viral abundance plots.")
    parser.add_argument("--use_precomputed_reports", action='store_true', help="Use precomputed Kraken reports instead of running Kraken2.")
    parser.add_argument("--col_filter", type=str, nargs='+', help="Bacteria or virus name to be removed")
    parser.add_argument("--pat_to_keep", type=str, nargs='+', help="Bacteria or virus name to be kept")
    parser.add_argument("--max_read_count", type=int, default=5000000000, help="Maximum number of read counts")
    parser.add_argument("--run_ref_base", action="store_true", help="Run the additional processing pipeline for each taxon (BWA, Samtools, BCFtools, iVar)")
    parser.add_argument("--run_deno_ref", action="store_true", help="Run the additional processing pipeline for each taxon (BWA, Samtools, BCFtools, iVar)")
    parser.add_argument("--process_all_ranks", action='store_true', help="Process and generate abundance plots for all taxonomic ranks (S, K, G, F).")
    parser.add_argument("--filtered_tsv", help="Path to the filtered merged Kraken output .tsv file for assembly (optional).")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Validate Kraken database
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)

    # Determine if Bowtie2 should run
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

    # Process each sample
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward).replace("_R1_001.fastq.gz", "").replace("_R1_001.fastq", "")
        base_name = base_name.replace("_R1.fastq.gz", "").replace("_R1.fastq", "").replace("R1.fastq.gz", "").replace("R1.fastq", "").replace("_R1_001", "").replace("_R1", "")

        reverse_candidates = [
            os.path.join(args.input_dir, f"{base_name}_R2_001.fastq.gz"),
            os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz"),
            os.path.join(args.input_dir, f"{base_name}_R2.fastq"),
        ]

        reverse = next((f for f in reverse_candidates if os.path.isfile(f)), None)

        if reverse:
            logging.info(f"Processing sample {base_name} with paired files.")
            process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads, run_bowtie, args.use_precomputed_reports)
        else:
            logging.warning(f"No matching R2 file found for {base_name}. Skipping.")

    # Initialize merged_tsv_path
    merged_tsv_path = None

    # Handle aggregation or use provided filtered TSV
    if not args.process_all_ranks:
        if args.filtered_tsv:
            # Validate provided filtered TSV
            if not os.path.isfile(args.filtered_tsv):
                logging.error(f"Provided filtered_tsv file '{args.filtered_tsv}' not found.")
                sys.exit(1)
            merged_tsv_path = args.filtered_tsv
            logging.info(f"Using provided filtered merged Kraken output: {merged_tsv_path}")
        else:
            # Aggregate Kraken results
            if args.no_metadata:
                sample_id_df = create_sample_id_df(args.input_dir)
                logging.info("Using sample IDs as metadata.")
                sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
                merged_tsv_path = aggregate_kraken_results(args.output_dir, sample_id_df=sample_id_df, read_count=args.read_count, max_read_count=args.max_read_count)
            else:
                if not args.metadata_file or not os.path.isfile(args.metadata_file):
                    logging.error(f"Metadata file '{args.metadata_file}' not found.")
                    sys.exit(1)
                merged_tsv_path = aggregate_kraken_results(args.output_dir, metadata_file=args.metadata_file, read_count=args.read_count, max_read_count=args.max_read_count)
    else:
        # Warn if filtered_tsv is provided with process_all_ranks
        if args.filtered_tsv:
            logging.warning("--filtered_tsv is provided but will be ignored since --process_all_ranks is set.")
        if args.no_metadata:
            sample_id_df = create_sample_id_df(args.input_dir)
            logging.info("Using sample IDs as metadata.")
            sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
            process_all_ranks(args.output_dir, sample_id_df=sample_id_df, read_count=args.read_count, max_read_count=args.max_read_count, top_N=args.top_N, col_filter=args.col_filter, pat_to_keep=args.pat_to_keep)
        else:
            if not args.metadata_file or not os.path.isfile(args.metadata_file):
                logging.error(f"Metadata file '{args.metadata_file}' not found.")
                sys.exit(1)
            process_all_ranks(args.output_dir, metadata_file=args.metadata_file, read_count=args.read_count, max_read_count=args.max_read_count, top_N=args.top_N, col_filter=args.col_filter, pat_to_keep=args.pat_to_keep)

    # Generate abundance plots for species level
    if not args.process_all_ranks and merged_tsv_path and os.path.isfile(merged_tsv_path):
        if args.virus:
            logging.info("Generating viral abundance plots for species level.")
            generate_abundance_plots(merged_tsv_path, args.top_N, args.col_filter, args.pat_to_keep)
        if args.bacteria:
            logging.info("Generating bacterial abundance plots for species level.")
            generate_abundance_plots(merged_tsv_path, args.top_N, args.col_filter, args.pat_to_keep)

    # Run additional reference-based processing for species level
    if not args.process_all_ranks and merged_tsv_path and os.path.isfile(merged_tsv_path):
        df = pd.read_csv(merged_tsv_path, sep='\t')
        df = df[df['Scientific_name'].str.contains('virus', case=False, na=False)]
        df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))

        if args.run_ref_base:
            logging.info("Starting reference-based pipeline.")
            ref_based(df, run_bowtie, args.output_dir)
        if args.run_deno_ref:
            logging.info("Starting de novo reference assembly pipeline.")
            deno_ref_based(df, args.output_dir, args.output_dir, run_bowtie)

if __name__ == "__main__":
    main()
