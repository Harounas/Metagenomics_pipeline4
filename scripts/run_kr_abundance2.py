import os
import glob
import argparse
import pandas as pd
import sys
from Metagenomics_pipeline.kraken_abundance_pipeline import process_sample, aggregate_kraken_results, generate_abundance_plots
import logging

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
        sample_id = sample_id.replace("_R1_001", "").replace("_R1", "")  # For cases without ".fastq"
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
    parser = argparse.ArgumentParser(description="Pipeline for Trimmomatic trimming, Bowtie2 host depletion (optional), and Kraken2 taxonomic classification.")
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", help="Path to Bowtie2 index (optional).")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use.")
    parser.add_argument("--metadata_file", help="Path to the metadata CSV file (optional).")
    parser.add_argument("--no_metadata", action='store_true', help="Use sample IDs as metadata instead of a metadata file.")
    parser.add_argument("--read_count", type=int, default=0, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=None, help="Select the top N most common viruses or bacteria.")
    parser.add_argument("--no_bowtie2", action='store_true', help="Skip Bowtie2 host depletion.")
    parser.add_argument("--bacteria", action='store_true', help="Generate bacterial abundance plots.")
    parser.add_argument("--virus", action='store_true', help="Generate viral abundance plots.")
    parser.add_argument("--use_precomputed_reports", action='store_true', help="Use precomputed Kraken reports instead of running Kraken2.")
    parser.add_argument("--col_filter", type=str,nargs='+', help="Bacteria or virus name to be removed")
    parser.add_argument("--pat_to_keep", type=str,nargs='+', help="Bacteria or virus name to be kept")
    parser.add_argument("--max_read_count", type=str,nargs='+', help="Maximum number of read counts")
    #parser.add_argument("--contigs_file", help="Path to a file containing paths to contig files for Kraken analysis.")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Check Kraken database existence
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)

      # Process contigs if provided
    """
    if args.contigs_file:
        print(f"Using contigs file: {args.contigs_file}")
        contig_paths = read_contig_files(args.contigs_file)
        for contig_file in contig_paths:
            if os.path.isfile(contig_file):
                # Extract base name from the file path
                cleaned_path = contig_file.rstrip("/")
                path_parts = cleaned_path.split(os.sep)
                base_name = path_parts[-2]#.replace("_denovo","")
  # Get the base name without the full path
                
                logging.info(f"Processing sample name: {base_name} for Kraken analysis.")
                logging.info(f"Processing contig file: {contig_file} for Kraken analysis.")
                
                # Directly process with Kraken without Trimmomatic or Bowtie2
                process_sample(
                    contig_file,  # Contig file as input for Kraken
                    None,  # No paired-end reads, so no reverse file
                    base_name,  # Sample base name
                    None,  # No Bowtie2 index
                    args.kraken_db,  # Kraken2 database
                    args.output_dir,  # Output directory
                    args.threads,  # Number of threads to use
                    False,  # Do not run Bowtie2
                    True,  # Process Kraken (skip trimming and Bowtie2)
                    args.use_precomputed_reports  # Use precomputed Kraken reports if specified
                )
            else:
                logging.warning(f"Contig file '{contig_file}' not found. Skipping.")
     
    else:
        # Normal processing for paired-end FASTQ files without contigs file
        run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

     """
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)

    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
            base_name = os.path.basename(forward)
            base_name = base_name.replace("_R1_001.fastq.gz", "").replace("_R1_001.fastq", "")
            base_name = base_name.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
            base_name = base_name.replace("R1.fastq.gz", "").replace("R1.fastq", "")
            base_name = base_name.replace("_R1_001", "").replace("_R1", "")

            # Define reverse file candidates
            reverse_candidates = [
                os.path.join(args.input_dir, f"{base_name}_R2_001.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2_001.fastq"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}_R2.fastq"),
                os.path.join(args.input_dir, f"{base_name}R2_001.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}R2_001.fastq"),
                os.path.join(args.input_dir, f"{base_name}R2.fastq.gz"),
                os.path.join(args.input_dir, f"{base_name}R2.fastq"),
            ]

            reverse = next((f for f in reverse_candidates if os.path.isfile(f)), None)

            if reverse:
                logging.info(f"Processing sample {base_name} with paired files.")
                process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads, run_bowtie, args.use_precomputed_reports)
            else:
                logging.warning(f"No matching R2 file found for {base_name}. Skipping.")

    # Metadata handling and abundance plot generation logic remains the same...
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        logging.info("Using sample IDs as metadata.")
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
        merged_tsv_path = aggregate_kraken_results(args.output_dir, sample_id_df=sample_id_df, read_count=args.read_count,args.max_read_count)
    else:
        if not args.metadata_file:
            raise ValueError("Metadata file must be provided if --no_metadata is not specified.")
        elif not os.path.isfile(args.metadata_file):
            logging.error(f"Metadata file '{args.metadata_file}' not found.")
            sys.exit(1)
        merged_tsv_path = aggregate_kraken_results(args.output_dir, metadata_file=args.metadata_file, read_count=args.read_count,max_read_count=args.max_read_count)

    # Generate abundance plots based on provided flags
    if merged_tsv_path and os.path.isfile(merged_tsv_path):
        if args.virus:
            logging.info("Generating viral abundance plots.")
            generate_abundance_plots(merged_tsv_path, args.top_N, args.col_filter,args.pat_to_keep)
        elif args.bacteria:
            logging.info("Generating bacterial abundance plots.")
            generate_abundance_plots(merged_tsv_path, args.top_N,args.col_filter,args.pat_to_keep)
        else:
            logging.warning("No plot type specified. Use --virus or --bacteria to generate plots.")

if __name__ == "__main__":
    main()
