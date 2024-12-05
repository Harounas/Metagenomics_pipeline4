def main():
    import os
    import glob
    import argparse
    import pandas as pd
    import sys
    from Metagenomics_pipeline.kraken_abundance_pipeline import process_sample, aggregate_kraken_results, generate_abundance_plots
    from Metagenomics_pipeline.ref_based_assembly import ref_based
    import logging

    # Configure logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    def create_sample_id_df(input_dir):
        sample_ids = [
            os.path.basename(f).split("_R1")[0]
            for f in glob.glob(os.path.join(input_dir, "*_R1*.fastq*"))
        ]
        return pd.DataFrame(sample_ids, columns=["Sample_IDs"])

    def read_contig_files(contig_file):
        try:
            with open(contig_file, 'r') as f:
                return [line.strip() for line in f if line.strip()]
        except FileNotFoundError:
            logging.error(f"Contig file '{contig_file}' not found.")
            sys.exit(1)

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Pipeline for Kraken2 analysis.")
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", help="Path to Bowtie2 index (optional).")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use.")
    parser.add_argument("--metadata_file", help="Path to the metadata CSV file (optional).")
    parser.add_argument("--no_metadata", action='store_true', help="Use sample IDs as metadata.")
    parser.add_argument("--read_count", type=int, default=1, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=10000, help="Select top N most common taxa.")
    parser.add_argument("--no_bowtie2", action='store_true', help="Skip Bowtie2 host depletion.")
    parser.add_argument("--bacteria", action='store_true', help="Generate bacterial abundance plots.")
    parser.add_argument("--virus", action='store_true', help="Generate viral abundance plots.")
    parser.add_argument("--use_precomputed_reports", action='store_true', help="Use precomputed Kraken reports.")
    parser.add_argument("--col_filter", type=str, nargs='+', help="List of taxa to exclude.")
    parser.add_argument("--pat_to_keep", type=str, nargs='+', help="List of taxa to include.")
    parser.add_argument("--max_read_count", type=int, default=5000000000, help="Maximum read count.")
    parser.add_argument("--run_ref_base", action="store_true", help="Run reference-based pipeline.")
    parser.add_argument("--fastq_path", type=str, help="Explicit path to FASTQ files for reference-based analysis.")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Check Kraken database existence
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)

    # Determine Bowtie2 usage
    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

    # Process samples
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1*.fastq*")):
        base_name = os.path.basename(forward).split("_R1")[0]
        reverse_candidates = [os.path.join(args.input_dir, f"{base_name}_R2{ext}")
                              for ext in ["_001.fastq.gz", "_001.fastq", ".fastq.gz", ".fastq"]]
        reverse = next((f for f in reverse_candidates if os.path.isfile(f)), None)

        if reverse:
            logging.info(f"Processing paired files for sample {base_name}.")
            process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db,
                           args.output_dir, args.threads, run_bowtie, args.use_precomputed_reports)
        else:
            logging.warning(f"Missing R2 file for sample {base_name}. Skipping.")

    # Aggregate Kraken results
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
        merged_tsv_path = aggregate_kraken_results(args.output_dir, sample_id_df=sample_id_df,
                                                   read_count=args.read_count, max_read_count=args.max_read_count)
    elif args.metadata_file and os.path.isfile(args.metadata_file):
        merged_tsv_path = aggregate_kraken_results(args.output_dir, metadata_file=args.metadata_file,
                                                   read_count=args.read_count, max_read_count=args.max_read_count)
    else:
        logging.error("Metadata file is missing. Exiting.")
        sys.exit(1)

    # Generate plots
    if os.path.isfile(merged_tsv_path):
        if args.virus:
            logging.info("Generating viral abundance plots.")
            generate_abundance_plots(merged_tsv_path, args.top_N, args.col_filter, args.pat_to_keep)
        if args.bacteria:
            logging.info("Generating bacterial abundance plots.")
            generate_abundance_plots(merged_tsv_path, args.top_N, args.col_filter, args.pat_to_keep)

    # Run reference-based pipeline
    # Run reference-based pipeline
    if args.run_ref_base:
        if args.fastq_path:
            logging.info(f"Using FASTQ files from explicitly provided path: {args.fastq_path}")
            forward = os.path.join(args.fastq_path, "*_R1*.fastq*")
            reverse = os.path.join(args.fastq_path, "*_R2*.fastq*")
        else:
            logging.info("Using FASTQ files from the input directory.")
            forward = os.path.join(args.input_dir, "*_R1*.fastq*")
            reverse = os.path.join(args.input_dir, "*_R2*.fastq*")

        # Verify if files exist
        forward_files = glob.glob(forward)
        reverse_files = glob.glob(reverse)
        if not forward_files or not reverse_files:
            logging.error("No FASTQ files found for reference-based pipeline. Exiting.")
            sys.exit(1)

        logging.info("Running reference-based pipeline.")
        df = pd.read_csv(merged_tsv_path, sep='\t')
        df = df[df['Scientific_name'].str.contains('virus', case=False, na=False)]
        df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
        ref_based(df, run_bowtie, args.output_dir, forward_files, reverse_files)

    if __name__ == "__main__":
        main()
