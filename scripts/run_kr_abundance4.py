# Initialize merged_tsv_path to None
merged_tsv_path = None

if not args.process_all_ranks:
    if args.filtered_tsv:
        # Check if the provided file exists
        if not os.path.isfile(args.filtered_tsv):
            logging.error(f"Provided filtered_tsv file '{args.filtered_tsv}' not found.")
            sys.exit(1)
        merged_tsv_path = args.filtered_tsv
        logging.info(f"Using provided filtered merged Kraken output: {merged_tsv_path}")
    else:
        # Proceed with aggregation if --filtered_tsv is not provided
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
    # If --process_all_ranks is set, assembly is skipped, so warn if --filtered_tsv is provided
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

# Generate abundance plots for species level if not processing all ranks and merged_tsv_path is available
if not args.process_all_ranks and merged_tsv_path and os.path.isfile(merged_tsv_path):
    if args.virus:
        logging.info("Generating viral abundance plots for species level.")
        generate_abundance_plots(merged_tsv_path, args.top_N, args.col_filter, args.pat_to_keep)
    if args.bacteria:
        logging.info("Generating bacterial abundance plots for species level.")
        generate_abundance_plots(merged_tsv_path, args.top_N, args.col_filter, args.pat_to_keep)

# Additional reference-based processing (only for species level)
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
