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
    run_multiqc,aggregate_kraken_results,process_kraken_reports,extract_domains_from_kraken_report,run_fastqc
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





def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None,
                             read_count=1, max_read_count=10**30, rank_code='S',
                             domain_filter=None):
    """
    Aggregates Kraken results at a specified rank code.
    If domain_filter is provided (e.g. "Bacteria"), only files whose names contain that domain
    (with spaces removed) will be processed. The output merged TSV filename reflects the domain.
    """
    try:
        # Load metadata
        if metadata_file:
            metadata = pd.read_csv(metadata_file, sep=",")
            logging.info("Using metadata from the provided metadata file.")
        elif sample_id_df is not None:
            metadata = sample_id_df
            logging.info("Using sample IDs as metadata.")
        else:
            raise ValueError("Either metadata_file or sample_id_df must be provided.")
        
        sample_id_col = metadata.columns[0]
        
        # Define rank codes and file suffixes (adjust if needed)
        rank_map = {
            'S': ('S', ['S', 'S1', 'S2', 'S3'], ""),
            'K': ('K', ['D', 'K'], "_K"),
            'G': ('G', ['G'], "_G"),
            'F': ('F', ['F'], "_F")
        }
        
        if rank_code not in rank_map:
            raise ValueError("Invalid rank_code. Use 'S', 'K', 'G', or 'F'.")
        
        selected_rank, rank_codes, file_suffix = rank_map[rank_code]
        
        aggregated_results = {}
        
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                # If a domain_filter is provided, check that the filename contains it
                if domain_filter and domain_filter.replace(' ', '') not in file_name:
                    continue
                
                with open(os.path.join(kraken_dir, file_name), 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        if len(fields) < 6:
                            continue
                        perc_frag_cover = fields[0]
                        nr_frag_cover = fields[1]
                        nr_frag_direct_at_taxon = int(fields[2])
                        rank_code_field = fields[3]
                        ncbi_ID = fields[4]
                        scientific_name = fields[5]
                        parts = file_name.split('_')
                        extracted_part = '_'.join(parts[:-2])
                        sampleandtaxonid = extracted_part + str(ncbi_ID)
                        
                        if rank_code_field in rank_codes and (read_count <= nr_frag_direct_at_taxon <= max_read_count):
                            if extracted_part in metadata[sample_id_col].unique():
                                sample_metadata = metadata.loc[metadata[sample_id_col] == extracted_part].iloc[0].to_dict()
                                aggregated_results[sampleandtaxonid] = {
                                    'Perc_frag_cover': perc_frag_cover,
                                    'Nr_frag_cover': nr_frag_cover,
                                    'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                                    'Rank_code': rank_code_field,
                                    'NCBI_ID': ncbi_ID,
                                    'Scientific_name': scientific_name,
                                    'SampleID': extracted_part,
                                    **sample_metadata
                                }
        
        domain_suffix = f"_{domain_filter.replace(' ', '')}" if domain_filter else ""
        merged_tsv_path = os.path.join(kraken_dir, f"merged_kraken{file_suffix}{domain_suffix}.tsv")
        with open(merged_tsv_path, 'w') as f:
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for data in aggregated_results.values():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")
        
        logging.info(f"Aggregated results saved to {merged_tsv_path}")
        return merged_tsv_path
    
    except Exception as e:
        logging.error(f"Error aggregating Kraken results: {e}")
        return None

def generate_unfiltered_merged_tsv(kraken_dir, metadata_file=None, sample_id_df=None):
    """
    Generates a merged TSV file containing Kraken report data for all ranks without filtering.
    """
    try:
        if metadata_file:
            metadata = pd.read_csv(metadata_file, sep=",")
            logging.info("Using metadata from the provided metadata file.")
        elif sample_id_df is not None:
            metadata = sample_id_df
            logging.info("Using sample IDs as metadata.")
        else:
            raise ValueError("Either metadata_file or sample_id_df must be provided.")
        
        sample_id_col = metadata.columns[0]
        unfiltered_results = {}
        
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                with open(os.path.join(kraken_dir, file_name), 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        if len(fields) < 6:
                            continue
                        perc_frag_cover = fields[0]
                        nr_frag_cover = fields[1]
                        nr_frag_direct_at_taxon = int(fields[2])
                        rank_code_field = fields[3]
                        ncbi_ID = fields[4]
                        scientific_name = fields[5]
                        parts = file_name.split('_')
                        extracted_part = '_'.join(parts[:-2])
                        sampleandtaxonid = extracted_part + str(ncbi_ID)
                        
                        if extracted_part in metadata[sample_id_col].unique():
                            sample_metadata = metadata.loc[metadata[sample_id_col] == extracted_part].iloc[0].to_dict()
                            unfiltered_results[sampleandtaxonid] = {
                                'Perc_frag_cover': perc_frag_cover,
                                'Nr_frag_cover': nr_frag_cover,
                                'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                                'Rank_code': rank_code_field,
                                'NCBI_ID': ncbi_ID,
                                'Scientific_name': scientific_name,
                                'SampleID': extracted_part,
                                **sample_metadata
                            }
        
        merged_tsv_path = os.path.join(kraken_dir, "merged_kraken_all_ranks_unfiltered.tsv")
        with open(merged_tsv_path, 'w') as f:
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for data in unfiltered_results.values():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")
        
        logging.info(f"Unfiltered merged Kraken results saved to {merged_tsv_path}")
        return merged_tsv_path
    
    except Exception as e:
        logging.error(f"Error generating unfiltered merged TSV: {e}")
        return None




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
    parser.add_argument("--top_N", type=int, default=10000, help="Select the top N most common taxa.")
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
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)

    run_fastqc(args.output_dir, args.threads)
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
                           args.output_dir, args.threads, run_bowtie,  args.use_precomputed_reports)
        else:
            logging.warning(f"No matching R2 file found for {base_name}. Skipping.")

    run_multiqc(args.output_dir)
    
    # Generate sample IDs CSV (if needed)
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
    else:
        sample_id_df = None

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
    
    for domain in domains_to_process:
        logging.info(f"Aggregating results for domain: {domain}")
        merged_tsv = aggregate_kraken_results(args.output_dir, args.metadata_file, sample_id_df,
                                               args.read_count, args.top_N, 'S', domain_filter=domain)
        if merged_tsv and os.path.isfile(merged_tsv):
            logging.info(f"Generating abundance plots for {domain}.")
            generate_abundance_plots(merged_tsv, args.top_N, None, None, 'S')
            if args.run_ref_base:
                df = pd.read_csv(merged_tsv, sep="\t")
                df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                logging.info(f"Starting reference-based assembly for {domain}.")
                ref_based(df, run_bowtie, args.output_dir)
            if args.run_deno_ref:
                df = pd.read_csv(merged_tsv, sep="\t")
                df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                logging.info(f"Starting de novo reference assembly for {domain}.")
                deno_ref_based(df, args.output_dir, args.output_dir, run_bowtie)
    
    # Optionally process all ranks (if --process_all_ranks is set)
    if args.process_all_ranks:
        if args.no_metadata:
            sample_id_df = create_sample_id_df(args.input_dir)
            sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
            process_all_ranks(args.output_dir, sample_id_df=sample_id_df,
                              read_count=args.read_count, max_read_count=args.top_N,
                              top_N=args.top_N)
        else:
            if not args.metadata_file or not os.path.isfile(args.metadata_file):
                logging.error(f"Metadata file '{args.metadata_file}' not found.")
                sys.exit(1)
            process_all_ranks(args.output_dir, metadata_file=args.metadata_file,
                              read_count=args.read_count, max_read_count=args.top_N,
                              top_N=args.top_N)

if __name__ == "__main__":
    main()
