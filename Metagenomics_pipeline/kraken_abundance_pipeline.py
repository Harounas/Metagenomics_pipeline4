import pandas as pd
import random
from collections import defaultdict
import plotly.express as px
import plotly.io as pio
import os
import logging
from .trimmomatic import run_trimmomatic
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2
import distinctipy
import numpy as np
import matplotlib.pyplot as plt
import subprocess

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads, run_bowtie, use_precomputed_reports):
    try:
        if not use_precomputed_reports:
            # Step 1: Run Trimmomatic
            trimmed_forward, trimmed_reverse = run_trimmomatic(forward, reverse, base_name, output_dir, threads)
            
            # Step 2: Optionally run Bowtie2 to deplete host genome reads
            if run_bowtie:
                unmapped_r1, unmapped_r2 = run_bowtie2(trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads)
            else:
                unmapped_r1, unmapped_r2 = trimmed_forward, trimmed_reverse
             
            # Step 3: Run Kraken2 with the reads
            kraken_report = run_kraken2(unmapped_r1, unmapped_r2, base_name, kraken_db, output_dir, threads)
        else:
            # Use the precomputed Kraken2 report
            kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
            if not os.path.exists(kraken_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report not found: {kraken_report}")
        
        return kraken_report
    
    except Exception as e:
        print(f"Error processing sample {base_name}: {e}")
        return None

def generate_sample_ids_csv(kraken_dir):
    """Generates a CSV file containing sample IDs extracted from Kraken report filenames."""
    try:
        sample_ids = []
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                sample_id = '_'.join(file_name.split('_')[:-2])
                sample_ids.append(sample_id)
        
        sampleid_df = pd.DataFrame(sample_ids, columns=['Sample_IDs'])
        sampleid_csv_path = os.path.join(kraken_dir, "sample_ids.csv")
        sampleid_df.to_csv(sampleid_csv_path, index=False)
        
        print(f"Sample IDs saved to {sampleid_csv_path}")
        return sampleid_csv_path
    
    except Exception as e:
        print(f"Error generating sample_ids.csv: {e}")
        return None


def run_fastqc(output_dir, threads):
    """
    Runs FastQC on all FASTQ files in the output directory.
    """
    try:
        fastq_files = glob.glob(os.path.join(output_dir, "*trimmed*.fastq*"))
        if not fastq_files:
            logging.warning("No FASTQ files found in the specified directory.")
            return
        for fastq_file in fastq_files:
            logging.info(f"Running FastQC on {fastq_file}")
            subprocess.run(
                ["fastqc", fastq_file, "-o", output_dir, "-t", str(threads)],
                check=True
            )
        logging.info("FastQC completed successfully.")
    except Exception as e:
        logging.error(f"Error running FastQC: {e}")

def extract_domains_from_kraken_report(kraken_report_path):
    """
    Extracts rows for each domain (Bacteria, Eukaryota, Archaea, Viruses) from a Kraken2 report.
    Returns a dictionary with keys as domain names and values as DataFrames.
    """
    columns = [
        "Percentage", "Reads_Covered", "Reads_Assigned", "Rank_Code",
        "NCBI_TaxID", "Scientific_Name"
    ]
    df = pd.read_csv(kraken_report_path, sep="\t", header=None, names=columns)
    domains = {}
    current_domain = None
    current_rows = []
    for _, row in df.iterrows():
        if row["Rank_Code"] == "D":
            if current_domain:
                domains[current_domain] = pd.DataFrame(current_rows)
            current_domain = row["Scientific_Name"]
            current_rows = [row]
        else:
            current_rows.append(row)
    if current_domain:
        domains[current_domain] = pd.DataFrame(current_rows)
    return domains

def process_kraken_reports(kraken_dir):
    """
    Processes Kraken2 report files in kraken_dir by splitting each into domain-specific files.
    The output filenames are of the form: {sample}_{DomainWithoutSpaces}_kraken_report.txt.
    """
    for file_name in os.listdir(kraken_dir):
        if file_name.endswith("_report.txt"):
            parts = file_name.split('_')
            parts = [item for item in parts if item not in ['Viruses', 'Eukaryota', 'Bacteria', 'Archaea']]
            extracted_part = '_'.join(parts[:-2])
            kraken_report_path = os.path.join(kraken_dir, file_name)
            domains = extract_domains_from_kraken_report(kraken_report_path)
            for domain, df in domains.items():
                output_path = os.path.join(
                    kraken_dir, f"{extracted_part}_{domain.replace(' ', '')}_kraken_report.txt"
                )
                df.to_csv(output_path, sep="\t", index=False,header=False)
                logging.info(f"Saved {domain} data to {output_path}")

def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None,
                             min_read_counts=None, max_read_counts=None, rank_code='S', domain_filter=None):
    """
    Aggregates Kraken results at a specified rank code, applying per-domain read count filtering.
    
    For example, at species level (default 'S'), rows with Rank_code in ['S', 'S1', 'S2', 'S3']
    are selected; at Family level ('F'), rows with Rank_code in ['F', 'F1', 'F2', 'F3'] are selected,
    and so on.
    
    Args:
        - min_read_counts (dict): Dictionary of minimum read count per domain.
        - max_read_counts (dict): Dictionary of maximum read count per domain.
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
        aggregated_results = {}
        
        # Define rank mapping: keys are desired rank level; values are acceptable rank codes.
        rank_mapping = {
            'S': ['S', 'S1', 'S2', 'S3'],
            'K': ['K', 'K1', 'K2', 'K3'],
            'F': ['F', 'F1', 'F2', 'F3'],
            'D': ['D', 'D1', 'D2', 'D3']
        }
        
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                if domain_filter and domain_filter.replace(' ', '').lower() not in file_name.lower():
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
                        extracted_part = '_'.join(parts[:-3])
                        sampleandtaxonid = extracted_part + str(ncbi_ID)

                        # Determine domain from file name using keys in min_read_counts
                        domain = None
                        for key in min_read_counts.keys():
                            if key.lower() in file_name.lower():
                                domain = key
                                break

                        # Apply domain-specific min/max read count
                        min_read = min_read_counts.get(domain, 1)
                        max_read = max_read_counts.get(domain, 10**30)

                        # Filter rows using the rank mapping
                        if rank_code_field in rank_mapping.get(rank_code, [rank_code]) and \
                           (min_read <= nr_frag_direct_at_taxon <= max_read):
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
        
        domain_suffix = f"_{domain_filter.replace(' ', '')}" if domain_filter else "_all"
        merged_tsv_path = os.path.join(kraken_dir, f"merged_kraken_{rank_code}{domain_suffix}.tsv")
        
        with open(merged_tsv_path, 'w') as f:
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 
                       'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
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

def generate_abundance_plots(merged_tsv_path, top_N, col_filter, pat_to_keep, rank_code='S'):
    """
    Generates abundance plots for classifications at the specified rank code,
    considering Viruses, Eukaryota, Bacteria, and Archaea.
    
    For example, at species level (default 'S'), rows with Rank_code in ['S','S1','S2','S3']
    are selected; similarly for other ranks.
    """
    try:
        # Get the directory where the aggregated TSV files are stored.
        kraken_dir = os.path.dirname(merged_tsv_path)
        
        # Define human-friendly rank titles and suffix for file naming.
        rank_titles = {
            'S': 'Species',
            'K': 'Kingdom',
            'G': 'Genus',
            'F': 'Family',
            'D': 'Domain'
        }
        rank_suffix = '' if rank_code == 'S' else f"_{rank_code}"
        
        # Mapping for rank-level filtering.
        rank_mapping = {
            'S': ['S', 'S1', 'S2', 'S3'],
            'K': ['K', 'K1', 'K2', 'K3'],
            'G': ['G', 'G1', 'G2', 'G3'],
            'F': ['F', 'F1', 'F2', 'F3'],
            'D': ['D', 'D1', 'D2', 'D3']
        }
        
        # Define the domain-specific TSV files based on the current rank.
        categories = [
            ('Viruses', os.path.join(kraken_dir, f"merged_kraken_{rank_code}_Viruses.tsv"), 'Viral'),
            ('Bacteria', os.path.join(kraken_dir, f"merged_kraken_{rank_code}_Bacteria.tsv"), 'Bacterial'),
            ('Archaea', os.path.join(kraken_dir, f"merged_kraken_{rank_code}_Archaea.tsv"), 'Archaeal'),
            ('Eukaryota', os.path.join(kraken_dir, f"merged_kraken_{rank_code}_Eukaryota.tsv"), 'Eukaryotic')
        ]
        
        for focus, file_name, plot_title in categories:
            try:
                # Load the domain-specific aggregated TSV file.
                df_focus = pd.read_csv(file_name, sep="\t")
                
                # Filter rows by Rank_code using the rank mapping.
                df_focus = df_focus[df_focus['Rank_code'].isin(rank_mapping.get(rank_code, [rank_code]))]
                
                # Rename the 'Scientific_name' column to the current domain name.
                df_focus = df_focus.rename(columns={'Scientific_name': focus})
                
                # Optionally restrict to top N categories.
                if top_N:
                    top_N_categories = df_focus[focus].value_counts().head(top_N).index
                    df_focus = df_focus[df_focus[focus].isin(top_N_categories)]
                
                # Optionally apply additional filtering.
                if col_filter:
                    df_focus = df_focus[~df_focus[focus].isin(col_filter)]
                if pat_to_keep:
                    df_focus = df_focus[df_focus[focus].isin(pat_to_keep)]
                
                # Identify categorical columns (excluding the domain column).
                categorical_cols = df_focus.select_dtypes(include=['object']).columns.tolist()
                if focus in categorical_cols:
                    categorical_cols.remove(focus)
                
                # Use groupby with as_index=False to avoid duplicate column insertion.
                for col in categorical_cols:
                    grouped_sum = df_focus.groupby([focus, col], as_index=False)['Nr_frag_direct_at_taxon'].mean()
                    
                    fig = px.bar(
                        grouped_sum,
                        x=col,
                        y='Nr_frag_direct_at_taxon',
                        color=focus,
                        title=f"{plot_title} {rank_titles[rank_code]} Abundance by {col}"
                    )
                    output_img = f"{plot_title}_Abundance_by_{col}{rank_suffix}.png"
                    fig.write_image(output_img, format='png', scale=3, width=1920, height=1080)
                    logging.info(f"Abundance plot saved to {output_img}")
            except FileNotFoundError:
                logging.warning(f"File {file_name} not found, skipping {focus} plot generation.")
    
    except Exception as e:
        logging.error(f"Error generating abundance plots for rank {rank_code}: {e}")


def process_all_ranks(kraken_dir, metadata_file=None, sample_id_df=None, read_count=1,
                      max_read_count=10**30, top_N=None, col_filter=None, pat_to_keep=None):
    """
    Processes Kraken results, generates abundance plots for all rank codes (S, K, G, F),
    and creates an unfiltered merged TSV.
    """
    unfiltered_tsv = generate_unfiltered_merged_tsv(kraken_dir, metadata_file, sample_id_df)
    rank_codes = ['S', 'K', 'G', 'F']
    for rank in rank_codes:
        merged_tsv = aggregate_kraken_results(kraken_dir, metadata_file, sample_id_df,
                                              read_count, max_read_count, rank)
        if merged_tsv:
            generate_abundance_plots(merged_tsv, top_N, col_filter, pat_to_keep, rank)
    return unfiltered_tsv

def run_multiqc(trimmomatic_output_dir):
    """Runs MultiQC on Trimmomatic output files."""
    try:
        subprocess.run(["multiqc", trimmomatic_output_dir], check=True)
        logging.info("MultiQC report generated successfully.")
    except Exception as e:
        logging.error(f"Error running MultiQC: {e}")
