import pandas as pd
import random
from collections import defaultdict
import plotly.express as px
import plotly.io as pio
import os
from .trimmomatic import run_trimmomatic
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2
import distinctipy
import numpy as np
import matplotlib.pyplot as plt

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
    """
    Generates a CSV file containing sample IDs extracted from Kraken report filenames.
    """
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

def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None, read_count=1, max_read_count=1000000000000000000000000000, rank_code='S'):
    """
    Aggregates Kraken results at a specified rank code, merging metadata or using sample IDs.

    Parameters:
    - kraken_dir (str): Path to the directory containing Kraken report files.
    - metadata_file (str, optional): Path to the metadata CSV file.
    - sample_id_df (DataFrame, optional): DataFrame of sample IDs.
    - read_count (int): Minimum read count threshold for filtering results.
    - max_read_count (int): Maximum read count threshold for filtering results.
    - rank_code (str): Taxonomic rank to aggregate ('S' for species, 'K' for kingdom, 'G' for genus, 'F' for family). Default is 'S'.
    
    Returns:
    - str: Path to the generated merged TSV file.
    """
    try:
        # Load metadata
        if metadata_file:
            metadata = pd.read_csv(metadata_file, sep=",")
            print("Using metadata from the provided metadata file.")
        elif sample_id_df is not None:
            metadata = sample_id_df
            print("Using sample IDs as metadata.")
        else:
            raise ValueError("Either metadata_file or sample_id_df must be provided.")

        sample_id_col = metadata.columns[0]

        # Define rank codes and corresponding file suffixes
        rank_map = {
            'S': ('S', ['S', 'S1', 'S2', 'S3'], ""),  # Species: no suffix
            'K': ('K', ['D', 'K'], "_K"),      # Kingdom: _K suffix
            'G': ('G', ['G'], "_G"),           # Genus: _G suffix
            'F': ('F', ['F'], "_F")            # Family: _F suffix
        }

        if rank_code not in rank_map:
            raise ValueError("Invalid rank_code. Use 'S', 'K', 'G', or 'F'.")

        selected_rank, rank_codes, file_suffix = rank_map[rank_code]

        # Dictionary to store aggregated results
        aggregated_results = {}

        # Iterate over each Kraken report file
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                with open(os.path.join(kraken_dir, file_name), 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        perc_frag_cover = fields[0]
                        nr_frag_cover = fields[1]
                        nr_frag_direct_at_taxon = int(fields[2])
                        rank_code_field = fields[3]
                        ncbi_ID = fields[4]
                        scientific_name = fields[5]
                        parts = file_name.split('_')
                        extracted_part = '_'.join(parts[:-2])
                        sampleandtaxonid = extracted_part + str(ncbi_ID)

                        # Check if rank code matches and meets read count threshold
                        if rank_code_field in rank_codes and (nr_frag_direct_at_taxon >= read_count and nr_frag_direct_at_taxon <= max_read_count):
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

        # Output aggregated results to a TSV file
        merged_tsv_path = os.path.join(kraken_dir, f"merged_kraken{file_suffix}.tsv")
        with open(merged_tsv_path, 'w') as f:
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for sampleandtaxonid, data in aggregated_results.items():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")

        print(f"Aggregated results saved to {merged_tsv_path}")
        return merged_tsv_path

    except Exception as e:
        print(f"Error aggregating Kraken results: {e}")
        return None

def generate_abundance_plots(merged_tsv_path, top_N, col_filter, pat_to_keep, rank_code='S'):
    """
    Generates abundance plots for viral and bacterial classifications at the specified rank code.

    Parameters:
    - merged_tsv_path (str): Path to the merged TSV file.
    - top_N (int): Number of top categories to plot.
    - col_filter (list): List of scientific names to exclude.
    - pat_to_keep (list): List of scientific names to include.
    - rank_code (str): Taxonomic rank for the plot ('S', 'K', 'G', 'F'). Default is 'S'.
    """
    try:
        df = pd.read_csv(merged_tsv_path, sep="\t")
        df.columns = df.columns.str.replace('/', '_').str.replace(' ', '_')
        df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
        df = df[df['Scientific_name'] != 'Homo sapiens']  # Remove human reads
        if col_filter:
            df = df[~df['Scientific_name'].isin(col_filter)]
        if pat_to_keep:
            df = df[df['Scientific_name'].isin(pat_to_keep)]

        # Define rank-specific titles and filters
        rank_titles = {
            'S': 'Species',
            'K': 'Kingdom',
            'G': 'Genus',
            'F': 'Family'
        }
        rank_suffix = '' if rank_code == 'S' else f"_{rank_code}"

        # Generate viral and bacterial abundance plots
        for focus, filter_str, plot_title in [
    (f'Virus_{rank_titles[rank_code]}', 
     'Virus' if rank_code in ['S', 'S1', 'S2'] else 'viridae' if rank_code in ['F', 'F1', 'F2'] else 'virus' if rank_code in ['G', 'G1', 'G2'] else None, 
     f'Viral_{rank_titles[rank_code]}'),
    (f'Bacteria_{rank_titles[rank_code]}', 
     'Virus' if rank_code in ['S', 'S1', 'S2'] else 'viridae' if rank_code in ['F', 'F1', 'F2'] else 'bacteria' if rank_code in ['G', 'G1', 'G2'] else None, 
     f'Bacterial_{rank_titles[rank_code]}')
]:
            if focus.startswith('Bacteria'):
                df_focus = df[~df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
            else:
                df_focus = df[df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
            df_focus = df_focus.rename(columns={'Scientific_name': focus})

            if top_N:
                top_N_categories = df_focus[focus].value_counts().head(top_N).index
                df_focus = df_focus[df_focus[focus].isin(top_N_categories)]

            categorical_cols = df_focus.select_dtypes(include=['object']).columns.tolist()
            categorical_cols.remove(focus)

            for col in categorical_cols:
                grouped_sum = df_focus.groupby([focus, col])['Nr_frag_direct_at_taxon'].mean().reset_index()
                colordict = defaultdict(int)
                random_colors1 = ['#000000', '#FF0000', '#556B2F', '#ADD8E6', '#6495ED', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000',
                                 '#008080', '#000080', '#CD5C5C', '#DAA520', '#FFA500', '#F0E68C', '#ADFF2F', '#2F4F4F', '#E0FFFF', '#4169E1', '#8A2BE2', '#4B0082', '#EE82EE', '#D2691E', '#BC8F8F', '#800080', '#DDA0DD', '#FF1493', '#8B4513', '#A0522D', '#708090', '#B0C4DE', '#FFFFF0', '#DCDCDC', '#FFEFD5', '#F5DEB3', '#7FFFD4', '#FFC0CB', '#A52A2A']
                random_colors0 = ["#{:06X}".format(random.randint(0, 0xFFFFFF)) for _ in range(len(grouped_sum[focus].unique()))]
                
                if len(grouped_sum[focus].unique()) <= len(random_colors1):
                    for target, color in zip(grouped_sum[focus].unique(), random_colors1[:len(grouped_sum[focus].unique())]):
                        colordict[target] = color
                else:
                    for target, color in zip(grouped_sum[focus].unique(), random_colors0):
                        colordict[target] = color

                plot_width = 1100 + 5 * len(grouped_sum[col].unique())
                plot_height = 800 + 5 * len(grouped_sum[col].unique())
                font_size = max(10, 14 - len(grouped_sum[col].unique()) // 10)

                fig = px.bar(
                    grouped_sum,
                    x=col,
                    y='Nr_frag_direct_at_taxon',
                    color=focus,
                    color_discrete_map=colordict,
                    title=f"{plot_title} Abundance by {col}"
                )
                summary_csv_path = os.path.join(f"{plot_title}_summary{rank_suffix}.csv")
                grouped_sum.to_csv(summary_csv_path, index=False)
                fig.update_layout(
                    xaxis=dict(tickfont=dict(size=font_size), tickangle=45),
                    yaxis=dict(tickfont=dict(size=font_size)),
                    title=dict(text=f'Average {plot_title} Abundance by {col}', x=0.5, font=dict(size=16)),
                    bargap=0.5,
                    legend=dict(
                        font=dict(size=font_size),
                        x=1,
                        y=1,
                        traceorder='normal',
                        orientation='v',
                        itemwidth=30,
                        itemsizing='constant',
                        itemclick='toggleothers',
                        itemdoubleclick='toggle'
                    ),
                    width=plot_width,
                    height=plot_height
                )

                fig.write_image(f"{plot_title}_Abundance_by_{col}{rank_suffix}.png", format='png', scale=3, width=1920, height=1080)

    except Exception as e:
        print(f"Error generating abundance plots for rank {rank_code}: {e}")

def process_all_ranks(kraken_dir, metadata_file=None, sample_id_df=None, read_count=1, max_read_count=1000000000000000000000000000, top_N=None, col_filter=None, pat_to_keep=None):
    """
    Processes Kraken results and generates abundance plots for all rank codes (S, K, G, F).
    """
    rank_codes = ['S', 'K', 'G', 'F']
    for rank in rank_codes:
        merged_tsv = aggregate_kraken_results(kraken_dir, metadata_file, sample_id_df, read_count, max_read_count, rank)
        if merged_tsv:
            generate_abundance_plots(merged_tsv, top_N, col_filter, pat_to_keep, rank)
