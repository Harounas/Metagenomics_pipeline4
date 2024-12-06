import os
import subprocess
from Bio import SeqIO
import pandas as pd
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pipeline.log"),
        logging.StreamHandler()
    ]
)

def download_and_index_reference(tax, scientific_name, tax_dir):
    """
    Download and index the reference genome for a given taxonomic ID.

    Parameters:
    tax (str): Taxonomic ID.
    scientific_name (str): Scientific name of the organism.
    tax_dir (str): Directory to store reference files.

    Returns:
    str: Path to the indexed FASTA file, or None if failed.
    """
    fasta_file = os.path.join(tax_dir, f"{scientific_name}.fasta")
    command = f'esearch -db nucleotide -query "txid{tax}[Organism]" | efilter -source refseq | efetch -format fasta > {fasta_file}'
    try:
        logging.info(f"Downloading reference genome for {scientific_name} (taxid {tax})...")
        subprocess.run(command, shell=True, check=True)
        subprocess.run(f"bwa index {fasta_file}", shell=True, check=True)
        logging.info(f"Reference genome downloaded and indexed: {fasta_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error preparing reference for {tax}: {e}")
        return None
    return fasta_file

def run_denovo_assembly(sample, sample_r1, sample_r2, output_dir):
    """
    Run de novo assembly using MetaSPAdes.

    Parameters:
    sample (str): Sample ID.
    sample_r1 (str): Path to forward reads.
    sample_r2 (str): Path to reverse reads.
    output_dir (str): Directory to store de novo assembly results.

    Returns:
    str: Path to the contigs file, or None if failed.
    """
    contigs_file = os.path.join(output_dir, "contigs.fasta")
    if os.path.exists(contigs_file):
        logging.info(f"De novo assembly already exists for {sample}. Skipping.")
        return contigs_file

    command = f"metaspades.py -1 {sample_r1} -2 {sample_r2} -o {output_dir} -t 32"
    try:
        logging.info(f"Running de novo assembly for {sample}: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in de novo assembly for {sample}: {e}")
        return None

    if not os.path.exists(contigs_file):
        logging.error(f"De novo assembly failed: contigs not generated for {sample}")
        return None

    logging.info(f"De novo assembly completed for {sample}: {contigs_file}")
    return contigs_file

def align_reads_to_reference(fasta_file, sample_r1, sample_r2, output_dir, sample):
    """
    Align reads to a reference genome using BWA-MEM.

    Parameters:
    fasta_file (str): Path to the reference genome.
    sample_r1 (str): Path to forward reads.
    sample_r2 (str): Path to reverse reads.
    output_dir (str): Directory to store BAM files.
    sample (str): Sample ID.

    Returns:
    str: Path to the sorted BAM file, or None if failed.
    """
    bam_file = os.path.join(output_dir, f"{sample}_mapped_reads.bam")
    command = f"bwa mem -t 16 {fasta_file} {sample_r1} {sample_r2} | samtools sort -o {bam_file}"
    try:
        logging.info(f"Aligning reads for {sample}: {command}")
        subprocess.run(command, shell=True, check=True)
        subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error aligning reads for {sample}: {e}")
        return None

    if not os.path.exists(bam_file):
        logging.error(f"BAM file not generated for {sample}")
        return None

    logging.info(f"Alignment completed for {sample}: {bam_file}")
    return bam_file

def calculate_completeness(fasta_file, consensus_file):
    """
    Calculate genome completeness based on reference and consensus lengths.

    Parameters:
    fasta_file (str): Path to the reference genome.
    consensus_file (str): Path to the consensus genome.

    Returns:
    tuple: Reference length, consensus length, and completeness percentage.
    """
    try:
        ref_len = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))
        consensus_len = sum(len(record.seq) for record in SeqIO.parse(consensus_file, "fasta"))
        completeness = (consensus_len / ref_len) * 100 if ref_len > 0 else 0
        return ref_len, consensus_len, completeness
    except Exception as e:
        logging.error(f"Error calculating completeness: {e}")
        return 0, 0, 0

def ref_based(df, run_bowtie, input_dir):
    """
    Ref-based pipeline with consensus genome polishing and de novo assembly.

    Parameters:
    df (pd.DataFrame): DataFrame containing sample metadata.
    run_bowtie (bool): Whether to process unmapped reads (True) or trimmed reads (False).
    input_dir (str): Path to input directory containing reads.

    Returns:
    pd.DataFrame: Updated DataFrame with genome statistics.
    """
    df['Ref_len'] = ""
    df['Consensus_len'] = ""
    df['Completeness(%)'] = ""

    dfs = []
    for tax in df['NCBI_ID'].unique():
        dftax = df[df['NCBI_ID'] == tax]
        scientific_name = dftax['Scientific_name'].iloc[0].replace(' ', '_')
        tax_dir = os.path.join("Fasta_files", f"{scientific_name}_txid{tax}")
        os.makedirs(tax_dir, exist_ok=True)

        fasta_file = download_and_index_reference(tax, scientific_name, tax_dir)
        if not fasta_file:
            continue

        for sample in dftax['SampleID']:
            sample_r1 = os.path.join(input_dir, f"{sample}_unmapped_1.fastq.gz") if run_bowtie else os.path.join(input_dir, f"{sample}_trimmed_R1.fastq.gz")
            sample_r2 = os.path.join(input_dir, f"{sample}_unmapped_2.fastq.gz") if run_bowtie else os.path.join(input_dir, f"{sample}_trimmed_R2.fastq.gz")

            sample_dir = os.path.join(tax_dir, f"{scientific_name}_assembled")
            os.makedirs(sample_dir, exist_ok=True)

            bam_file = align_reads_to_reference(fasta_file, sample_r1, sample_r2, sample_dir, sample)
            if not bam_file:
                continue

            denovo_output_dir = os.path.join(input_dir, f"{sample}_denovo")
            denovo_contigs = run_denovo_assembly(sample, sample_r1, sample_r2, denovo_output_dir)
            if not denovo_contigs:
                continue

            # Further steps: Polishing, consensus genome creation, etc.
            consensus_file = os.path.join(sample_dir, f"{sample}_consensus_genome.fa")
            ref_len, consensus_len, completeness = calculate_completeness(fasta_file, consensus_file)
            logging.info(f"Sample {sample}: Ref_len={ref_len}, Consensus_len={consensus_len}, Completeness={completeness:.2f}%")

            dftax.loc[dftax['SampleID'] == sample, 'Ref_len'] = ref_len
            dftax.loc[dftax['SampleID'] == sample, 'Consensus_len'] = consensus_len
            dftax.loc[dftax['SampleID'] == sample, 'Completeness(%)'] = f"{completeness:.2f}"

        dfs.append(dftax)

    merged_df = pd.concat(dfs).reset_index(drop=True)
    
    merged_df.to_csv(f"{input_dir}/all-summary_len.csv", index=False)
    return merged_df
