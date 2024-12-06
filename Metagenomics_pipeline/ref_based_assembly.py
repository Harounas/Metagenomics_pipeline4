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

def calculate_genome_length(fasta_file):
    """
    Calculate the genome length based on valid ACTG nucleotides in a FASTA file.
    """
    command = f"grep -v '^>' {fasta_file} | tr -d '\\n' | tr -cd 'ACTG' | wc -c"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())

def ref_based(df, run_bowtie, input_dir):
    """
    Run a reference-based pipeline with consensus genome creation.
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

            consensus_file = os.path.join(sample_dir, f"{sample}_consensus_genome.fa")
            try:
                ref_len = calculate_genome_length(fasta_file)
                consensus_len = calculate_genome_length(consensus_file)
                completeness = (consensus_len / ref_len) * 100 if ref_len > 0 else 0

                dftax.loc[dftax['SampleID'] == sample, 'Ref_len'] = ref_len
                dftax.loc[dftax['SampleID'] == sample, 'Consensus_len'] = consensus_len
                dftax.loc[dftax['SampleID'] == sample, 'Completeness(%)'] = f"{completeness:.2f}"
            except Exception as e:
                logging.error(f"Error calculating lengths for {sample}: {e}")

        dfs.append(dftax)

    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df.to_csv(f"{input_dir}/summary_genome_stats.csv", index=False)
    logging.info(f"Pipeline completed. Results saved to: {input_dir}/summary_genome_stats.csv")
