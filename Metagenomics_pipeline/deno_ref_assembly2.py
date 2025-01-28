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

def run_de_novo_assembly(sample, sample_r1, sample_r2, output_dir):
    """
    Run de novo assembly using MetaSPAdes.
    """
    contigs_file = os.path.join(f"{output_dir}/{sample}_denovo", "contigs.fasta")
    if os.path.exists(contigs_file):
        logging.info(f"De novo assembly already exists for {sample}. Skipping.")
        return contigs_file

    command = f"metaspades.py -1 {sample_r1} -2 {sample_r2} -o {output_dir}/{sample}_denovo -t 32"
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

def extract_first_contig_id(fasta_file, output_file):
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                contig_id = line[1:].strip()  # Remove the '>' and any leading/trailing whitespace
                with open(output_file, 'w') as out_file:
                    out_file.write(contig_id + '\n')  # Write contig ID to the output file
                break  # Exit after writing the first contig ID




def deno_ref_based(df, input_dir, output_dir, run_bowtie):
    base_dir = "Fasta_files"
    """
    Process samples, run de novo assembly, and calculate completeness.
    """
    taxids = df['NCBI_ID'].unique()
    df['Ref_len'] = ""
    df['Consensus_len'] = ""
    df['Completeness(%)'] = ""
    df['sequence'] = ""
    dfs = []

    for tax in taxids:
        dftax = df[df['NCBI_ID'] == tax].copy()
        scientific_name = dftax['Scientific_name'].iloc[0].replace(' ', '_')
        tax_dir = os.path.join(base_dir, f"{scientific_name}_txid{tax}")
        os.makedirs(tax_dir, exist_ok=True)

        fasta_file = download_and_index_reference(tax, scientific_name, tax_dir)
        if not fasta_file:
            continue

        # Run de novo assembly for each sample
        for sample in dftax['SampleID']:

            sample_r1 = os.path.join(input_dir, f"{sample}_unmapped_1.fastq.gz" if run_bowtie else f"{sample}_trimmed_R1.fastq.gz")
            sample_r2 = os.path.join(input_dir, f"{sample}_unmapped_2.fastq.gz" if run_bowtie else f"{sample}_trimmed_R2.fastq.gz")
            
            contigs_file = run_de_novo_assembly(sample, sample_r1, sample_r2, output_dir)
            if not contigs_file:
                logging.error(f"Skipping sample {sample} due to assembly failure.")
                continue

            rag_file = os.path.join(f"{output_dir}/{sample}_{scientific_name}_rag", "ragtag.scaffold.fasta")
            # Perform further processing on the assembly
            #rag_file = os.path.join(tax_dir, "ragtag.scaffold.fasta")
            if os.path.exists(rag_file):
                logging.info(f"RagTag output already exists for {scientific_name}. Skipping.")
            else:
                command = f"ragtag.py scaffold {fasta_file} {contigs_file} -o {output_dir}/{sample}_{scientific_name}_rag"
                subprocess.run(command, shell=True, check=True)
            
            # Extract first contig
            #seqid_file = os.path.join(f"{output_dir}/{sample}_rag", "seqid.txt")
            seqid_file  = os.path.dirname(rag_file)
            logging.info(f"RagTag output path {seqid_file }.")
    
            output_file = f"Output_{sample}_{scientific_name}_seqid.txt"
            extract_first_contig_id(rag_file, output_file)
            first_contig = os.path.join(f"{seqid_file}", "first_contig.fasta")
            logging.info(f"RagTag output path {first_contig }.")
            subprocess.run(f"seqtk subseq {rag_file} {output_file}> {first_contig}", shell=True, check=True)
            subprocess.run(f"bwa index {first_contig}", shell=True, check=True)

            sample_dir = os.path.join(base_dir, f"{scientific_name}_assembled1")
            os.makedirs(sample_dir, exist_ok=True)

            bam_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_mapped_reads.bam")
            vcf_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_variants.vcf")
            consensus_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_consensus_genome.fa")
            
            bwa_command = f"bwa mem -a -t 16 {first_contig} {sample_r1} {sample_r2} | samtools view -u - | samtools sort -o {bam_file}"
            subprocess.run(bwa_command, shell=True, check=True)
            subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
            subprocess.run(f"bcftools mpileup -f {first_contig} {bam_file} | bcftools call -c --ploidy 1 -v -o {vcf_file}", shell=True, check=True)
            subprocess.run(f"samtools mpileup -aa -A -d 0 -Q 0 -f {first_contig} {bam_file} | ivar consensus -p {consensus_file}", shell=True, check=True)
            
            def calculate_length(fasta_file):
                command = f"grep -v '^>' {fasta_file} | tr -d '\n' | tr -cd 'ACTG' | wc -c"
                result = subprocess.run(command, shell=True, capture_output=True, text=True)
                return int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0
            
            try:
                ref_len = calculate_length(first_contig)
                consensus_len = calculate_length(consensus_file)
                completeness = round((consensus_len / ref_len) * 100, 2) if ref_len > 0 else 0
                sequence = SeqIO.read(consensus_file, "fasta").seq if os.path.exists(consensus_file) else ""
                logging.info(f"Consensus length {consensus_len}.")
                logging.info(f"Reference length {ref_len}.")
                logging.info(f"sequence {sequence}.")
                dftax.loc[dftax['SampleID'] == sample, ['Ref_len', 'Consensus_len', 'Completeness(%)', 'sequence']] = [ref_len, consensus_len, completeness, sequence]
            except Exception as e:
                logging.error(f"Error processing sample {sample}: {e}")
        
        dfs.append(dftax)
    
    merged_df = pd.concat(dfs, ignore_index=True)
    filtered_df = merged_df[pd.to_numeric(merged_df['Completeness(%)'], errors='coerce') >= 60]
    filtered_df.to_csv("Output-summary_complete.csv", index=False)
    merged_df.drop(columns=['sequence'], inplace=True)
