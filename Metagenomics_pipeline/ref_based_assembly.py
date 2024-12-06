import os
import subprocess
from Bio import SeqIO
import pandas as pd

def ref_based(df, run_bowtie, input_dir):
    """
    Ref-based pipeline with consensus genome polishing and denovo assembly.

    Parameters:
    df (pd.DataFrame): DataFrame containing sample metadata (columns include 'NCBI_ID', 'Scientific_name', 'SampleID').
    run_bowtie (bool): Whether to process unmapped reads (True) or trimmed reads (False).
    input_dir (str): Path to input directory containing reads.

    Returns:
    pd.DataFrame: Updated DataFrame with genome statistics.
    """
    # Initialize output columns
    df['Ref_len'] = ""
    df['Consensus_len'] = ""
    df['Completeness(%)'] = ""

    dfs = []  # List to store updated DataFrame chunks

    for tax in df['NCBI_ID'].unique():
        dftax = df[df['NCBI_ID'] == tax]
        scientific_name = dftax['Scientific_name'].iloc[0].replace(' ', '_')
        tax_dir = os.path.join("Fasta_files", f"{scientific_name}_txid{tax}")
        os.makedirs(tax_dir, exist_ok=True)

        fasta_file = os.path.join(tax_dir, f"{scientific_name}.fasta")
        command = f'esearch -db nucleotide -query "txid{tax}[Organism]" | efilter -source refseq | efetch -format fasta > {fasta_file}'

        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running EDirect: {e}")
            continue

        if not os.path.exists(fasta_file):
            print(f"FASTA file not created for {tax}")
            continue

        bwa_ref = f"bwa index {fasta_file}"
        try:
            subprocess.run(bwa_ref, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error indexing reference: {e}")
            continue

        for sample in dftax['SampleID']:
            # Determine input files based on run_bowtie flag
            sample_r1 = os.path.join(input_dir, f"{sample}_unmapped_1.fastq.gz") if run_bowtie else os.path.join(input_dir, f"{sample}_trimmed_R1.fastq.gz")
            sample_r2 = os.path.join(input_dir, f"{sample}_unmapped_2.fastq.gz") if run_bowtie else os.path.join(input_dir, f"{sample}_trimmed_R2.fastq.gz")

            sample_dir = os.path.join(tax_dir, f"{scientific_name}_assembled")
            os.makedirs(sample_dir, exist_ok=True)

            bam_file = os.path.join(sample_dir, f"{sample}_mapped_reads.bam")
            consensus_file = os.path.join(sample_dir, f"{sample}_consensus_genome.fa")

            # Denovo assembly using metaspades
            denovo_output_dir = os.path.join(input_dir,f"{sample}_denovo")
            denovo_command = f"metaspades.py -1 {sample_r1} -2 {sample_r2} -o {denovo_output_dir} -t 32"
            try:
                print(f"Running denovo assembly: {denovo_command}")
                subprocess.run(denovo_command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error in denovo assembly for {sample}: {e}")
                continue

            denovo_contigs = os.path.join(denovo_output_dir, "contigs.fasta")
            if not os.path.exists(denovo_contigs):
                print(f"Denovo assembly failed: contigs not generated for {sample}")
                continue

            # BWA MEM: Align reads to the reference genome
            bwa_command = f"bwa mem -t 16 {fasta_file} {sample_r1} {sample_r2} | samtools sort -o {bam_file}"
            try:
                subprocess.run(bwa_command, shell=True, check=True)
                subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error processing sample {sample}: {e}")
                continue

            # Generate initial consensus genome using iVar
            ivar_command = f"samtools mpileup -aa -A -d 0 -Q 0 -f {fasta_file} {bam_file} | ivar consensus -p {consensus_file}"
            try:
                subprocess.run(ivar_command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running iVar for {sample}: {e}")
                continue

            # Polishing step using denovo contigs
            vcf_contig_consensus_file = os.path.join(sample_dir, f"{sample}_consensus_variants.vcf")
            consensus_contig_file = consensus_file
            consensus_contig_polished_file = os.path.join(sample_dir, f"{sample}_consensus_polished_genome.fa")
            bam_contig_consensus_file = os.path.join(sample_dir, f"{sample}_consensus_reads.bam")

            # Index the consensus genome
            bwa_contig = f"bwa index {consensus_contig_file}"
            try:
                subprocess.run(bwa_contig, shell=True, check=True)

                # Align the consensus genome against the denovo contigs
                bwa_command_consensus = f"bwa mem -a -t 16 {consensus_contig_file} {denovo_contigs} | samtools sort -o {bam_contig_consensus_file}"
                subprocess.run(bwa_command_consensus, shell=True, check=True)

                # Index the BAM file
                samtools_index_consensus_command = f"samtools index {bam_contig_consensus_file}"
                subprocess.run(samtools_index_consensus_command, shell=True, check=True)

                # Generate VCF file for consensus variants
                bcftools_command_consensus = f"bcftools mpileup -f {consensus_contig_file} {bam_contig_consensus_file} | bcftools call -c --ploidy 1 -v -o {vcf_contig_consensus_file}"
                subprocess.run(bcftools_command_consensus, shell=True, check=True)

                # Generate polished consensus genome using iVar
                ivar_command_polished = f"samtools mpileup -aa -A -d 0 -Q 0 -f {consensus_contig_file} {bam_contig_consensus_file} | ivar consensus -p {consensus_contig_polished_file}"
                subprocess.run(ivar_command_polished, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error during polishing for sample {sample}: {e}")
                continue

            # Calculate genome statistics and update DataFrame
            try:
                ref_len = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))
                consensus_len = sum(len(record.seq) for record in SeqIO.parse(consensus_contig_polished_file, "fasta"))
                completeness = (consensus_len / ref_len) * 100 if ref_len > 0 else 0
                print(f"Sample {sample}: Ref_len={ref_len}, Consensus_len={consensus_len}, Completeness={completeness:.2f}%")

                dftax.loc[dftax['SampleID'] == sample, 'Ref_len'] = ref_len
                dftax.loc[dftax['SampleID'] == sample, 'Consensus_len'] = consensus_len
                dftax.loc[dftax['SampleID'] == sample, 'Completeness(%)'] = f"{completeness:.2f}"
            except Exception as e:
                print(f"Error calculating statistics for {sample}: {e}")

        dfs.append(dftax)

    merged_df=pd.concat(dfs).reset_index(drop=True)
    merged_df.to_csv("all-summary_len.csv", index=False)
