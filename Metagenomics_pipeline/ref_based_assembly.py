import os
import subprocess
import pandas as pd

def ref_based(df,run_bowtie,input_dir):
    """Run the additional processing pipeline for BWA, Samtools, BCFtools, and iVar."""
    
    # Get unique tax IDs
    taxids = df['NCBI_ID'].unique()
    os.makedirs(base_dir, exist_ok=True)  # Create the base directory if it doesn't exist
    base_dir="Fasta_files"
    # Iterate over tax IDs
    for tax in taxids:
        # Filter DataFrame for the current tax ID
        dftax = df[df['NCBI_ID'] == tax]

        # Get the first scientific name for the tax ID (assumes consistency)
        scientific_name = dftax['Scientific_name'].iloc[0].replace(' ', '_')  # Replace spaces with underscores
        tax_dir = os.path.join(base_dir, f"{scientific_name}_txid{tax}")
        os.makedirs(tax_dir, exist_ok=True)

        # File path for the FASTA file
        fasta_file = os.path.join(tax_dir, f"{scientific_name}.fasta")
        Refname = scientific_name  # Reference name
        Ref = fasta_file  # Reference file path

        # Construct the EDirect command
        command = f'esearch -db nucleotide -query "txid{tax}[Organism]" | efilter -source refseq | efetch -format fasta > {fasta_file}'
        print(f"Running command: {command}")  # For debugging

        # Run the EDirect command using subprocess
        subprocess.run(command, shell=True, check=True)

        # Check if the FASTA file was created successfully
        if os.path.exists(fasta_file):
            print(f"FASTA file saved to: {fasta_file}")
        else:
            print(f"FASTA file {fasta_file} was not created.")
            continue  # Skip this tax ID if the FASTA file is missing

        # Get the list of sample IDs for the current tax ID
        samplelist = dftax['SampleID'].tolist()
        print(f"Sample list for tax ID {tax} ({scientific_name}): {samplelist}")

        # BWA and other commands for each sample
        for sample in samplelist:
       
            if run_bowtie:
               sample_r1 = os.path.join(input_dir, f"{sample}_unmapped_1.fastq.gz")  # Path to R1 FASTQ file
               sample_r2 = os.path.join(input_dir, f"{sample}_unmapped_2.fastq.gz")  # Path to R2 FASTQ file
            else:
               sample_r1 = os.path.join(input_dir, f"{sample}_trimmed_R1.fastq.gz")  # Path to trimmed R1 FASTQ file
               sample_r2 = os.path.join(input_dir, f"{sample}_trimmed_R2.fastq.gz")  # Path to trimmed R2 FASTQ file    
            
            
            print(sample_r1)
            print(f'{input_dir} input dir') 
            sample_dir = f"{sample}_assembled1"
            os.makedirs(sample_dir, exist_ok=True)
            bam_file = os.path.join(sample_dir, f"{sample}_{Refname}_mapped_reads.bam")
            vcf_file = os.path.join(sample_dir, f"{sample}_{Refname}_variants.vcf")
            consensus_file = os.path.join(sample_dir, f"{sample}_{Refname}_consensus_genome.fasta")

            # BWA MEM: Align reads to the reference genome
            bwa_command = f"bwa mem -a -t 16 {Ref} {sample_r1} {sample_r2} | samtools view -u -@ 3 - | samtools sort -@ 16 -o {bam_file}"
            print(f"Running BWA command: {bwa_command}")
            subprocess.run(bwa_command, shell=True, check=True)

            # Index the BAM file
            samtools_index_command = f"samtools index {bam_file}"
            print(f"Running Samtools Index command: {samtools_index_command}")
            subprocess.run(samtools_index_command, shell=True, check=True)

            # Generate VCF file
            bcftools_command = f"bcftools mpileup -f {Ref} {bam_file} | bcftools call -c --ploidy 1 -v -o {vcf_file}"
            print(f"Running BCFtools command: {bcftools_command}")
            subprocess.run(bcftools_command, shell=True, check=True)

            # Generate consensus genome using iVar
            ivar_command = f"samtools mpileup  -aa -A -d 0 -Q 0 -f {Ref} {bam_file} | ivar consensus -p {consensus_file}"
            print(f"Running iVar command: {ivar_command}")
            subprocess.run(ivar_command, shell=True, check=True)
