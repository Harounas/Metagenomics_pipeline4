import os
import subprocess
import pandas as pd

def ref_based(df, run_bowtie, input_dir):
    """
    Run the additional processing pipeline for BWA, Samtools, BCFtools, and iVar.

    Parameters:
    - df: pandas.DataFrame containing information about samples and tax IDs.
    - run_bowtie: bool, whether to use Bowtie-derived FASTQ files or trimmed FASTQ files.
    - input_dir: str, the directory containing the input FASTQ files.
    """
    
    # Define base directory for output FASTA files
    base_dir = "Fasta_files"
    os.makedirs(base_dir, exist_ok=True)  # Create the base directory if it doesn't exist

    # Get unique tax IDs
    taxids = df['NCBI_ID'].unique()

    # Iterate over tax IDs
    df['Ref_len'] =""
    df['Consensus_len']=""
    df['Completness(%)']=""
    dfs = []

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
        print(f"DEBUG: Sample list contains {len(samplelist)} samples.")
        bwa_ref= f"bwa index {Ref}"
        print(f"Running BWA for reference indexing: {bwa_ref}")
        subprocess.run(bwa_ref, shell=True, check=True)
      # BWA and other commands for each sample
        for sample in samplelist:
            print(f"DEBUG: Processing sample: {sample}")  
            # Construct file paths for R1 and R2 based on `run_bowtie` flag
            if run_bowtie:
                sample_r1 = os.path.join(input_dir,f"{sample}_unmapped_1.fastq.gz")
                sample_r2 =  os.path.join(input_dir,f"{sample}_unmapped_2.fastq.gz")
            else:
                sample_r1 = os.path.join(input_dir,f"{sample}_trimmed_R1.fastq.gz")
                sample_r2 = os.path.join(input_dir,f"{sample}_trimmed_R2.fastq.gz")

            # Debugging: Print paths to ensure correctness
            print(f"Sample R1: {sample_r1}")
            print(f"Sample R2: {sample_r2}")

            sample_dir = f"{scientific_name}_assembled1"
            os.makedirs(sample_dir, exist_ok=True)
            bam_file = os.path.join(sample_dir, f"{sample}_{Refname}_mapped_reads.bam")
            vcf_file = os.path.join(sample_dir, f"{sample}_{Refname}_variants.vcf")
            consensus_file = os.path.join(sample_dir, f"{sample}_{Refname}_consensus_genome.fa")

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
            def calculate_genome_length(fasta_file):
           # Command to calculate genome length (valid nucleotides only: A, C, T, G)
               command = f"grep -v '^>' {fasta_file} | tr -d '\\n' | tr -cd 'ACTG' | wc -c"
               result = subprocess.run(command, shell=True, capture_output=True, text=True)
               return int(result.stdout.strip())
             # Calculate lengths
            try:
             ref_genome_len = calculate_genome_length(Ref)
             consensus_genome_len = calculate_genome_length(consensus_file)

             print(f"{consensus_file} my consensus")    
   # Save lengths to a text file
             with open(f"{sample_dir}/{sample}_{Refname}.txt", "w") as f:
              f.write(f"Reference Genome Length: {ref_genome_len} bp\n")
              f.write(f"Consensus Genome Length: {consensus_genome_len} bp\n")

        # print(f"Genome lengths saved to: {output_file}")
            except Exception as e:
             print(f"Error calculating genome lengths: {e}")
 
        # Add a new column and conditionally set values
        dftax['Ref_len'] = dftax.apply(lambda row: ref_genome_len if row['NCBI_ID'] == tax and row['SampleID'] == sample else row['Ref_len'],axis=1)
        dftax['Consensus_len'] = dftax.apply(lambda row: consensus_genome_len if row['NCBI_ID'] == tax and row['SampleID'] == sample else row['Consensus_len'], axis=1)
        dftax['Completness%'] = dftax.apply(lambda row:f"{consensus_genome_len / ref_genome_len:.2f}" if row['NCBI_ID'] == tax and row['SampleID'] == sample else row['Completness%'], axis=1) 

merged_df = pd.concat(dfs, ignore_index=True)
merged_df.to_csv("Output-summary.csv", index=False)
