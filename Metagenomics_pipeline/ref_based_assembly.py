import os
import subprocess
import pandas as pd

def extract_sequence(fasta_file):
    """Extracts the sequence from a FASTA file as a single string."""
    try:
        with open(fasta_file, "r") as f:
            lines = f.readlines()
            sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
        return sequence
    except Exception as e:
        print(f"Error extracting sequence from {fasta_file}: {e}")
        return ""

def ref_based(df, run_bowtie, input_dir):
    base_dir = "Fasta_files"
    os.makedirs(base_dir, exist_ok=True)
    
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
        fasta_file = os.path.join(tax_dir, f"{scientific_name}.fasta")
        
        command = f'esearch -db nucleotide -query "txid{tax}[Organism]" | efilter -source refseq | efetch -format fasta > {fasta_file}'
        subprocess.run(command, shell=True, check=True)
        
        if not os.path.exists(fasta_file):
            print(f"FASTA file {fasta_file} was not created.")
            continue
        
        bwa_ref = f"bwa index {fasta_file}"
        subprocess.run(bwa_ref, shell=True, check=True)
        
        for sample in dftax['SampleID']:
            sample_r1 = os.path.join(input_dir, f"{sample}_unmapped_1.fastq.gz" if run_bowtie else f"{sample}_trimmed_R1.fastq.gz")
            sample_r2 = os.path.join(input_dir, f"{sample}_unmapped_2.fastq.gz" if run_bowtie else f"{sample}_trimmed_R2.fastq.gz")
            
            sample_dir = os.path.join(base_dir, f"{scientific_name}_assembled1")
            os.makedirs(sample_dir, exist_ok=True)
            
            bam_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_mapped_reads.bam")
            vcf_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_variants.vcf")
            consensus_file = os.path.join(sample_dir, f"{sample}_{scientific_name}_consensus_genome.fa")
            
            bwa_command = f"bwa mem -a -t 16 {fasta_file} {sample_r1} {sample_r2} | samtools view -u -@ 3 - | samtools sort -@ 16 -o {bam_file}"
            subprocess.run(bwa_command, shell=True, check=True)
            subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
            subprocess.run(f"bcftools mpileup -f {fasta_file} {bam_file} | bcftools call -c --ploidy 1 -v -o {vcf_file}", shell=True, check=True)
            subprocess.run(f"samtools mpileup -aa -A -d 0 -Q 0 -f {fasta_file} {bam_file} | ivar consensus -p {consensus_file}", shell=True, check=True)
            
            def calculate_length(fasta_file):
                command = f"grep -v '^>' {fasta_file} | tr -d '\n' | tr -cd 'ACTG' | wc -c"
                result = subprocess.run(command, shell=True, capture_output=True, text=True)
                return int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0
            
            try:
                ref_len = calculate_length(fasta_file)
                consensus_len = calculate_length(consensus_file)
                completeness = round((consensus_len / ref_len) * 100, 2) if ref_len > 0 else 0
                sequence = extract_sequence(consensus_file)
                
                dftax.loc[dftax['SampleID'] == sample, ['Ref_len', 'Consensus_len', 'Completeness(%)', 'sequence']] = [ref_len, consensus_len, completeness, sequence]
            except Exception as e:
                print(f"Error processing sample {sample}: {e}")
        
        dfs.append(dftax)
    
    merged_df = pd.concat(dfs, ignore_index=True)
    
    
    filtered_df = merged_df[pd.to_numeric(merged_df['Completeness(%)'], errors='coerce') >= 60]
    filtered_df.to_csv("Output-summary_complete.csv", index=False)
    merged_df.drop(columns=['sequence'], inplace=True)
    merged_df.to_csv("Output-summary1.csv", index=False)
  
