import os
import subprocess
import pandas as pd
from pathlib import Path



def split_fasta(input_file, output_dir):
    """
    Splits a multi-sequence FASTA file into individual FASTA files, 
    each corresponding to a separate sequence, and returns a list of 
    those file paths.

    Args:
    - input_file (str): Path to the input FASTA file.
    - output_dir (str): Directory to store the individual FASTA files.

    Returns:
    - List of paths to the generated FASTA files.
    """
    Path(output_dir).mkdir(exist_ok=True)  # Ensure the output directory exists
    fasta_files = []  # List to store the paths of the generated FASTA files

    with open(input_file, "r") as fasta:
        sequence = []
        accession = None
        for line in fasta:
            if line.startswith(">"):  # New sequence header
                if accession:  # Save the previous sequence
                    file_path = f"{output_dir}/{accession}.fasta"
                    with open(file_path, "w") as out_fasta:
                        out_fasta.write("".join(sequence))
                    fasta_files.append(file_path)  # Add to list of FASTA file paths
                # Extract accession number from header
                accession = line.split()[0][1:]  # Remove '>' and take first part
                sequence = [line]  # Start a new sequence
            else:
                sequence.append(line)
        # Save the last sequence
        if accession:
            file_path = f"{output_dir}/{accession}.fasta"
            with open(file_path, "w") as out_fasta:
                out_fasta.write("".join(sequence))
            fasta_files.append(file_path)  # Add to list of FASTA file paths

    return fasta_files





def get_best_reference(sample_r1, sample_r2, reference_list):
    """
    Align paired-end FASTQ files to a list of reference FASTA files using BWA
    and return the best reference based on alignment scores.

    Parameters:
        sample_r1 (str): Path to the first paired-end FASTQ file.
        sample_r2 (str): Path to the second paired-end FASTQ file.
        reference_list (list): List of paths to reference FASTA files.

    Returns:
        str: The best reference file with the highest alignment score.
    """
    # Dictionary to store alignment scores
    alignment_scores = {}
   
    for fasta in reference_list:
        index_base = Path(fasta).stem  # Extract the base name for the index
        output_sam = f"{index_base}_aligned.sam"  # SAM file for the output

        # Check if the BWA index exists; if not, create it
        if not Path(f"{fasta}.bwt").exists():
            print(f"Index for {fasta} not found. Creating index...")
            try:
                subprocess.run(["bwa", "index", fasta], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error in BWA index creation for {fasta}: {e}")
                continue

        # Run BWA for alignment
        bwa_command = [
            "bwa", "mem", fasta, sample_r1, sample_r2
        ]
        try:
            with open(output_sam, "w") as sam_file:
                subprocess.run(bwa_command, check=True, stdout=sam_file)
        except subprocess.CalledProcessError as e:
            print(f"Error in BWA alignment for {fasta}: {e}")
            continue

        # Parse alignment score from the SAM file
        score = 0
        try:
            with open(output_sam, "r") as sam_file:
                for line in sam_file:
                    if not line.startswith("@"):  # Skip header lines
                        fields = line.strip().split("\t")
                        try:
                            # Use the AS:i:<score> tag for alignment score (if available)
                            score += int([tag for tag in fields if tag.startswith("AS:i:")][0].split(":")[2])
                        except IndexError:
                            continue
        except FileNotFoundError:
            print(f"Failed to open SAM file: {output_sam}")
            continue

        alignment_scores[fasta] = score

    # Pick the reference with the highest alignment score
    if alignment_scores:
        best_reference = max(alignment_scores, key=alignment_scores.get)
        print(f"The best reference is {best_reference} with a score of {alignment_scores[best_reference]}")
        return best_reference
    else:
        print("No alignments were successful.")
        return None

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
    df["Accession_number"]=""
    df['sequence'] = ""
    dfs = []

    for tax in taxids:
        dftax = df[df['NCBI_ID'] == tax].copy()
        #scientific_name = dftax['Scientific_name'].iloc[0].replace(' ', '_')
        scientific_name=dftax['Scientific_name'].iloc[0].replace(' ', '_').replace('/', '_')
        
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
                 # Example usage
            input_files = fasta_file
            #output_dir = f"{output_dir}/{sample}_{scientific_name}"
            reference_list = split_fasta(input_files, sample_dir)
            fasta_file=get_best_reference(sample_r1, sample_r2, reference_list)
            cmd = f"grep '^>' {fasta_file} | cut -d ' ' -f1 | sed 's/^>//'"
            acc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

            # To get the output as a list of sequence headers:
            acc_ids = acc.stdout.strip().split("\n")[0]
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
                
                dftax.loc[dftax['SampleID'] == sample, ['Ref_len', 'Consensus_len', 'Completeness(%)', 'Accession', 'sequence']] = [ref_len, consensus_len, completeness, acc, sequence]
            except Exception as e:
                print(f"Error processing sample {sample}: {e}")
        
        dfs.append(dftax)
    
    merged_df = pd.concat(dfs, ignore_index=True)
    
    
    filtered_df = merged_df[pd.to_numeric(merged_df['Completeness(%)'], errors='coerce') >= 60]
    filtered_df.to_csv("Output-summary_complete.csv", index=False)
    merged_df.drop(columns=['sequence'], inplace=True)
    merged_df.to_csv("Output-summary1.csv", index=False)
  
