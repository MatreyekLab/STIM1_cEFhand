import os
import pandas as pd
import csv
import codecs
import sys
import time
import warnings
from multiprocessing import Pool, freeze_support

genetic_code = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}
domain_away_aa = 31

warnings.simplefilter(action='ignore', category=FutureWarning)

WT_dna = 'GCGACAGGAACCAGCTCGGGGGCCAACTCTGAGGAGTCCACTGCAGCAGAGTTTTGCCGAATTGACAAGCCCCTGTGTCACAGTGAGGATGAGAAACTCAGCTTCGAGGCAGTCCGTAACATCCACAAACTGATGGACGATGATGCCAATGGTGATGTGGATGTGGAAGAAAGTGATGAGTTCCTGAGGGAAGACCTCAATTACCATGACCCAACAGTGAAACACAGCACCTTCCATGGTGAGGATAAGCTCATCAGCGTGGAGGACCTGTGGAAGGCATGGAAGTCATCAGAAGTATACAATTGGACCGTGGATGAGGTGGTA'

def process_fastq_file(file_path):
    index = 1
    read = 0
    seq_index = 1  # Initialize the sequence index
    processed_seq_indices = set()  # Set to store processed seq_index values
    if os.path.exists(file_path):
        print(f'Processing file: {file_path}')
        outfile_name = str(file_path[:len(file_path)-6])+"_codon.tsv"
        outfile = codecs.open(outfile_name, "w", "utf-8", "replace")
        outfile.write("read_class\tvariant\tseq_index\n")  # Add "seq_index" header
        with open(file_path, 'rt') as csvfile:
            csvfile = csv.reader(csvfile, delimiter=' ', quotechar='|')
            for row in csvfile:
                if index == 1:
                    current_header = row[0]
                    index += 1
                    continue
                if index == 2:
                    current_sequence = row[0]
                    index += 1
                    continue
                if index == 3:
                    index += 1
                    continue
                if index == 4:
                    current_qscore = row[0]
                    read += 1
                    index = 1
                    if current_sequence == WT_dna:
                        if seq_index not in processed_seq_indices:
                            outfile.write("WT\tZ0Z\t{}\n".format(seq_index))
                            processed_seq_indices.add(seq_index)
                        seq_index += 1
                    else:
                        read_class = "Variant"
                        sequence = current_sequence
                        for i in range(0, len(sequence), 3):
                            codon = sequence[i:i+3]
                            reference_codon = WT_dna[i:i+3]
                            position = (i // 3) + domain_away_aa
                            if codon != reference_codon:
                                amino_acid = genetic_code.get(codon, '?')
                                reference_amino_acid = genetic_code.get(reference_codon, '?')
                                variant = reference_amino_acid + str(position) + amino_acid
                                outfile.write("Variant\t{}\t{}\n".format(variant, seq_index))
                                if seq_index not in processed_seq_indices:
                                    processed_seq_indices.add(seq_index)
                        seq_index += 1
                        continue

        outfile.close()
        print("Starting codon and variant count")
        outfile2_name = str(file_path[:len(file_path)-6])+"_summary.tsv"
        #outfile2 = codecs.open(outfile2_name, "w", "utf-8", "replace")
        results_df = pd.read_csv(outfile_name, delimiter='\t', header=0)
        fil_results = results_df.drop_duplicates(subset='seq_index', keep=False)
        variant_count = fil_results['variant'].value_counts().reset_index()
        variant_count.columns = ['variant', 'count']
        variant_count.to_csv(outfile2_name, sep='\t') 
        print('Nucleotide and Variant counts done, file in directory')


    else:
        print(f'File {file_path} does not exist')
        return

if __name__ == '__main__':
    freeze_support()

    # Enter the directory path where your FASTQ files are located
    fastq_dir = '/Users/nishakamath/Downloads/eflib_amp'

    # Get the list of FASTQ files in the directory
    fastq_files = [os.path.join(fastq_dir, file) for file in os.listdir(fastq_dir) if file.endswith('.fastq')]

    # Set the number of processes to run in parallel (adjust as needed)
    num_processes = 5

    # Start the timer
    start_time = time.time()

    # Create a pool of worker processes
    pool = Pool(processes=num_processes)

    # Process each FASTQ file in parallel
    pool.map(process_fastq_file, fastq_files)

    # Close the pool
    pool.close()
    pool.join()

    end_time = time.time()
    execution_time = (end_time - start_time) / 60
    print('Execution Time:', execution_time, 'min')

