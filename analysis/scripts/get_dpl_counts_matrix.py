import os
import pandas as pd

#initialize dictionary to hold
tpm_matrix = {}
raw_counts_matrix = {}

#load kallisto output directories
kallisto_dir = '/home/jdubos2/mtstp/data/dpl_kallisto_quantifications'
samples = os.listdir(kallisto_dir)

for sample in samples:

    tpm_matrix[sample] = {}
    raw_counts_matrix[sample] = {}

    abundance_files = f'{kallisto_dir}/{sample}/abundance.tsv'
    with open(abundance_files, 'r') as infile:
        lines = infile.readlines()
    
    for line in lines[1:]:
        gene_id = line.split('\t')[0]
        raw_estimated_count = line.split('\t')[3]
        tpm_count = line.split('\t')[4].strip()

        tpm_matrix[sample][gene_id] = tpm_count
        raw_counts_matrix[sample][gene_id] = raw_counts_matrix

tpm_df = pd.DataFrame(tpm_matrix)
transposed_tpm_df = tpm_df.transpose()
transposed_tpm_df.to_csv('/home/jdubos2/mtstp/data/counts_tables/dpl_tpm_counts_kallisto.csv')

raw_df = pd.DataFrame(raw_counts_matrix)
transposed_raw_df = raw_df.transpose()
transposed_raw_df.to_csv('/home/jdubos2/mtstp/data/counts_tables/dpl_raw_counts_kallisto.csv')
