import os
import json

#get dpl reads list files
dpl_reads_list_path = '/home/jdubos2/mtstp/data/dpl_aligned_reads'
dpl_read_list_files = os.listdir(dpl_reads_list_path)
dpl_read_list_files = [file_name for file_name in dpl_read_list_files if 'ii' in file_name or 'ci' in file_name]

#get oe reads list files
oe_reads_list_path = '/home/jdubos2/mtstp/data/oe_aligned_reads'
oe_read_list_files = os.listdir(oe_reads_list_path)

#get base names for each oe reads file
samples = []
for file in oe_read_list_files:
    base_name = file.split('_')[0]
    samples.append(base_name)


#define dictionary to store ambiguous reads
ambiguous_reads = {}

#iterate through each pair of oe and dpl files and store contents in a list
for sample in samples:
    print(f"Comparing read ids from {sample}")
    ambiguous_reads[sample] = []

    #open files
    oe_file = f"{oe_reads_list_path}/{sample}_reads_mapped.txt"
    print(f"Reading: {oe_file}")
    with open(oe_file, 'r') as infile:
        oe_read_ids = infile.readlines()

    dpl_file = f"{dpl_reads_list_path}/{sample}_reads_mapped.txt"
    print(f"Reading: {dpl_file}")
    with open(dpl_file, 'r') as infile:
        dpl_read_ids = infile.readlines()
    
    #iterate through each oe read id and check if it is also in the monarch read ids
    for read_id in oe_read_ids:
        if read_id in dpl_read_ids or f"{read_id}\n" in dpl_read_ids or read_id.strip() in dpl_read_ids:
            read_id = read_id.strip()
            ambiguous_reads[sample].append(read_id)

    print(f"{sample}: {len(ambiguous_reads[sample])} ambiguous reads identified")

#define output file
outfile_name = '/home/jdubos2/mtstp/data/ambiguous_reads.json'
#write results
with open(outfile_name, "w") as outfile:
    json.dump(ambiguous_reads, outfile)