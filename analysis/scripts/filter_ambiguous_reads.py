#!/usr/bin/env python3

import pyfastx
import os

def get_unique_ids(file1, file2):

    #dictionary to hold ids
    id_dict = {}
    
    #process file 1
    with open(file1, 'r') as infile:
        lines_file1 = infile.readlines()
    lines_file1[-1] = f"{lines_file1[-1]}\n"

    #process file 2
    with open(file2, 'r') as infile:
        lines_file2 = infile.readlines()
    lines_file2[-1] = f"{lines_file2[-1]}\n"

    #add all ids from file 1 (dpl)
    for id in lines_file1:
        id = id.strip()
        id_dict[id] = ''

    #if id in set 2 is not in the current id_dict, add it
    #if id in set 2 is in current id_dict, subtract id from dict
    for id in lines_file2:
        id = id.strip()
        if id not in id_dict:
            id_dict[id] = ''
        elif id in id_dict:
            del id_dict[id]
    
    print(f"Total unique ids: {len(id_dict)}")

    return id_dict
    
def filter_sequences(sequences_ids_to_keep, infile, outfile):
    sequences_to_keep = []

    #load file
    fq_file = pyfastx.Fastx(infile, comment=True)

    filtered = 0
    kept = 0
    #get sequences to keep
    for name, seq, qual, comment in fq_file:
        if name in sequences_ids_to_keep:
            seq_entry = f"@{name}\n{seq}\n+\n{qual}\n"
            sequences_to_keep.append(seq_entry)
            kept += 1
        else:
            filtered += 1

    print(f"sample:{outfile} kept:{kept} ambiguous:{filtered}")
    
    #remove newline on end of last record
    sequences_to_keep[-1] = sequences_to_keep[-1].strip()

    #write to outfile
    print(f"Writing output")
    with open(outfile, 'w') as out:
        for sequence in sequences_to_keep:
            out.write(sequence)

    print("Compressing output")
    #zip outfile
    os.system(f"gzip {outfile}")


def run_filtering(oe_aligned_reads_dir):
    #load files
    oe_files = os.listdir(oe_aligned_reads_dir)
    
    for oe_file in oe_files:

        sample_id = oe_file.split('_')[0]

        #check to see if sample is done or is being worked on
        working_or_complete = 0
        done_files = os.listdir('/home/jdubos2/mtstp/data/unambiguous_reads')

        for file in done_files:
            if sample_id in file:
                working_or_complete += 1

        #if sample is not running or completed, run
        if working_or_complete == 0:
            print(f"Working on {sample_id}")

            #add placeholder in outdir to indicate the sample is being worked on
            os.system(f"touch /home/jdubos2/mtstp/data/unambiguous_reads/{sample_id}")

            #get uniquely mapped reads
            oe_mapped_reads = f'/home/jdubos2/mtstp/data/oe_aligned_reads/{sample_id}_reads_mapped.txt'
            dpl_mapped_reads = f'/home/jdubos2/mtstp/data/dpl_aligned_reads/{sample_id}_reads_mapped.txt'
            unique_ids = get_unique_ids(dpl_mapped_reads, oe_mapped_reads)

            #assemble files
            forward_seq_file = f'/home/jdubos2/mtstp/data/raw_sequences/{sample_id}_1.fq.gz'
            reverse_seq_file = f'/home/jdubos2/mtstp/data/raw_sequences/{sample_id}_2.fq.gz'
            forward_outfile = f'/home/jdubos2/mtstp/data/unambiguous_reads/{sample_id}_unambiguous_1.fq'
            reverse_outfile = f'/home/jdubos2/mtstp/data/unambiguous_reads/{sample_id}_unambiguous_2.fq'

            #filter sequences
            filter_sequences(unique_ids, forward_seq_file, forward_outfile)
            filter_sequences(unique_ids, reverse_seq_file, reverse_outfile)

            #remove placeholder when sample is done
            os.system(f"rm /home/jdubos2/mtstp/data/unambiguous_reads/{sample_id}")

run_filtering('/home/jdubos2/mtstp/data/oe_aligned_reads')