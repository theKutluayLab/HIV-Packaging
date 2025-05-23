import glob
import sys
import pandas as pd
import os
import re

# Parses BBDuk and Bowtie output files and produces summary statistics

direc = '/path/to/experiment_folder'   # <<< directory of the whole experiment folder, with slash
direc_BBDuk = '/path/to/BBDuk_output_folder' #with slash

df = pd.DataFrame(columns=['Input Reads', 'After Adapter Trimming'])
bbduk_output = direc_BBDuk + 'logs/ATlog.txt'

input_reads = ""
output_reads = ""
BC_files = []
total_counts = []
aligned = []

with open(bbduk_output, 'r') as f:
    for line in f:
        if 'Input:' in line:
            input_reads = int(re.findall(r"[\w']+", line)[1])
        if 'Result:' in line:
            output_reads = int(re.findall(r"[\w']+", line)[1])
            break

direc = direc + "Bowtie_output_folder/" #change, with slash, hgenome

BC_folders = [o for o in os.listdir(direc) if os.path.isdir(os.path.join(direc, o))]
BC_names = []

for folder in BC_folders:
    full_path = os.path.join(direc, folder)
    BC_files.append(folder)
    print(folder + "\n")
    
    BC_subfolders = [o for o in os.listdir(full_path) if o.startswith("BC") and os.path.isdir(os.path.join(full_path, o))]
    
    for bc_folder in BC_subfolders:
        rrna_alignment_paths = glob.glob(os.path.join(full_path, bc_folder, "**", "bowtie_mapping_output.txt"), recursive=True)
        BC_names.append(bc_folder)
        
        if rrna_alignment_paths:
            rrna_alignment_path = rrna_alignment_paths[0] # assuming there's only one such file
            with open(rrna_alignment_path, 'r') as f:
                for line in f:
                    if 'reads processed' in line:
                        print(line)
                        total_counts.append(int(line.split(" ")[3]))
                    elif 'at least one alignment' in line:
                        print(line)
                        aligned.append(int(line.split(" ")[7]))
                        break
        else:
            print(f"No bowtie_mapping_output.txt in {os.path.join(full_path, bc_folder)}")

percentage_aligned = [aligned[x] / total_counts[x] for x in range(len(total_counts))]


df['Files'] = BC_names
df.iat[0, 0] = input_reads
df.iat[0, 1] = output_reads
df['Post BBDuk'] = total_counts
df['Aligned Reads'] = aligned
df['Aligned Percentage'] = percentage_aligned


df.to_csv(direc + 'Read_Count_Summary.txt', sep='\t') #change
