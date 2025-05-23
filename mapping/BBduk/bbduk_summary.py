import glob
import sys
import pandas as pd
import os

direc = "/storage1/fs1/kutluay/Active/data/Outputs/Abhishek_Proj_3_dedup/Expts-5-CLIP-seq-of-Proj_3/mapping_Dedup/BBDuk"    # <<< BBDuk folder, no slash


os.chdir(direc + "/logs")

df = pd.DataFrame(columns=['Input Reads', 'Adapter Removed', 'Adapter Removed Percentage'])

adapter_removal = 'ATlog.txt'

input_reads = 0
adapter_removed = 0
adapter_removed_precent = 0
BC_files = []
filtered_reads = []
filtered_percent = []
trimmed_reads = []
trimmed_percent = []


with open(adapter_removal, 'r') as f:
    for line in f:
        if "Input:" in line:
            input_reads = int(line.split("\t")[1].split(" ")[0])
        elif "Result:" in line:
            adapter_removed = int(line.split("\t")[1].split(" ")[0])
            break

adapter_removed_precent = adapter_removed / input_reads


BC_files = [x.split("r")[1].split(".")[0] for x in glob.glob("filter*")]
BC_files = sorted(BC_files)
print(BC_files)

for filter_log in ["filter" + x + ".txt" for x in BC_files]:
    print(filter_log)

    with open(filter_log, 'r') as f:
        for line in f:
            if "Total Removed:" in line:
                fr = int(line.split("\t")[1].split(" ")[0])
                filtered_reads.append(fr)
                break

for trim_log in ["trim" + x + ".txt" for x in BC_files]:
    print(trim_log)

    with open(trim_log, 'r') as f:
        for line in f:
            if "Result:" in line:
                tr = int(line.split("\t")[1].split(" ")[0])
                trimmed_reads.append(tr)
                break

filtered_percent = [x / adapter_removed for x in filtered_reads]
trimmed_percent = [trimmed_reads[x] / filtered_reads[x] for x in range(len(filtered_reads))]


df["BC Files"] = BC_files
df.iat[0, 0] = input_reads
df.iat[0, 1] = adapter_removed
df.iat[0, 2] = adapter_removed_precent
df["Filtered Reads"] = filtered_reads
df["Filtered Percentage"] = filtered_percent
df["Trimmed Reads"] = trimmed_reads
df["Trimmed Percentage"] = trimmed_percent


df.to_csv(direc + "/" + "BBduk_summary.txt", sep="\t")
