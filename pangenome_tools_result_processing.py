#Consolidating core genome detected from PIRATE pangenome tool* with result from other pangenome tool(s)*
##by Krisna, M & Monteith, W 2023

## *Also applicable for result from panaroo, PEPPAN, and chewBBACA with modifications

import pandas as pd
import os

infile = sys.argv[1]
  #infile should be formatted to a format detailed in https://www.protocols.io/private/EF6DB7FE429311EEB8630A58A9FEAC02
outfile = sys.argv[2]

df_p2 = pd.read_csv(infile)
core_gene_list_p2 = df_p2['locus_id'].to_list()

os.chdir('/path/to/modified/gffs/result/')
files_p2 = os.listdir()

result_p2 = {}
for locus in core_gene_list_p2:
    for filename in files_p2:
        with open(filename, 'r') as f:
            content = f.readlines()
        for line in content:
            if locus in line:
                if locus not in result_p2:
                    result_p2[locus] = line 
            
if len(result_p2) != len(core_gene_list_p2):
    for locus in core_gene_list_p2:
        if locus not in result_p2:
            print(locus, 'not found in gffs')

start_loc_p2 = []
end_loc_p2 = []
note = []
for locus in core_gene_list_p2:
    result = result_p2[locus] 
    result = result.split('\t')
    start_loc_p2.append(result[3])
    end_loc_p2.append(result[4])
    note.append(result[8])

note2 = []
prev_locus = []
for data in note:
    split_data = data.split('prev_locus=',1)[1]
    note2.append(split_data)

for data in note2:
    split_data1 = data.split(';')
    prev_locus.append(split_data1[0])

pirate_output = pd.DataFrame (
    {
        'gene_family' : df_p2['gene_family'],
        'locus_id' : df_p2['locus_id'],
        'isolate_id' : df_p2['isolate_id'],
        'start_loc' : start_loc_p2,
        'end_loc' : end_loc_p2,
        'prokka_locus_name' : prev_locus
    }
)

pirate_output.to_csv(outfile)
