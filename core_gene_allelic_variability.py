import pandas as pd
import os
import Bio.SeqIO as IO
import statistics

os.chdir('/path/to/folder/containing/fasta/file/of/core/genes')

infile = sys.argv[1] #a text file listing filenames of core gene fasta files
outfile = sys.argv[2]

filenames = []

with open(infile, 'r') as n:
    lines = n.readlines()

for line in lines:
    filenames.append(line.split('\n')[0])

med = []
lociname = []
count = []

for file in filenames:
    record_dict = IO.to_dict(IO.parse(file, 'fasta'))
    length = []
    for record in record_dict.items():
        length.append(len(record[1].seq))
    count.append(len(length))
    lociname.append(str(file).split(".")[0])
    med.append(statistics.median(length))

core_loci_var = pd.DataFrame(
    {
        'core_locus': lociname,
        'locus_length_median': med,
        'num_alleles': count
    }
)

core_loci_var.to_csv(outfile)
