from Bio import SeqIO
import pandas as pd

seq_stats=pd.DataFrame(columns=['Sample','Ntargets','SumContigLength'])
idx=0
for filename in snakemake.input:
	print(filename)
	records = list(SeqIO.parse(filename, "fasta"))
	seq_stats.loc[idx,'Sample']=filename
	seq_stats.loc[idx,'Ntargets']=len(records)
	seq_length=0
	for record in records:
		seq_length += len(record)
	seq_stats.loc[idx,'SumContigLength']=seq_length
	print(seq_stats.loc[idx,:])
	idx += 1
seq_stats.to_csv(snakemake.output[0],index=False)