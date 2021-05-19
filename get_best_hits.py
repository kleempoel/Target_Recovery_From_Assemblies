import sys
import pandas as pd

# Main
input_blast=sys.argv[1]
output_targets=sys.argv[2]
output_contigs=sys.argv[3]
output_blast=sys.argv[4]
print(input_blast, output_targets, output_contigs)

#load blast result
sample_blast=pd.read_table(input_blast, sep='\t', header=None) 
sample_blast.columns=['qseqid', 'sseqid', 'pident', 'length', 'bitscore']
sample_blast[['ref_sp', 'target']] = sample_blast['qseqid'].str.split('-', n=1, expand=True)
print(sample_blast.shape[0], sample_blast.target.nunique())

#Keep highest bitscore per target
target_blast = sample_blast.sort_values('bitscore',ascending=False).groupby('target').head(1).reset_index().drop(columns='index')
try:
	target_blast.sseqid=target_blast.sseqid.str.split('|',expand=True)[1]
except:
	print('no change to sseqid')
print(target_blast.shape[0], target_blast.target.nunique())

#Save sequence lists
target_blast.qseqid.drop_duplicates().to_csv(output_targets,header=None,index=False)
target_blast.sseqid.drop_duplicates().to_csv(output_contigs,header=None,index=False)
target_blast.to_csv(output_blast,index=False)