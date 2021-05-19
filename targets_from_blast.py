import sys
import pandas as pd

def blast_results(input_blast,qs=True):
    if qs:
        qseqid='qseqid'; sseqid='sseqid'; qlen='qlen'; qstart='qstart'; 
        qend='qend'; slen='slen'; sstart='sstart'; send='send'
        
    else:
        qseqid='sseqid'; sseqid='qseqid'; qlen='slen'; qstart='sstart'; 
        qend='send'; slen='qlen'; sstart='qstart'; send='qend'
        
    #load blast result
    sample_blast_raw=pd.read_table(input_blast, sep='\t', header=None) 
    sample_blast_raw.columns=[qseqid, sseqid, 'pident', 'length', 'mismatch', 'gapopen', qlen, qstart, 
                              qend, slen, sstart, send, 'evalue', 'bitscore','qseq']
    print(sample_blast_raw.shape[0], sample_blast_raw.sseqid.nunique())

    #keep only per gene based on highest evalue
    sample_blast=sample_blast_raw.sort_values('bitscore',ascending=False).groupby(qseqid).head(1)\
                    .sort_values(['qseqid']).reset_index().drop(columns='index')
    
    #split sseqid to get ref species and ref gene
    sample_blast[['ref_sp', 'starget']] = sample_blast['sseqid'].str.split('-', n=1, expand=True)
    #remove 'g' in gene name for compatibility with paftol
    #try: 
    #    sample_blast['starget']=sample_blast['starget'].str.replace('g', '')
    #except:
    #    sample_blast['starget']=sample_blast['starget']
    #try: 
    #    sample_blast['qseqid']=sample_blast['qseqid'].str.replace('g', '')
    #except:
    #    sample_blast['qseqid']=sample_blast['qseqid']
    
    #calculate overlap
    sample_blast['cover_reference_pc']=round((abs(sample_blast.sstart-sample_blast.send)+1)/sample_blast.slen*100,1)
    sample_blast['cover_retrieved_pc']=round((abs(sample_blast.qstart-sample_blast.qend)+1)/sample_blast.qlen*100,1)
    sample_blast['seqid_match']=sample_blast.qseqid==sample_blast.starget
    
    return sample_blast[['qseqid', 'sseqid', 'starget', 'seqid_match', 'ref_sp', 'pident', 'evalue', 
                         'bitscore','length', 'qlen', 'slen', 'qstart', 'qend', 'sstart','send',
                        'cover_reference_pc', 'cover_retrieved_pc', 'qseq']]
						

def writefasta(output_fasta, df_sample, keep_align=False):
    file1 = open(output_fasta,"w")
    for index, iseq in df_sample.iterrows():
        file1.write('>' + iseq.starget + ' ' + iseq.ref_sp.split('_')[0] 
                    + ' organism-gene:' + iseq.ref_sp.split('_')[0] + '-' + iseq.starget
                    + ' originalID:' + iseq.qseqid
                    + ' length:' + str(iseq.length)
                    + ' pident:' + str(round(iseq.pident,2)) + ' evalue=' + str(iseq.evalue)
                    + ' cover_reference_pc:' + str(iseq.cover_reference_pc)
                    + ' cut2codon:' + str(iseq.cut2codon) 
                    + ' sstart:' + str(iseq.sstart) 
                    + ' send:' + str(iseq.send) 
                    + '\n') 
        
		#Sequences are aligned, choice to removed alignment information
        if keep_align==False:
            file1.write(iseq.qseq.replace('-', '') + '\n') 
        else:
            file1.write(iseq.qseq + '\n') 
    file1.close()
    
def filter_to_fasta(output_fasta, df_sample, 
                    min_identity=80, min_ref_overlap=0, min_length=50, write_fasta=True, keep_align=False):
    
	#trim out-of-frame nucleotides
    df_sample['cut2contig']=0
    for idx in range(0,df_sample.shape[0]):
        #if sequence does not cover the reference entirely
        if df_sample.loc[idx,'sstart']>df_sample.loc[idx,'send']:
            df_sample.loc[idx,'qseq']=str(Seq(df_sample.loc[idx,'qseq']).reverse_complement())
        if df_sample.loc[idx,'cover_reference_pc']<100:
            ##assumes all sequences are in the same 5'-3' direction as the reference
            seq_ini=min([df_sample.loc[idx,'sstart'],df_sample.loc[idx,'send']])
            #e.g. if ini is 4 then its fine as 4-1%3=1
            to_contig=seq_ini
            diff_ini=to_contig-seq_ini
            while (to_contig-1)%3!=0:
                to_contig=to_contig+1
                diff_ini=to_contig-seq_ini
            qseq=str(df_sample.loc[idx,'qseq'])
            df_sample.loc[idx,'qseq']=qseq[diff_ini:]
            df_sample.loc[idx,'length']=len(qseq[diff_ini:])
            df_sample.loc[idx,'cut2codon']=diff_ini
#             print(seq_ini, to_contig,diff_ini, len(qseq),df_sample.loc[idx,'length'])
    
    #filter
    Filtered_sample=df_sample[(df_sample.pident>min_identity) & (df_sample.cover_reference_pc>min_ref_overlap) 
                                 & (df_sample.length>min_length) & (df_sample.seqid_match==True)]
    
    #save to fasta
    if write_fasta:
        if Filtered_sample.shape[0]>0:
            writefasta(output_fasta=output_fasta, df_sample=Filtered_sample, keep_align=keep_align)
    
    return Filtered_sample

# Main
input_blast=sys.argv[1]
output_fasta=sys.argv[2]
print(input_blast, output_fasta)

df_blast=blast_results(input_blast=input_blast, qs=False)
df_blast['seqid_match']=True
print(df_blast.shape[0], df_blast.starget.nunique())
df_blast = df_blast.sort_values('bitscore',ascending=False).groupby('starget').head(1).reset_index().drop(columns='index')
print(df_blast.shape[0], df_blast.starget.nunique())
df_filter=filter_to_fasta(output_fasta=output_fasta, df_sample=df_blast, write_fasta=True, 
						 min_identity=70, min_ref_overlap=0, min_length=50,keep_align=False)
print(input_blast,', Ntargets:',df_blast.shape[0],'SumContigLength:',df_blast.length.sum())