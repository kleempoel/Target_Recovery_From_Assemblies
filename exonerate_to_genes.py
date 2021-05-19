#!/usr/bin/env python
# coding: utf-8

# In[33]:


from Bio import SearchIO, SeqIO
import pandas as pd
import os
import sys


# In[34]:


in_alignment=sys.argv[1]
in_targets=sys.argv[2]
out_fasta=sys.argv[3]
out_stats=sys.argv[4]
# in_alignment='Alignments/GCA_900406275.1.aln'
# in_targets='../../GitHub_Repositories/Angiosperms353_targetSequences_organism-gene_format_corrected.fasta'
# out_fasta='test.fasta'
# out_stats='test_stats.csv'


# ## Functions

# In[35]:


def get_align(query_res):
    rd={}
    rd['Nhits']=len(query_res)
    # If more than one match, keep the one with the highest score
    if len(query_res)>1:
        hsp_scores=[]
        for ihsp in query_res:
            hsp_scores.append(ihsp[0].score)
        hsp=query_res[hsp_scores.index(max(hsp_scores))][0]
    else:
        hsp=query_res[0][0]
        
    rd['qseqid']=hsp.query_id
    rd['sseqid']=hsp.hit_id
    rd['Score']=hsp.score
    full_seq='';seq_sim='';strand=[]
    for ihsp in hsp:
        full_seq += ihsp.hit.seq
        seq_sim += str(ihsp.aln_annotation).replace("{'similarity': '",'').replace("'}",'')
        strand.append(ihsp.hit_strand)
    rd['Sequence']=str(full_seq).replace('-','').upper()
    rd['Ngaps']=str(full_seq).count('-')
    rd['length']=len(rd['Sequence'])
    rd['Nucl_same']=seq_sim.count('|')
    rd['Nucl_diff']=seq_sim.count(' ')
    rd['Nucl_minus']=seq_sim.count('-')
    rd['Similarity']=round(rd['Nucl_same']/(rd['Nucl_same']+rd['Nucl_diff']),2)
    rd['Strand']=max(set(strand), key = strand.count)
    inter_ranges=pd.DataFrame(hsp.hit_inter_ranges,columns=['ini','end'])
    inter_ranges['inter_len']=inter_ranges.end-inter_ranges.ini
    rd['inter_length']=inter_ranges.inter_len.sum()
    rd['Nfragments']=len(hsp)
    return rd


# In[36]:


def writefasta(out_file, df_target):
    file1 = open(out_file,"w")
    for idx, iseq in df_target.iterrows():
        file1.write('>' + iseq.target + ' ' + iseq.organism
                    + ' organism-gene:' + iseq.qseqid
                    + ' originalID:' + iseq.sseqid
                    + ' length=' + str(iseq.length)
                    + ' Score=' + str(iseq.Score)
                    + ' Similarity=' + str(iseq.Similarity)
                    + ' Ngaps=' + str(iseq.Ngaps) 
                    + ' Coverage=' + str(iseq.Coverage)
                    + ' Nfragments=' + str(iseq.Nfragments) 
                    + ' exon_intron_len_ratio=' + str(iseq.len_ratio) 
                    + ' Nhits=' + str(iseq.Nhits) 
                    + ' Strand=' + str(iseq.Strand) 
                    + '\n') 

        #Sequences are aligned, choice to removed alignment information
        file1.write(iseq.Sequence + '\n') 
    file1.close()


# ## Main

# In[43]:


if os.path.isfile(in_alignment):
    # Load data and ref targets
    all_qresult = list(SearchIO.parse(in_alignment, 'exonerate-text')) #Doesn't seem to mind the extra text
    print(in_alignment, len(all_qresult),end=', ')
    
    t353={}
    with open(in_targets) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            t353[record.id]=len(record.seq)
    t353 = pd.DataFrame.from_dict(t353,orient='index').reset_index().rename(columns={'index':'qseqid',0:'Ref_length'})

    # Process alignments
    dc_res=[get_align(x) for x in all_qresult]
    df_res=pd.DataFrame(dc_res)
    df_res[['organism','target']]=df_res.qseqid.str.split('-',expand=True)
    df_res['len_ratio']=round(df_res.length/df_res.inter_length,6)
    df_res['hsp_idx']=df_res.index
    print(df_res.shape[0],df_res.target.nunique(),df_res.organism.nunique())

    # Add target length and coverage
    df_res=pd.merge(df_res,t353[['qseqid','Ref_length']],how='left',on='qseqid')
    print('N missing data:',df_res.isna().sum().to_dict())
    df_res['Coverage']=round(df_res.length/df_res.Ref_length,2)

    # Filter
    col_keep=['target','sseqid','Nhits','qseqid','Strand',
          'Score','length', 'inter_length','len_ratio','Similarity','Coverage','Nfragments','Ngaps', 
       'Nucl_same','Nucl_diff','Nucl_minus','organism','hsp_idx','Ref_length','Sequence']
    print('\nOriginal:',df_res.shape[0],df_res.target.nunique())
    df_target=df_res[df_res['Coverage']>0.25]
    print('Filter by coverage:',df_target.shape[0],df_target.target.nunique())
    df_target=df_target[df_target['Similarity']>0.7]
    print('Filter by Similarity:',df_target.shape[0],df_target.target.nunique())
    df_target=df_target[df_target['len_ratio']>0.1]
    print('Filter by len_ratio:',df_target.shape[0],df_target.target.nunique())
    # df_target=df_target[df_target['Nhits']<=2]
    # print(df_target.shape[0],df_target.target.nunique())
    df_target=df_target.sort_values('Score',ascending=False).groupby('target').head(1)
    df_target=df_target[col_keep]
    print('N entries:',df_target.shape[0],', Nseq:',df_target.target.nunique(),
          ', SumContigsLength:',df_target.length.sum(), 'bp , Sum Score:', df_target.Score.sum(),'\n')
    print('Ref organisms:',df_target.groupby('organism').size().sort_values(0,ascending=False)[:20].to_dict(),'\n')
    print(df_target.describe())
#     df_target.length.hist(bins=30);

    # Save stats
    df_target.drop(columns=['Sequence']).to_csv(out_stats,index=False)

    # Write fasta
    writefasta(out_file=out_fasta, df_target=df_target)
else:
    print('no file',fname)

