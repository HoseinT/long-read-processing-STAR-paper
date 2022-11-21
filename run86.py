#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from align import aligner
import numpy as np
from pandarallel import pandarallel
import pickle
import re
import plotly_express as px
import seaborn as sns
import plotly.graph_objects as go
from subprocess import Popen
import os


# In[2]:


sheet1 =pd.read_csv('data/20210321_Indexes_Hybcap2.csv',sep='\t', header=1)


# In[3]:


sheet1


# In[4]:


sheet2 = pd.read_csv('data/STAR_Hybcapno3_Sample lists_forHosein.csv', sep='\t', skiprows=7, nrows=8).iloc[:,1:]


# In[5]:


sheet2.rename(columns={'Patient/Sample':'Hybrid capture sample list'}, inplace=True)


# In[6]:


sheet1['run'] = 'pt_072'
sheet2['run'] = 'pt_086'


# In[7]:


sheet = pd.concat([sheet1,sheet2]).drop(columns='Sample no.')


# In[8]:


sheet


# In[2]:


sample ='BCSA3'


# In[9]:


sheet = sheet[[sample in x for x in [re.sub('-','',y) for  y in sheet['Hybrid capture sample list']]]]


# In[10]:


sheet


# In[11]:


sheet = sheet.reset_index()


# In[12]:


sheet = sheet.assign(aux_id = ['A' if os.path.exists('data/{}/demultiplexed/demultiplex.bc{}_BAK8A_OA--bc{}_BAK8A_OA.hifi_reads.fastq'.format(x['run'],
x['PacBio index ID'],x['PacBio index ID'])) else 'B' for i,x in sheet.iterrows()]) 
sheet = sheet.assign(fastq_file = ['data/{}/demultiplexed/demultiplex.bc{}_BAK8{}_OA--bc{}_BAK8{}_OA.hifi_reads.fastq'.format(x['run'],
x['PacBio index ID'],x.aux_id,x['PacBio index ID'],x.aux_id) for i,x in sheet.iterrows()])


# In[13]:


sheet.fastq_file.apply(os.path.exists)


# In[14]:


def fastq_to_df(fastq_file, library_name):
    reads2 = []
    i = 0
    for record in SeqIO.parse(open(fastq_file,'rt'), 'fastq'):
        i+=1
        reads2.append({'seq':str(record.seq), 'qid':record.id, 'qual':record.letter_annotations['phred_quality']})

    df = pd.DataFrame(reads2)
    df['library'] = library_name
    return df


# In[15]:


readtab2 = pd.concat((fastq_to_df(x.fastq_file,x['Hybrid capture sample list']) for i,x in sheet.iterrows()))


# In[16]:


readtab2= readtab2.assign(rcseq = [str(Seq(x).reverse_complement()) for x in readtab2.seq])


# In[17]:



def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_adapter(readstring):
    A_seq = 'CGACGCTCTTCCGATCT'
    A_len = len(A_seq)
    distances = [(offset+A_len,hamming_distance(readstring[offset:offset+A_len], A_seq)) for offset in range(7)]
    argmin   = [x[1] for x in distances].index(min([x[1] for x in distances]))
    return distances[argmin]

def find_TSO(readstring):
    TSO_rc_seq = 'CCCATGTACTCTGCGTTGATACCACTGCTT'
    T_seq = 'TCTGCGTTGATACCACT'
    T_len = len(T_seq)
    distances = [(-(offset+T_len+10),hamming_distance(readstring[-(offset+T_len):-offset], T_seq)) for offset in range(7)]
    argmin   = [x[1] for x in distances].index(min([x[1] for x in distances]))
    return distances[argmin]

readtab2 = readtab2.assign(forward_distance = [find_adapter(x) for x in readtab2.seq])
readtab2 = readtab2.assign(reverse_distance = [find_adapter(x) for x in readtab2.rcseq])

readtab2 = readtab2.assign(direction = ['R' if x.reverse_distance[1]<x.forward_distance[1] else 'F' for i,x in readtab2.iterrows()])

readtab2 = readtab2.assign(min_distance = [x.reverse_distance[1] if x.direction =='R' else x.forward_distance[1] for i,x in readtab2.iterrows()])


# In[18]:


readtab2.min_distance.value_counts()


# In[19]:


readtab2 = readtab2.assign(barcode_pos = [x.reverse_distance[0] if x.direction =='R' else x.forward_distance[0] for i,x in readtab2.iterrows()])

readtab2 = readtab2[readtab2.min_distance ==0]
readtab2 = readtab2.assign(cseq = [x.seq if x.direction=='F' else x.rcseq for i,x in readtab2.iterrows()])

readtab2 = readtab2.drop(columns=['seq','rcseq','reverse_distance', 'forward_distance'])


# In[20]:


temp = [find_TSO(x) for x in readtab2.cseq]


# In[21]:


pd.Series(x[1] for x in temp).value_counts()


# In[22]:


readtab2 = readtab2.assign(st_barcode = [x.cseq[x.barcode_pos:x.barcode_pos+16] for i,x in readtab2.iterrows()])
readtab2 = readtab2.assign(st_umi = [x.cseq[x.barcode_pos+16:x.barcode_pos+28] for i,x in readtab2.iterrows()])


# In[23]:


readtab2 = readtab2.assign(st_polyT = [x.cseq[x.barcode_pos+28:x.barcode_pos+32] for i,x in readtab2.iterrows()])


# In[24]:


readtab2 = readtab2[(readtab2.st_umi!='TTTTTTTTTTTT') & (readtab2.st_polyT=='TTTT')]


# In[25]:


readtab2 = readtab2.reset_index()


# In[26]:


readtab2 = readtab2[[re.search('[^T]T{0,2}[^T]T{0,2}[^T]', readtab2.cseq[i][readtab2.barcode_pos[i]+28:])!=None for i in range(len(readtab2))]]


# In[27]:


readtab2 = readtab2.assign(polyT_end_pos = [re.search('[^T]T{0,2}[^T]T{0,2}[^T]', x.cseq[x.barcode_pos+28:]).span()[0]+x.barcode_pos+28 for i,x in readtab2.iterrows()])


# In[28]:


spottab = pickle.load(open('pickeld_results/visium-v1-coordiantes.pk','rb'))


# In[29]:


bset = set(spottab.index)

readtab2['y_coor'] = [spottab.loc[x,'y']-1 if x in bset else None for x in readtab2.st_barcode]
readtab2['x_coor'] = [spottab.loc[x,'x']-1 if x in bset else None for x in readtab2.st_barcode]


# In[30]:


os.makedirs(f'MIXCR_results/run86/IGH/{sample}', exist_ok=True)


# In[31]:


os.makedirs(f'MIXCR_results/run86/TCR/{sample}', exist_ok=True)


# In[32]:


readtab2['nqid'] = readtab2['qid']+'-'+readtab2['library']


# In[33]:


readtab2 = readtab2.set_index('nqid')


# In[112]:


with open('data/trimmed_fastq_files/'+sample+'-trimmed.fastq','w') as fq:
    for i,row in sheet.iterrows():
        fastq_file_1 = row['fastq_file']
        for record in SeqIO.parse(open(fastq_file_1,'rt'), 'fastq'):
            qid = record.id + '-' + row['Hybrid capture sample list']
            if qid in readtab2.index:
                record.id = record.id + '-' + row['Hybrid capture sample list']
                if readtab2.loc[qid,'direction']=='F':
                    SeqIO.write(record[readtab2.loc[qid,'polyT_end_pos']:-31], fq, 'fastq')
                if readtab2.loc[qid,'direction']=='R':
                    SeqIO.write(record.reverse_complement(id=True, name=True, description=True)[readtab2.loc[qid,'polyT_end_pos']:-31], fq, 'fastq')


# In[150]:


process = Popen(f'mixcr align  -f -s hs -O saveOriginalReads=true data/trimmed_fastq_files/{sample}-trimmed.fastq MIXCR_results/run86/IGH/{sample}/{sample}IGH.vdjca --report MIXCR_results/run86/IGH/{sample}/{sample}IGH-align.report', shell=True)
process.wait()


# In[151]:


process = Popen(f'mixcr assemble -f --write-alignments MIXCR_results/run86/IGH/{sample}/{sample}IGH.vdjca MIXCR_results/run86/IGH/{sample}/{sample}IGH.clna --report MIXCR_results/run86/IGH/{sample}/{sample}IGH-assemble.report', shell=True)
process.wait()


# In[152]:


process = Popen(f'mixcr  exportClones -f MIXCR_results/run86/IGH/{sample}/{sample}IGH.clna MIXCR_results/run86/IGH/{sample}/{sample}IGH.txt', shell=True)
process.wait()


# In[153]:


process = Popen(f'mixcr exportAlignments  -f -cloneIdWithMappingType -cloneId -readIds -descrsR1  MIXCR_results/run86/IGH/{sample}/{sample}IGH.clna MIXCR_results/run86/IGH/{sample}/{sample}IGH-alignments.txt', shell=True)
process.wait()


# In[154]:


alignmnet_report = f'MIXCR_results/run86/IGH/{sample}/{sample}IGH-alignments.txt'


# In[155]:


altable = pd.read_csv(alignmnet_report, sep = '\t')


# In[156]:


altable.descrsR1[0]


# In[157]:


altable.head()


# In[158]:


altable['qid'] = [x.split()[0] for x in altable.descrsR1]


altable = altable.set_index('qid').merge(readtab2, left_index=True, right_index=True)


# In[159]:


altable['coordinates'] = altable.x_coor.astype(str) + '-' + altable.y_coor.astype(str)


# In[160]:


altable.library


# In[161]:


altable[altable.cloneId>-1].groupby(['library','cloneId']).apply(lambda x: x['coordinates'].value_counts()).to_csv(f'tabular_results/run86-{sample}-clonotype_coordinates.csv')


# In[162]:


clonotypes_csv = f'MIXCR_results/run86/IGH/{sample}/{sample}IGH.txt'


# In[163]:


clonetable = pd.read_csv(clonotypes_csv, sep = '\t')


# In[164]:


clone_count_by_library = altable[altable.cloneId>-1].groupby(['cloneId','library']).agg(read_count = pd.NamedAgg(column='descrsR1',aggfunc='count'))


# In[165]:


pickle.dump(altable,open(f'pickeld_results/run86-{sample}-altable.pk','wb'))


# In[3]:


import pickle
altable = pickle.load(open(f'pickeld_results/run86-{sample}-altable.pk','rb'))


# In[6]:


altable.to_csv('tabular_results/{}-altable.tsv.gz'.format(sample), sep='\t', index=False, compression='gzip')


# In[166]:


for lib in sheet['Hybrid capture sample list'].unique():
    clonetable[lib+'_readcounts'] =  [clone_count_by_library.loc[(x,lib),'read_count'] if (x,lib) in clone_count_by_library.index else 0 for x in clonetable.cloneId]


# In[167]:


clonetable.to_csv(f'tabular_results/Pacbio_run86_{sample}-Mixcr_all.tsv', sep ='\t', index=False)


# In[168]:


umi_clone_table = altable[(altable.cloneId>-1) & (altable.coordinates != 'nan-nan')].groupby(
    ['coordinates','library','st_umi','x_coor','y_coor','st_barcode']).agg(NumberofClones = pd.NamedAgg(
    column='cloneId', aggfunc=lambda x:len(x.unique())),
    Clones = pd.NamedAgg(column='cloneId', aggfunc= lambda x:tuple(x.unique())),
    Readcount = pd.NamedAgg(column='cseq', aggfunc='count'),
    clone_read_counts = pd.NamedAgg(
    column='cloneId', aggfunc=lambda x:x.value_counts().to_dict()))


# In[169]:


rn_dict = umi_clone_table.groupby(['Readcount','NumberofClones']).apply(lambda x: len(x)).to_dict()
umi_clone_table['Log2Frequencyx+1'] = [np.log2(rn_dict[(x.Readcount,x.NumberofClones)]+1) for i,x in umi_clone_table.iterrows()]

px.scatter(data_frame=umi_clone_table, x='Readcount', y='NumberofClones', size='Log2Frequencyx+1', color = 'Log2Frequencyx+1', color_continuous_scale='viridis',color_continuous_midpoint=0)


# In[170]:


umi_clone_table_filtered_readcount = [umi_clone_table[(umi_clone_table.NumberofClones==1) & (umi_clone_table.Readcount>i)].reset_index() for i in range(4)]


# In[186]:


def get_clonetable_with_counts(df,library):
    test = df[df.library==library].copy()
    
    test.Clones = [x[0] for x in test.Clones]


    test.clone_read_counts = [x.clone_read_counts[x.Clones] for i,x in test.iterrows()]

    test2 = test.groupby('Clones').agg(read_count = pd.NamedAgg(column ='clone_read_counts', aggfunc = 'sum'),
                                umi_count = pd.NamedAgg(column = 'clone_read_counts', aggfunc = 'count')).reset_index()
    clonetab = pd.read_csv(clonotypes_csv, sep = '\t')
    return test2.rename(columns={'Clones':'cloneId'}).merge(clonetab)


# In[171]:


def get_umi_count_matrix(df, library):

    """gets the UMI clone tables and a library, filters out other libraries in the df and converts it to a UMI count Dataframe"""

    result = df[df.library==library]

    result.Clones = [x[0] for x in result.Clones]

    return result.groupby(['st_barcode','Clones']
                         ).agg({'st_umi':lambda x:x.nunique()}
                              ).reset_index().pivot(
        index='st_barcode', columns='Clones', values='st_umi').fillna(0)


# In[187]:


for i in range(4):
    for library in umi_clone_table_filtered_readcount[i].library.unique():
        get_umi_count_matrix(umi_clone_table_filtered_readcount[i], library).to_csv(f'tabular_results/run86-{library}-filtered-filter-readcount-{i}-umi-count.tsv',sep='\t')
        get_clonetable_with_counts(umi_clone_table_filtered_readcount[i], library).to_csv(f'tabular_results/run86-{library}-filtered-filter-readcount-{i}-clonotypes.tsv',sep='\t')


# In[174]:


import seaborn as sns

altable = altable.assign(seqlen = [len(x) for x in altable.cseq])

altable = altable.assign(hasClone = [x>-1 for x in altable.cloneId])


plot = sns.displot(altable, x="seqlen", hue="hasClone", kind="kde", fill=True, aspect=4,bw_adjust=.15)

plot.fig.suptitle(f'{sample} fragment length by whether the fragment is assigned to any clone')


plot.savefig(f"run86-{sample}-hasClone.pdf")


# In[175]:


alclonetable = altable[altable.cloneId>-1].merge(clonetable, left_on='cloneId' , right_on='cloneId')
alclonetable = alclonetable.assign(Gene = [y[:3]  if isinstance(y,str) else y for y in alclonetable.allCHitsWithScore])


# In[176]:





plot = sns.displot(alclonetable[~pd.isna(alclonetable.Gene)], x="seqlen", hue="Gene", kind="kde", fill=True, aspect=4,bw_adjust=.15)

plot.fig.suptitle(f'{sample} fragment length by aligned gene distribution')


plot.savefig(f"run86-{sample}.pdf")


# In[177]:


for gene in alclonetable.Gene.unique():
    if not pd.isna(gene):
        plot = sns.displot(alclonetable[alclonetable.Gene==gene], x="seqlen", kind="kde", fill=True, aspect=4,bw_adjust=.15)

        plot.fig.suptitle(f'{sample} fragment length distribution for gene {gene}')


        plot.savefig(f"run86-{sample}-{gene}.pdf")


# In[178]:


mixcr_reads = set(altable.index)

readtab2 = readtab2.assign(In_MIXCR = [x in mixcr_reads for x in readtab2.index])

readtab2 = readtab2.assign(seqlen = [len(x) for x in readtab2.cseq])


# In[179]:


plot = sns.displot(readtab2, x="seqlen", hue="In_MIXCR", kind="kde", fill=True, aspect=4,bw_adjust=.15)

plot.fig.suptitle(f'{sample} fragment length by whether the fragment is in MIXCR output')


plot.savefig(f"image_results/run86-{sample}-IN_MIXCR.pdf")


# In[180]:


with open('data/trimmed_fastq_files/'+sample+'-trimmed-background.fastq','w') as fq:
    for fastq_file_1 in sheet.fastq_file:
        for record in SeqIO.parse(open(fastq_file_1,'rt'), 'fastq'):
            qid = record.id
            if qid in readtab2.index:
                if readtab2.loc[qid,"In_MIXCR"]==False:
                    if readtab2.loc[qid,'direction']=='F':
                        SeqIO.write(record[readtab2.loc[qid,'polyT_end_pos']:-31], fq, 'fastq')
                    if readtab2.loc[qid,'direction']=='R':
                        SeqIO.write(record.reverse_complement(id=True, name=True, description=True)[readtab2.loc[qid,'polyT_end_pos']:-31], fq, 'fastq')


# In[181]:


p =Popen(f'minimap2 -ax splice:hq -uf /home/hosein.toosi/exome-analysis/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa  data/trimmed_fastq_files/{sample}-trimmed-background.fastq > data/Bam/run86-{sample}-background.sam',shell=True)


# In[182]:


p.wait()


# In[ ]:




