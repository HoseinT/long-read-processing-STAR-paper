#!/usr/bin/env python
# coding: utf-8

# In[1]:


sample ='BCSA3'


# In[2]:


import pickle


# In[3]:


# load altable from the pickle file
altable  = pickle.load(open(f'pickeld_results/run86-{sample}-altable.pk','rb'))


# In[4]:


altable = altable[altable.cloneId>-1]


# In[5]:


## get the clonotype table output from MIXCR run
clonotypes_csv = f'MIXCR_results/run86/IGH/{sample}/{sample}IGH.txt'


# In[6]:


import pandas as pd
clonetable = pd.read_csv(clonotypes_csv, sep = '\t')


# In[7]:


clonetable


# In[8]:


altable = altable.assign(combined_barcode_umi = altable.st_barcode +'-' + altable.st_umi )


# In[9]:


altable = altable[altable.coordinates != 'nan-nan']


# In[10]:


## get a table with read counts of each clone in each library
clone_count_by_library = altable.groupby(['cloneId','library']).agg(read_count = pd.NamedAgg(column='descrsR1',aggfunc='count'),
                                                                                       umi_count = pd.NamedAgg(column='combined_barcode_umi',aggfunc='nunique'))


# In[11]:


import re


# In[12]:


altable = altable.assign(library2 = [re.sub('_\dT:\dB','',x) for x in altable.library])


# In[13]:


# filters for reads assigned to  a clone that also have a valid sport barcode, groups on the same UMI and collects info abut read count and clone count
umi_clone_table = altable.groupby(
    ['coordinates','library2','st_umi','x_coor','y_coor','st_barcode']).agg(NumberofClones = pd.NamedAgg(
    column='cloneId', aggfunc=lambda x:len(x.unique())),
    Clones = pd.NamedAgg(column='cloneId', aggfunc= lambda x:tuple(x.unique())),
    Readcount = pd.NamedAgg(column='cseq', aggfunc='count'),
    clone_read_counts = pd.NamedAgg(
    column='cloneId', aggfunc=lambda x:x.value_counts().to_dict()))


# In[14]:


umi_clone_table = umi_clone_table.reset_index()


# In[15]:


umi_clone_table = umi_clone_table.assign(combined_barcode_umi = umi_clone_table.st_barcode +'-' + umi_clone_table.st_umi)


# In[16]:


umi_clone_table.head()


# In[17]:


umis_to_delete = set(umi_clone_table[umi_clone_table.NumberofClones>1].combined_barcode_umi)


# In[20]:


(umi_clone_table.library2+'-'+umi_clone_table.coordinates+'-'+umi_clone_table.st_umi).nunique() , umi_clone_table.combined_barcode_umi.nunique()


# In[22]:


len(umis_to_delete)


# In[23]:


altable[~altable.combined_barcode_umi.isin(umis_to_delete)].to_csv(f'tabular_results/{sample}-after-umi-filteration.tsv',sep='\t',compression='gzip')


# In[ ]:




