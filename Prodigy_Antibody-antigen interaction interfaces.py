#!/usr/bin/env python
# coding: utf-8

# # Import and Install

# In[2]:


# Import modules
import numpy as np
import pandas as pd
import os


# In[3]:


# Install Biopython and freesasa scripts
get_ipython().system('pip install biopython')
get_ipython().system('pip install freesasa')


# In[4]:


# Install PRODIGY Software
get_ipython().system('git clone http://github.com/haddocking/prodigy')
get_ipython().system('cd prodigy')
get_ipython().system('pip install prodigy')


# In[5]:


get_ipython().system('prodigy --help')


# In[6]:


get_ipython().system('prodigy "C:\\Users\\lucia\\OneDrive\\Desktop\\Python\\AMPM\\filtered_pdb\\1kb5_HLA_B\\com.pdb" --selection H,L A,B')


# In[7]:


# Import Biopython modules
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import seq3
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


# # Read data

# In[8]:


reference = pd.read_csv("data_team.csv")
reference


# In[9]:


reference.info()


# In[10]:


sample = pd.read_csv("CDRdata_team1.csv")
sample


# In[11]:


sample.info()


# # Using “CDRdata_team1.csv” to map back CDR (H1, H2, H3, L1, L2, L3) sequences with respect to pdb-ID

# In[12]:


ref_lst = reference['pdb'].values.tolist()
filtered_data = sample[sample['PDB'].apply(lambda x: str(x).lower() in ref_lst)]
filtered_data = filtered_data.reset_index()
filtered_data


# # Calculate length

# In[13]:


filtered_data['length'] = filtered_data['seq'].str.len()
filtered_data = filtered_data.reset_index(drop=True) 
filtered_data


# # Count the number of each aa type in each of the CDR sequences

# In[14]:


count = []
i = 0
while i < len(filtered_data) : 
    x = ProteinAnalysis(filtered_data.loc[i,'seq'])
    aminoacidlist = x.count_amino_acids()
    pdb_id = filtered_data.iloc[i,1]
    aminoacidlist['PDB_ID'] = pdb_id
    count.append(aminoacidlist)
    i += 1
count = pd.DataFrame(count)
count


# In[15]:


# Change Column name
count.rename(columns={'A': seq3('A'), 'C': seq3('C'), 'D': seq3('D'), 'E': seq3('E'), 'F': seq3('F'), 'G': seq3('G'), 'H': seq3('H'), 'I': seq3('I'), 'K': seq3('K'), 'L': seq3('L'), 'M': seq3('M'), 'N': seq3('N'), 'P': seq3('P'), 'Q': seq3('Q'), 'R': seq3('R'), 'S': seq3('S'), 'T': seq3('T'), 'V': seq3('V'), 'W': seq3('W'), 'Y': seq3('Y')}, inplace = True)
count.info()


# In[16]:


count.columns
count


# In[17]:


result = pd.concat([filtered_data, count],axis=1,sort=True)
result


# In[18]:


result.to_csv("result.csv", index=False)


# In[19]:


ref_2 = result['PDB'].values.tolist()
filtered_prodigy = reference[reference['pdb'].apply(lambda x: str(x).upper() in ref_2)]
filtered_prodigy = filtered_prodigy.reset_index()
filtered_prodigy


# In[20]:


filtered_prodigy.info()


# In[21]:


# Change " | " to "_" in PDB_CHAIN_ANTI
filtered_prodigy['PDB_CHAIN_ANTI']= filtered_prodigy['PDB_CHAIN_ANTI'].apply(lambda x: x.replace(" | ","_"))
filtered_prodigy


# In[22]:


# List of all folders:
path = "filtered_pdb"
all_folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
#os.listdir to get a list of all files and folders 
#os.path.isdir to filter the list and only include folders
len(all_folders)


# In[23]:


# List of folders to find
find_folders = filtered_prodigy['PDB_CHAIN_ANTI'].values.tolist()
len(find_folders)



# In[24]:


filtered_result = [x for x in all_folders if x in find_folders]
len(filtered_result)


# In[25]:


# Change " | " to " , " in antigen_chain
filtered_prodigy['antigen_chain'] = filtered_prodigy['antigen_chain'].apply(lambda x: x.replace(" | ",","))
filtered_prodigy['antigen_chain']


# In[35]:


# Prodigy function
def result_prodigy(i,error_pdb,error_runprodigy):
    pdbname = filtered_prodigy.loc[i,'PDB_CHAIN_ANTI']
    pdbfile = "filtered_pdb/" + pdbname 
    print(pdbfile)
    hchain = filtered_prodigy.loc[i,'Hchain']
    lchain = filtered_prodigy.loc[i,'Lchain']
    antigen = filtered_prodigy.loc[i,'antigen_chain']
    try : 
        if hchain == lchain :
            error_pdb.append(pdbname)
            get_ipython().system('prodigy "$pdbfile/com.pdb" --selection $hchain $antigen --contact_list > "$pdbfile/output.txt"')
        else :
            get_ipython().system('prodigy "$pdbfile/com.pdb" --selection $hchain,$lchain $antigen --contact_list > "$pdbfile/output.txt"')
    except :
        error_runprodigy.append(pdbname)


# In[36]:


# Run prodigy for whole data:
i = 0
error_pbd = []
error_runprodigy = []
while i < len(filtered_result) : 
    result_prodigy(i,error_pdb,error_runprodigy)
    i += 1


# In[30]:


get_ipython().system('prodigy "filtered_pdb/7l1v_HHA_B/com.pdb" --selection H A,B --contact_list > output.txt')


# In[39]:


pd.DataFrame(error_pdb)


# In[40]:


print(error_runprodigy)


# In[ ]:




