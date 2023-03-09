#!/usr/bin/env python
# coding: utf-8

# # Import and Install

# Import modules
import numpy as np
import pandas as pd
import os # Manipulate files and folders
import re # Extract positive/negative floating-point numbers from a string

# Install Biopython and freesasa scripts
!pip install biopython
!pip install freesasa

# Install PRODIGY Repository
!git clone http://github.com/haddocking/prodigy
!cd prodigy
!pip install prodigy

!prodigy --help

# Test sample
!prodigy "filtered_pdb/1kb5_HLA_B/com.pdb" --selection H,L A,B

# Import Biopython modules
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import seq3
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# # Read data
reference = pd.read_csv("data_team.csv")
reference

reference.info()

sample = pd.read_csv("CDRdata_team1.csv")
sample

sample.info()

# # Using 'CDRdata_team1.csv' to map back CDR (H1, H2, H3, L1, L2, L3) sequences with respect to pdb-ID

ref_lst = reference['pdb'].values.tolist()
filtered_data = sample[sample['PDB'].apply(lambda x: str(x).lower() in ref_lst)]
filtered_data = filtered_data.reset_index()
filtered_data

# # Calculate length

filtered_data['length'] = filtered_data['seq'].str.len()
filtered_data = filtered_data.reset_index(drop=True) 
filtered_data

# # Count frequency of amino acids type in each of the CDR sequences

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

# Change Column name : convert one-letter amino-acid code into the three-letter one
count.rename(columns={'A': seq3('A'), 'C': seq3('C'), 'D': seq3('D'), 'E': seq3('E'), 'F': seq3('F'), 'G': seq3('G'), 'H': seq3('H'), 'I': seq3('I'), 'K': seq3('K'), 'L': seq3('L'), 'M': seq3('M'), 'N': seq3('N'), 'P': seq3('P'), 'Q': seq3('Q'), 'R': seq3('R'), 'S': seq3('S'), 'T': seq3('T'), 'V': seq3('V'), 'W': seq3('W'), 'Y': seq3('Y')}, inplace = True)
count.info()

count.columns
count

result = pd.concat([filtered_data, count],axis=1,sort=True) # Map back amino-acid-counting DataFrame to 'result' DataFrame
result

# Edit result DataFrame
result = result.drop('index', axis=1) # Delete 'index' column
result = result.drop('PDB_ID', axis=1) # Delete 'PDB_ID' column
result.to_csv("result.csv", index=False) # Save 'result' DataFrame to csv file
result

# # Predict binding affinity using PRODIGY

# ### Prepare data

ref_2 = result['PDB'].values.tolist()
filtered_prodigy = reference[reference['pdb'].apply(lambda x: str(x).upper() in ref_2)]
filtered_prodigy = filtered_prodigy.reset_index()
filtered_prodigy

filtered_prodigy.info()

# Change " | " to "_" in PDB_CHAIN_ANTI
filtered_prodigy['PDB_CHAIN_ANTI']= filtered_prodigy['PDB_CHAIN_ANTI'].apply(lambda x: x.replace(" | ","_"))
filtered_prodigy

# List of all folders
path = "filtered_pdb"
#os.listdir to get a list of all files and folders 
#os.path.isdir to filter the list and only include folders
all_folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
len(all_folders)

# List of finding folders 
find_folders = filtered_prodigy['PDB_CHAIN_ANTI'].values.tolist()
len(find_folders)

# How many finding folders are there in available folders?
filtered_result = [x for x in all_folders if x in find_folders]
len(filtered_result)

# Change " | " to " , " in antigen_chain
filtered_prodigy['antigen_chain'] = filtered_prodigy['antigen_chain'].apply(lambda x: x.replace(" | ",","))
filtered_prodigy['antigen_chain']

# ### Run PRODIGY 

# Build a function of PRODIGY
def result_prodigy(i,error_pdb,error_runprodigy):
    pdbname = filtered_prodigy.loc[i,'PDB_CHAIN_ANTI']
    pdbfile = "filtered_pdb/" + pdbname 
    print(pdbfile)
    hchain = filtered_prodigy.loc[i,'Hchain']
    lchain = filtered_prodigy.loc[i,'Lchain']
    antigen = filtered_prodigy.loc[i,'antigen_chain']
    try : # test a block of code for errors
        if hchain == lchain :
            error_pdb.append(pdbname)
            !prodigy "$pdbfile/com.pdb" --selection $hchain $antigen --contact_list > "$pdbfile/output.txt"
        else :
            !prodigy "$pdbfile/com.pdb" --selection $hchain,$lchain $antigen --contact_list > "$pdbfile/output.txt"
    except :
        error_runprodigy.append(pdbname)

# Run PRODIGY for whole data
i = 0
error_pdb = [] # Errors due to pdb problems
error_runprodigy = [] # Errors due to PRODIGY-code problems
while i < len(filtered_prodigy) : 
    result_prodigy(i,error_pdb,error_runprodigy)
    i += 1
print(i) # Find out how many samples are covered?

pd.DataFrame(error_pdb) # Print pdbs' errors

error_runprodigy # Print prodigys' errors

# ### Export output

# Build a function of extraction information "Predicted binding affinity" from output.txt
def extract_output(i):
    pdbname = filtered_prodigy.loc[i,'PDB_CHAIN_ANTI']
    pdbfile = "filtered_pdb/" + pdbname + "/output.txt"
    file_exists = os.path.isfile(pdbfile)
    if file_exists:
        with open(pdbfile, mode = "r") as f:
            # "Predicted binding affinity (kcal.mol-1)" is the last line in file.txt
            lastline = f.readlines()[-1] # read the last line
            dGprodigy = re.findall(r'[-+]?\d*\.\d+', lastline) # class of output is list, eg: ['-9.9']
            dGprodigy = ''.join(dGprodigy) # convert to class 'str'
            filtered_prodigy.at[i, 'dG_prodigy'] = dGprodigy # insert results into 'dGprodigy' column in DataFrame
            f.close() # close the file
    else:
        print(pdbname)

# Extract information for whole data
for i in range(len(filtered_prodigy)):
    extract_output(i)

filtered_prodigy

filtered_prodigy.info()

filtered_prodigy.to_csv("filtered_prodigy.csv", index=False) # Save results to csv file 

# # Final step: Map back PRODIGY's results to 'result' DataFrame based on PDB and original_chain
filtered_prodigy['pdb']=filtered_prodigy['pdb'].str.upper() # Uppercase 'pdb' column in 'filtered_prodiy' DataFrame to compare to 'PDB' column in 'result' DataFrame
filtered_prodigy

for i in range(len(filtered_prodigy)):
    for j in range(len(result)):
        if result.loc[j,'PDB'] == filtered_prodigy.loc[i,'pdb'] :
            if result.loc[j, 'original_chain'] == filtered_prodigy.loc[i,'Hchain'] :
                result.loc[j,'dG_prodigy'] = filtered_prodigy.loc[i,'dG_prodigy']
            elif result.loc[j, 'original_chain'] == filtered_prodigy.loc[i,'Lchain'] :
                result.loc[j,'dG_prodigy'] = filtered_prodigy.loc[i, 'dG_prodigy']   
      
result.insert(4, 'dG_prodigy',result.pop('dG_prodigy')) # Change the order of 'dG_prodigy' columns
result.to_csv('finalresult.csv', index=False) # Save results to csv file
result

# Check for NaN in 'result' DataFrame
result[result.isnull().any(axis=1)]

