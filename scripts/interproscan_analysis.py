# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 00:51:42 2024

@author: Kausthubh R
"""

import os
from Bio import SeqIO
from pandas import read_csv, concat
import pandas as pd


# getting query sequences

print("Please ensure that this code is run with the folder containing your bbh run as a starting point\n\n")
og_seqs = []
os.chdir("query_sequences")
temp = os.listdir()
for i in temp:    
    og_seqs.append([x.id for x in SeqIO.parse(i, 'fasta')][0])
del i, temp

os.chdir(os.path.dirname(os.getcwd()))

dirname = [x for x in os.listdir() if x.endswith("new_species")]
if not dirname:
    dirname = [x for x in os.listdir() if x.endswith("final")]
    
os.chdir(dirname[0])
base_list = os.listdir()
base_path = os.getcwd()

interpro_files = {}
sequences = {}

for i in base_list:
    os.chdir(i)
    os.chdir("interproscan")
    interpro_files[i] = os.path.join(os.getcwd(), [x for x in os.listdir() if x.endswith("tsv")][0])
    sequences[i] = [x for x in SeqIO.parse([x for x in os.listdir() if x.endswith("fasta")][0], 'fasta')]
    os.chdir(base_path)
del i

interpro_files_read = {x : read_csv(interpro_files[x], sep="\t", index_col=None, names=["protein ID", "seq MD5 digest", "seq length", "analysis", "sign accession", "sign desc", "start", "stop", "score", "status", "date", "interpro accession", "interpro desc", "go annot", "pathway annot"]) for x in interpro_files}

#%%

selected_seqs = {}
selected_seqs_write = {}
# log_text = {}
list_of_analysis = list(dict.fromkeys(concat([interpro_files_read[x]["analysis"] for x in interpro_files_read])))
analysis = "SUPERFAMILY"
alternate_analysis = "PANTHER"

for i in interpro_files_read:
    # input("you should probably wait and see")
    print(i)
    flag_sf = 1
    flag_interpro = 1
    ref_prot = [x for x in [y.id for y in sequences[i]] if any(z in x for z in og_seqs)]
    temp = interpro_files_read[i]
    ref_analysis = temp[temp["protein ID"].isin(ref_prot)]
    ref_analysis_sf = ref_analysis[ref_analysis["analysis"].isin([analysis])]
    ref_analysis_pan = ref_analysis[ref_analysis["analysis"].isin([alternate_analysis])]
    
    if len(ref_analysis_sf) == 0:
        flag_sf = 0
    
    if any([sorted(list(set(ref_analysis_pan[ref_analysis_pan["protein ID"] == x]["interpro accession"]))) == ["-"] for x in ref_prot]):
        flag_interpro = 0
    # if len(ref_analysis_pan) == 0:
        # flag_pan = 0
        
    prot_list = list(set(temp["protein ID"]))
    prot_list_selected = []
    selected_data = []
    # if i == "Orthogroup_16":
    #     input("hi")
    for j in prot_list:
        # input()
        # if j.startswith("Ihof_5892"):
        #     input("hi")
        temp_analysis = temp[temp["protein ID"].isin([j])]
        temp_sf = temp_analysis[temp_analysis["analysis"].isin([analysis])]
        temp_pan = temp_analysis[temp_analysis["analysis"].isin([alternate_analysis])]
        
        if flag_sf == 1 and flag_interpro == 1:
            if len(list(temp_sf["sign accession"])) > max([len(list(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["sign accession"])) for x in ref_prot]):
                continue       
        
            if not any([sorted(list(set(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["sign accession"]))) == sorted(list(set(temp_sf["sign accession"]))) for x in ref_prot]):
                if not any([sorted(list(set(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["interpro accession"]))) == sorted(list(set(temp_sf["interpro accession"]))) for x in ref_prot]):
                    if not any([sorted(list(set(ref_analysis_pan[ref_analysis_pan["protein ID"] == x]["interpro accession"]))) == sorted(list(set(temp_pan["interpro accession"]))) for x in ref_prot]):
                        continue
        
        if flag_sf == 1 and flag_interpro != 1:
            if len(list(temp_sf["sign accession"])) > max([len(list(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["sign accession"])) for x in ref_prot]):
                continue       
        
            if not any([sorted(list(set(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["sign accession"]))) == sorted(list(set(temp_sf["sign accession"]))) for x in ref_prot]):                
                continue
        
        pan_score = 0
        sf_score = 0
        
        if any([sorted(list(set(ref_analysis_pan[ref_analysis_pan["protein ID"] == x]["sign accession"]))) == sorted(list(set(temp_pan["sign accession"]))) for x in ref_prot]):
            pan_score = 1
            
        if pan_score == 0:
            if any([sorted(list(set(ref_analysis_pan[ref_analysis_pan["protein ID"] == x]["interpro accession"]))) == sorted(list(set(temp_pan["interpro accession"]))) for x in ref_prot]):
                pan_score = 1
        
        if any([sorted(list(set(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["sign accession"]))) == sorted(list(set(temp_sf["sign accession"]))) for x in ref_prot]):
            if all([len(list(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["sign accession"])) == len(list(temp_sf["sign accession"])) for x in ref_prot]):
                sf_score = 1
            
        if sf_score == 0:
            if any([sorted(list(set(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["interpro accession"]))) == sorted(list(set(temp_sf["interpro accession"]))) for x in ref_prot]):
                if all([len(list(ref_analysis_sf[ref_analysis_sf["protein ID"] == x]["interpro accession"])) == len(list(temp_sf["interpro accession"])) for x in ref_prot]):
                    sf_score = 1
                
        if any([x == 1 for x in [pan_score, sf_score]]):
            prot_list_selected.append(j)
            selected_data.append(pd.concat([temp_sf, temp_pan]))
        
    selected_seqs[i] = pd.concat(selected_data)
    selected_seqs_write[i] = [x for x in sequences[i] if x.id in prot_list_selected]
    # selected_seqs_write[i] = [x for x in prot_list_selected]
    # input("you should probably wait and see")
        
del i, j, flag_sf, pan_score, sf_score, ref_analysis, ref_analysis_pan, ref_analysis_sf, ref_prot, temp, temp_analysis, temp_pan, temp_sf, selected_data, prot_list, prot_list_selected, flag_interpro
    

# del i, j, temp, temp1, ref_prot, ref_prot_data, temp_df, max_count, analysis_used, prot_list, prot_list_filter_length, ref_analysis, list_of_analysis

#%%
# temp1 = []
# for i in selected_seqs:
#     ref_prot = [x for x in [y.id for y in sequences[i]] if any(z in x for z in og_seqs)]
#     temp = interpro_files_read[i]
#     ref_analysis = temp[temp["protein ID"].isin(ref_prot)]
#     ref_analysis_sf = ref_analysis[ref_analysis["analysis"].isin([analysis])]
#     ref_analysis_pan = ref_analysis[ref_analysis["analysis"].isin([alternate_analysis])]
    
#     temp1.append(ref_analysis_sf)
#     temp1.append(ref_analysis_pan)
# del i, ref_analysis, ref_analysis_sf, ref_prot, ref_analysis_pan, temp

# temp1 = pd.concat(temp1)
    
# temp1.to_csv("reference_protein_annotations.csv")
#%%
# os.chdir(os.path.dirname(os.getcwd()))
# os.mkdir("orthogroups_post_intepro_analysis_3")
# os.chdir("orthogroups_post_intepro_analysis_3")
for i in selected_seqs_write:
    # os.mkdir(i)
    os.chdir(i)
    with open(f"{i}_post_interpro.fasta", 'w') as f:
        [SeqIO.write(x, f, 'fasta') for x in selected_seqs_write[i]]
        f.close()
    selected_seqs[i].to_csv(f"{i}_intepro_analysis.csv")
    os.chdir(os.path.dirname(os.getcwd()))
del i, f

#%%

