# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 19:44:58 2024

@author: Kausthubh R
"""

import os
# =============================================================================
# Needed for KMeans as it is to have memory leakege on Windows with MKL when 
# there are less chunks than available threads.

# This environment variable initialization is needed BEFORE sklearn is imported.
# This is necessary as the variables are checked only once at first module 
# import.

nthreads = 1
os.environ["OMP_NUM_THREADS"] = str(nthreads) 
os.environ["OPENBLAS_NUM_THREADS"] = str(nthreads) 
os.environ["MKL_NUM_THREADS"] = str(nthreads)
del nthreads

# =============================================================================

import pandas as pd
from Bio import SeqIO
import jenkspy as jpy
from math import dist
from collections import OrderedDict
import statistics as stat
# from sklearn.cluster import KMeans
import numpy as np
from datetime import datetime


if os.name == 'nt':
    def winpath2linuxpath(obj1:str):
        temp = obj1
        temp = temp.replace("\\", "/")
        temp = "/mnt/" + temp[0].lower() + "/" + temp[3:]
        return temp
    
    def serverpath2winpath(obj1:str):
        temp = obj1
        temp = r"Z:" + temp[6:]
        temp = temp.replace("/", r"\\")
        return temp
    
    def winpath2serverpath(obj1:str):
        temp = obj1
        temp = temp.replace("Z:", "/g/dey")
        temp = temp.replace("\\", "/")
        temp = temp.replace("//", "/")
        return temp

def getlistofspeciesindb(obj1:str):
    listsp = []
    spdict = []
    for record in SeqIO.parse(obj1, 'fasta'):
        x = record.description
        y = record.id
        ind1 = x.rfind("]")
        ind2 = x.rfind("[")
        if(ind1<10):
            continue
        else:
            listsp.append(x[(ind2+1):(ind1)])
            spdict.append(y[:y.find("_")])
    # print(len(listsp))
    # print(len(spdict))
    spdict = OrderedDict.fromkeys(zip(spdict, listsp))
    listsp = list(dict.fromkeys(listsp))
    spdict = dict(OrderedDict(spdict))
    spdict = dict([(i[0], i[1]) for i in spdict])
    return listsp, spdict

def estimate_len_covered(obj2:list):
    # hmm_cov = sorted(zip(list(obj1)))
    # hmm_cov = sorted(hmm_cov)
    obj1 = sorted(obj2)

    check0 = [x for x in range(len(obj1)) if x > 0 and obj1[x-1][1] > obj1[x][0]]
    # if stat.stdev([x[1] for x in obj1])/stat.mean([x[1] for x in obj1]) < 0.1:
    repeats = len([x for x in check0 if (max(obj1[x-1]) - min(obj1[x]))/(max(obj1[x-1]) - min(obj1[x-1])) > 0.4])
    # else:
        # repeats
    check1 = [x for x in check0 if obj1[x][1] < obj1[x-1][1]]
    hmm_cov = [obj1[x] for x in range(len(obj1)) if x not in check1]
    check0 = [x for x in range(len(hmm_cov)) if x > 0 and hmm_cov[x-1][1] > hmm_cov[x][0]]
    while check0:    
        # for i in check0:
        i=check0[0]
        hmm_cov = sorted([hmm_cov[x] for x in range(len(hmm_cov)) if x != i and x != i-1] + [(hmm_cov[i-1][0], hmm_cov[i][1])])
        check0 = [x for x in range(len(hmm_cov)) if x > 0 and hmm_cov[x-1][1] > hmm_cov[x][0]]

    return sum([x[1] - x[0] + 1 for x in hmm_cov]), repeats


def coverage_calc(obj1:list, obj2:pd.DataFrame(), obj3:int):
    cov_df = {}
    
    for j in obj1:
        temp_df = obj2[obj2["target name"] == j]
        hmm_cov, repeats = estimate_len_covered(sorted(zip(list(temp_df["hmm coord from"]), list(temp_df["hmm coord to"]))))
        hmm_cov = round(hmm_cov/obj3, 2)
        
        
        if len(temp_df) == 0:
            cov_df[j] = [] 
            continue
        
        target_cov = round(estimate_len_covered(sorted(zip(list(temp_df["ali coord from"]), list(temp_df["ali coord to"]))))[0]/stat.mean(temp_df["target length"]),2)
        
        cov_df[j] = [hmm_cov, target_cov, repeats, stat.mean(temp_df["target length"])]
    return cov_df

#%%
assigned_ogs_updated = {"Orthogroup_0":["AGRV1_HUMAN"],
               "Orthogroup_1":["ANS4B_HUMAN","USH1G_HUMAN"],
               "Orthogroup_2":["BI2L2_HUMAN","BAIP2_HUMAN"],
               "Orthogroup_3":["CAD23_HUMAN","CDHR2_HUMAN","PCD15_HUMAN","CDHR5_HUMAN"],               
               "Orthogroup_4":["PLSI_HUMAN"],
               "Orthogroup_5":["CIB2_HUMAN"],
               "Orthogroup_6":["CLRN2_HUMAN"],
               "Orthogroup_7":["COBL_HUMAN"],
               "Orthogroup_8":["ES8L2_HUMAN","EPS8_HUMAN"],
               "Orthogroup_9":["ESPN_HUMAN"],
               "Orthogroup_10":["EZRI_HUMAN","MOES_HUMAN","RADI_HUMAN"],
               "Orthogroup_11":["FSCN2_HUMAN"],
               "Orthogroup_12":["GRCR1_HUMAN"],
               "Orthogroup_13":["MYH14_HUMAN","MYO1C_HUMAN","MYO3A_HUMAN","MYO3B_HUMAN","MYO6_HUMAN","MYO10_HUMAN","MYO15_HUMAN","MYO7A_HUMAN"],
               "Orthogroup_14":["NHRF1_HUMAN"],
               "Orthogroup_15":["PACN2_HUMAN"],
               "Orthogroup_16":["PDZD7_HUMAN","USH1C_HUMAN","WHRN_HUMAN"],
               "Orthogroup_17":["PKHG6_HUMAN"],
               "Orthogroup_18":["S100P_HUMAN"],
               "Orthogroup_19":["STRC_HUMAN"],
               "Orthogroup_20":["TMC1_HUMAN","TMC2_HUMAN"],
               "Orthogroup_21":["TWF2_HUMAN"],
               "Orthogroup_22":["USH2A_HUMAN"],
               "Orthogroup_23":["VILI_HUMAN"]}

path2db = r"/g/dey/Comparativegenomics/Mylan/Opis_Choano_92/opis_choanos.fasta"
path2db = r"/g/dey/Comparativegenomics/Mylan/Species_database/new_species_added_db_18032024.fsa"
if os.name == "nt":
    path2db = serverpath2winpath(path2db)
    listofspecies, speciesdict = getlistofspeciesindb(path2db)
else:
    listofspecies, speciesdict = getlistofspeciesindb(path2db)
#%%
base_list = os.listdir()
base_path = os.getcwd()

results = {}
results_dom = {}

for i in base_list:
    os.chdir(i)
    # check = max([int(x[-1]) for x in os.listdir() if x.endswith(".csv")])
    # os.chdir(f"hmm_searches_{check}")
    # results_hmm_len[f"{i}_{check}"] = hmm_len
    temp = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch.csv")][0])
    temp1 = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch_doms.csv")][0])
    del temp["Unnamed: 0"]
    del temp1["Unnamed: 0"]
    temp1 = temp1.astype(object)
    temp.fillna("None", inplace=True)
    temp1.fillna("None", inplace=True)
    results[f"{i}"] = temp
    results_dom[f"{i}"] = temp1
    
    os.chdir(base_path)
del i, temp, temp1, 

results_existing = {}
results_dom_existing = {}

os.chdir(os.path.dirname(os.getcwd()))
os.chdir([x for x in os.listdir() if x.startswith("verified")][0])
temp_path = os.getcwd()
for i in base_list:
    os.chdir(i)
    os.chdir("HMM_searches")
    check = max([int(x[-1]) for x in os.listdir()])
    os.chdir(f"hmm_searches_{check}")
    # results_hmm_len[f"{i}_{check}"] = hmm_len
    temp = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch.csv")][0])
    temp1 = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch_doms.csv")][0])
    del temp["Unnamed: 0"]
    del temp1["Unnamed: 0"]
    temp1 = temp1.astype(object)
    temp.fillna("None", inplace=True)
    temp1.fillna("None", inplace=True)
    results_existing[f"{i}"] = temp
    results_dom_existing[f"{i}"] = temp1
    
    os.chdir(temp_path)
del i, temp, temp1, check, temp_path

os.chdir(base_path)
#%%
for i in range(len(speciesdict.keys())):
    print(f"{i} - {list(speciesdict.keys())[i]} ({speciesdict[list(speciesdict.keys())[i]]})\n")
del i
existing_sp = input("\nWhile adding new species to an existing analysis, some species (preferably model systems) from the existing analysis should have been included as a reference. Please select the species abbreviations for those species by entering them in this format - 1, 2, 3:\n").split(", ")

existing_sp = [list(speciesdict.keys())[x] for x in range(len(list(speciesdict.keys()))) if str(x) in existing_sp]

os.chdir(os.path.dirname(os.getcwd()))
os.chdir([x for x in os.listdir() if x.startswith("verified")][0])
temp_base = os.getcwd()
temp_list = os.listdir()
og_seqs = {}
for i in temp_list:
    os.chdir(i)
    og_seqs[i] = [x.id for x in SeqIO.parse([y for y in os.listdir() if y.endswith("_updated.fasta")][0], 'fasta')]
    # input()
    # og_seqs[i] = [x for x in og_seqs[i] if any([x.startswith(z) for z in existing_sp])]
    os.chdir(temp_base)
del i, temp_base, temp_list
os.chdir(base_path)
#%%
final_seqs_for_tree = dict([(x, []) for x in assigned_ogs_updated])
added_seqs = []
thresholds_existing = {}
coverage_existing = {}
og_seqs_missing_in_hmmsearch = {}
threshold_used_new = {}
coverage_stats_new = {} 
for i in assigned_ogs_updated:
    to_add = []
    missing = {}
    # threshols = {}
    for j in assigned_ogs_updated[i]:
        hmm_len = stat.mean([int(x) for x in results_dom_existing[j]["query length"] if x != "None"])
        coverage = coverage_calc(og_seqs[j], results_dom_existing[j], hmm_len)
        coverage_existing[j] = coverage
        to_remove = [x for x in coverage if not coverage[x]]
        
        if to_remove:
            missing[j] = to_remove
        coverage = {x: y for x,y in coverage.items() if y}    
        min_hmm_cov = [x[0] for x in coverage.values()]
        min_target_cov = [x[1] for x in coverage.values()]
        length = [x[3] for x in coverage.values()]
            # hmm_cov_stdev = round(stat.stdev(min_hmm_cov), 3)
            # target_cov_stdev = round(stat.stdev(min_target_cov), 3)
        thresholds_existing[j] = [(max(min_hmm_cov), min(min_hmm_cov)), (max(min_target_cov), min(min_target_cov)), (max(length), min(length))]
        [to_add.append(x) for x in og_seqs[j]]
    final_seqs_for_tree[i] = to_add
    added_seqs += to_add
    # thresholds_used[i] = threshols
    og_seqs_missing_in_hmmsearch[i] = missing
del i, min_hmm_cov, min_target_cov, length, coverage, to_remove, to_add, j, missing
#
jnb = jpy.JenksNaturalBreaks(2)
for i in results:
    print(f"\n{i}")
    threshold = {}
    covs = {}
    check = [x for x in assigned_ogs_updated.keys() if i in assigned_ogs_updated[x]][0]
    rejections = pd.DataFrame(columns=["Sequence ID", "Reason"])
    rejections1 = {x:[] for x in speciesdict}
    ind1 = len(rejections)
    for j in speciesdict:
        # if j in existing_sp:
        #     continue
        # if j=="Hsap":
        #     print("hi")
        #     input()
    
        temp = results[i][results[i]["target name"].str.startswith(j)]
        temp1 = list(map(sum, zip(list(temp["full sequence score"]), list(temp["best 1 domain score"]))))
        
        if len(list(dict.fromkeys(temp1))) > 2:
            jnb.fit(temp1)
            threshold[f"bitscore_threshold_{j}"] = [max(jnb.groups_[0])]
            temp_ind = np.add(temp["full sequence score"], temp["best 1 domain score"]) > max(jnb.groups_[0])
            temp = [temp["target name"][x] for x in temp_ind.index if temp_ind[x]]
            
        else:
            temp = list(temp["target name"])
        # temp = [x for x in temp if x not in added_seqs]
        if not temp:
            rejections.loc[ind1] = [[], f"No sequence found in {j}"]
            ind1 += 1
            rejections1[j] += []
            covs[j] = []
            # threshold[j]
            continue
        temp_cov = coverage_calc(temp, results_dom[i], hmm_len)
        to_remove = [x for x in temp_cov if not temp_cov[x] and x not in added_seqs]
        rejections.loc[ind1] = [to_remove, f"Sequences were not found in the hmmsearch for {j}"]
        ind1 += 1 
        rejections1[j] += to_remove
        
        coverage = {x: y for x,y in temp_cov.items() if y} 
        
        covs[j] = coverage
        temp = list(coverage.keys())
        
        for k in temp:
            if k in added_seqs or k in rejections1[j]:
                continue
            # if k in og_seqs[i]:
            #     added_seqs += [k]
            #     final
            
            redundancy_check = [x for x in results if len(results[x][results[x]["target name"] == k]) == 1 and x != i]
            if not redundancy_check:
                if k in added_seqs or k in rejections1[j]:
                    continue
                hmm_cov_check = thresholds_existing[i][0][0] >= temp_cov[k][0] >= thresholds_existing[i][0][1]
                if not hmm_cov_check:
                    hmm_cov_check = thresholds_existing[i][0][0] < temp_cov[k][0]
                
                tgt_cov_check = thresholds_existing[i][1][0] >= temp_cov[k][1] >= thresholds_existing[i][1][1]
                if not tgt_cov_check:
                    tgt_cov_check = thresholds_existing[i][1][0] < temp_cov[k][1]
                
                if hmm_cov_check and tgt_cov_check:
                    # # if temp_cov[k][2] > repeat_threshold:
                    # #     # rejections[j] += [k]
                    # #     rejections.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                    #     rejections1[j] += [k]
                    #     ind1 += 1 
                    #     continue                    
                    # else:
                    added_seqs += [k]
                    final_seqs_for_tree[check] += [k]
                    continue     
                else:
                    rejections.loc[ind1] = [k, f"Either the HMM coverage ({temp_cov[k][0]}) or the target coverage ({temp_cov[k][1]}) fell below the thresholds established by the reciprocal best hit sequences. (hmm coverage threshold range - {thresholds_existing[i][0][0]}, {thresholds_existing[i][0][1]}, tgt coverage threshold range - {thresholds_existing[i][1][0]}, {thresholds_existing[i][1][1]})"]
                    rejections1[j] += [k]
                    ind1 += 1 
                    continue
            
            
            paralog_check = [x for x in redundancy_check if x not in assigned_ogs_updated[check]]
            if not paralog_check:
                if k in added_seqs or k in rejections1[j]:
                    continue
                # input()
                hmm_cov_check = thresholds_existing[i][0][0] >= temp_cov[k][0] >= thresholds_existing[i][0][1]
                if not hmm_cov_check:
                    hmm_cov_check = thresholds_existing[i][0][0] < temp_cov[k][0]
                
                tgt_cov_check = thresholds_existing[i][1][0] >= temp_cov[k][1] >= thresholds_existing[i][1][1]
                if not tgt_cov_check:
                    tgt_cov_check = thresholds_existing[i][1][0] < temp_cov[k][1]
                
                if hmm_cov_check and tgt_cov_check:
                    # if temp_cov[k][2] > repeat_threshold:
                    #     rejections.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                    #     rejections1[j] += [k]
                    #     ind1 += 1 
                    #     continue
                    # else:
                    added_seqs += [k]
                    final_seqs_for_tree[check] += [k]
                    continue 
                else:
                    rejections.loc[ind1] = [k, f"Either the HMM coverage ({temp_cov[k][0]}) or the target coverage ({temp_cov[k][1]}) fell below the thresholds established by the reciprocal best hit sequences (hmm coverage threshold range - {thresholds_existing[i][0][0]}, {thresholds_existing[i][0][1]}, tgt coverage threshold range - {thresholds_existing[i][1][0]}, {thresholds_existing[i][1][1]})"]
                    rejections1[j] += [k]
                    ind1 += 1 
                    continue
            
            temp_og = results[i][results[i]["target name"] == k]
            dist_og_cov_score = dist((0,0,0),temp_cov[k][:-2]+list(temp_og["best 1 domain score"]))
            # input()
            # dist_og_scores = dist((0,0), list(temp_og["full sequence score"]) + list(temp_og["best 1 domain score"]))
            
            for l in paralog_check:  
                # input()                
                temp_new = results[l][results[l]["target name"] == k]
                temp_cov_new = coverage_calc(list(temp_new["target name"]), results_dom[l], hmm_len)
                dist_new_cov_score = dist((0,0,0),temp_cov_new[k][:-2]+list(temp_new["best 1 domain score"]) )
                # dist_new_scores = dist((0,0), list(temp_new["full sequence score"]) + list(temp_new["best 1 domain score"]))
                
                # if dist_og_cov >= dist_new_cov:
                if dist_og_cov_score >= dist_new_cov_score:
                    if k in added_seqs or k in rejections1[j]:
                        continue                    
                #     if temp_cov[k][2] > repeat_threshold:
                #         rejections.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                #         ind1 += 1
                #         rejections1[j] += [k]
                #         continue      
                #     # if list(temp_og["best 1 domain score"][0]) >= list(temp_new["best 1 domain score"][0]):
                #     final_seqs_for_tree[check] += [k]
                #     added_seqs.append(k)
                # else:
                    if k in added_seqs or k in rejections1[j]:
                        continue
                    # if temp_cov_new[k][2] > repeat_threshold:
                    #     rejections.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                    #     rejections1[j] += [k]
                    #     ind1 += 1
                    #     continue   
                    check_new = [x for x in assigned_ogs_updated.keys() if l in assigned_ogs_updated[x]][0]                    
                    added_seqs.append(k)
                    final_seqs_for_tree[check_new] += [k]
        
    coverage_stats_new[i] = covs
    threshold_used_new[i] = threshold
    # rejected_seqs[i] = rejections
    # reject_seqs[i] = rejections1
del covs, threshold, rejections, i, j, k, l, temp, temp1, temp_cov, temp_cov_new, temp_ind, temp_new, temp_og, hmm_len, hmm_cov_check, tgt_cov_check, to_remove, check, check_new, dist_og_cov_score, dist_new_cov_score, jnb, coverage, rejections1, ind1, paralog_check
#%%
os.chdir(os.path.dirname(os.getcwd()))
os.chdir("post_hmmsearch_orthogroups_final")
seqdb = [x for x in SeqIO.parse(path2db, 'fasta')]
seqdb_ids = [x.id for x in seqdb]

temp = os.getcwd()
seqs = {}
for i in assigned_ogs_updated:
    os.chdir(i)
    temp1 = [x for x in os.listdir() if x.endswith("updated.fasta")]
    temp1 = [x for x in SeqIO.parse(temp1[0], 'fasta')]
    seqs[i] = temp1
    os.chdir(temp)
del temp, i, temp1
#%%
# seqs_to_write = {}
temp_seqs = [x for x in og_seqs.values()]
temp_seqs = [x for y in temp_seqs for x in y]
new_seqs = {}
for i in seqs:
    temp = [x.id for x in seqs[i]]
    temp1 = [x for x in final_seqs_for_tree[i] if x not in temp]
    temp1 = [x for x in temp1 if x not in temp_seqs]
    new_seqs[i] = temp1
    # seqs_to_write[i] = [x.id for x in seqs[i]] + temp1
del temp, temp1, i
#%%

# os.chdir(os.path.dirname(os.getcwd()))
os.mkdir("post_hmmsearch_orthogroups_with_new_species")
os.chdir("post_hmmsearch_orthogroups_with_new_species")
rejects = {}
temp = os.getcwd()
for i in assigned_ogs_updated:
    os.mkdir(i)
    # print(f"Adding sequences to {i}\n")
    os.chdir(i)    
    rejects[i] = [x for x in new_seqs[i] if x not in seqdb_ids]
    print(f"Writing sequences for {i}")
    with open(f"{i}.fasta", 'a') as f:
        [SeqIO.write(seqdb[seqdb_ids.index(x)], f, "fasta") for x in new_seqs[i]]
        [SeqIO.write(x, f, "fasta") for x in seqs[i]]
        f.close()
    os.chdir(temp)
del temp, i, f
#%%
import shutil
temp = os.getcwd()
for i in assigned_ogs_updated:
    os.chdir(i)
    os.mkdir("interproscan")
    shutil.copy(os.path.join(os.getcwd(), f"{i}.fasta"), os.path.join(os.path.join(os.getcwd(), "interproscan"), f"{i}.fasta"))
    os.chdir(temp)
del temp, i
#%%
temp = os.getcwd()
for i in assigned_ogs_updated:
    os.chdir(i)
    os.chdir("interproscan")
    [shutil.rmtree(x) for x in os.listdir() if x == "temp"]
    [os.remove(x) for x in os.listdir() if not x.endswith("fasta")]
    seqtemp = [x for x in SeqIO.parse([x for x in os.listdir() if x.endswith(".fasta")][0], 'fasta')]
    seqtemp_txt = ""
    seqtemp_txt = [seqtemp_txt.join(f">{x.id}\n{str(x.seq).replace('*', '')}\n") for x in seqtemp]
    with open([x for x in os.listdir() if x.endswith(".fasta")][0], 'w') as f:
        [f.write(x) for x in seqtemp_txt]
        f.close()
    os.chdir(temp)
del i, seqtemp, seqtemp_txt, f
