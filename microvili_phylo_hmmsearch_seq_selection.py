# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 22:37:15 2024

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


def coverage_calc(obj0:str, obj1:list, obj2:pd.DataFrame(), obj3:pd.DataFrame()):
    cov_df = {}
    
    for j in obj1:
        temp_df = obj2[f"{obj0}_0"][obj2[f"{obj0}_0"]["target name"] == j]
        hmm_cov, repeats = estimate_len_covered(sorted(zip(list(temp_df["hmm coord from"]), list(temp_df["hmm coord to"]))))
        hmm_cov = round(hmm_cov/obj3[f"{obj0}_0"], 2)
        
        
        if len(temp_df) == 0:
            cov_df[j] = [] 
            continue
        
        target_cov = round(estimate_len_covered(sorted(zip(list(temp_df["ali coord from"]), list(temp_df["ali coord to"]))))[0]/stat.mean(temp_df["target length"]),2)
        
        cov_df[j] = [hmm_cov, target_cov, repeats, stat.mean(temp_df["target length"])]
    return cov_df
# =============================================================================
# Orthogroups are input based on the existing literature
# =============================================================================

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
#%% =============================================================================
# This block of code serves to import all the speciws-wise HMM search results 
# which have been processed into a single CSV file. The results are imported 
# into a dictionary
# =============================================================================
base_list = os.listdir()
base_path = os.getcwd()
path2db = r"/g/dey/Comparativegenomics/Mylan/Opis_Choano_92/opis_choanos.fasta"

if os.name == "nt":
    path2db = serverpath2winpath(path2db)
    listofspecies, speciesdict = getlistofspeciesindb(path2db)
else:
    listofspecies, speciesdict = getlistofspeciesindb(path2db)


#
results = {}
results_dom = {}
results_hmm_len = {}
original_seqs = {}
added_sp = []
for i in base_list:
    os.chdir(i)
    original_seqs[i] = [x.id for x in SeqIO.parse([x for x in os.listdir() if x.endswith("updated.fasta")][0],'fasta')]
    hmm_len = stat.mean([len(x.seq) for x in SeqIO.parse([x for x in os.listdir() if x.endswith("trim.phy")][0], 'fasta')])
    
    try:
        os.chdir("HMM_searches")
    except FileNotFoundError:
        print(f"HMM searches seem to not have been completed for {i}. Please check")
        os.chdir(base_path)
        continue
    # results.append(i)
    check = max([int(x[-1]) for x in os.listdir()])
    os.chdir(f"hmm_searches_{check}")
    results_hmm_len[f"{i}_{check}"] = hmm_len
    temp = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch.csv")][0])
    temp1 = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch_doms.csv")][0])
    del temp["Unnamed: 0"]
    del temp1["Unnamed: 0"]
    temp1 = temp1.astype(object)
    temp.fillna("None", inplace=True)
    temp1.fillna("None", inplace=True)
    results[f"{i}_{check}"] = temp
    results_dom[f"{i}_{check}"] = temp1
    [added_sp.append(y) for y in list(dict.fromkeys([x[:x.find("_")] for x in results[f"{i}_{check}"]["target name"] if x != "None"]))]
    os.chdir(base_path)
del i, temp, hmm_len, temp1

added_sp = list(dict.fromkeys(added_sp))
assigned_ogs_updated = {x: [f"{z}_{check}" for z in y] for x, y in assigned_ogs_updated.items()}
os.chdir(os.path.dirname(base_path))
if [x for x in os.listdir() if x.startswith("new")]:
    path2db = r"/g/dey/Comparativegenomics/Mylan/Species_database/new_species_added_db_18032024.fsa"
    if os.name == "nt":
        path2db = serverpath2winpath(path2db)
        listofspecies1, speciesdict1 = getlistofspeciesindb(path2db)
    else:
        listofspecies1, speciesdict1 = getlistofspeciesindb(path2db)
    os.chdir([x for x in os.listdir() if x.startswith("new")][0])
    base_list = os.listdir()
    base_path = os.getcwd()
    
    for i in base_list:
        os.chdir(i)        
        temp = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch.csv")][0])
        temp1 = pd.read_csv([x for x in os.listdir() if x.endswith("hmmsearch_doms.csv")][0])
        del temp["Unnamed: 0"]
        del temp1["Unnamed: 0"]
        temp1 = temp1.astype(object)
        temp.fillna("None", inplace=True)
        temp1.fillna("None", inplace=True)
        temp = temp[temp["target name"].str.contains("|".join([x for x in speciesdict1  if x not in added_sp]) + "|None")]
        temp1 = temp1[temp1["target name"].str.contains("|".join([x for x in speciesdict1 if x not in added_sp]) + "|None")]
        results[f"{i}_{check}"] = pd.concat([results[f"{i}_{check}"], temp])
        results_dom[f"{i}_{check}"] = pd.concat([results_dom[f"{i}_{check}"], temp])
        os.chdir(base_path)
    del i, temp, check, temp1

listofspecies = list(dict.fromkeys(listofspecies + listofspecies1))
speciesdict.update(speciesdict1)
del listofspecies1, added_sp, speciesdict1
#%%
added_seqs = []
rejected_seqs = {}
og_seqs_missing_in_hmmsearch = {}
thresholds_used = {}
coverage_stats = {}
final_seqs_for_tree = dict([(x, []) for x in assigned_ogs_updated])
jnb = jpy.JenksNaturalBreaks(2)

reject_seqs = {}
print("Assigning sequenced identified in previous iteration with BBH approach to their orthogroups\n")
for i in assigned_ogs_updated:
    to_add = []
    missing = {}
    # threshols = {}
    for j in assigned_ogs_updated[i]:
        coverage = coverage_calc(j[:-2], original_seqs[j[:-2]], results_dom, results_hmm_len)
        coverage_stats[j] = coverage
        to_remove = [x for x in coverage if not coverage[x]]
        
        if to_remove:
            missing[j] = to_remove
        coverage = {x: y for x,y in coverage.items() if y}    
        min_hmm_cov = [x[0] for x in coverage.values()]
        min_target_cov = [x[1] for x in coverage.values()]
        length = [x[3] for x in coverage.values()]
            # hmm_cov_stdev = round(stat.stdev(min_hmm_cov), 3)
            # target_cov_stdev = round(stat.stdev(min_target_cov), 3)
        thresholds_used[j] = [(max(min_hmm_cov), min(min_hmm_cov)), (max(min_target_cov), min(min_target_cov)), (max(length), min(length))]
        [to_add.append(x) for x in original_seqs[j[:-2]]]
    final_seqs_for_tree[i] = to_add
    added_seqs += to_add
    # thresholds_used[i] = threshols
    og_seqs_missing_in_hmmsearch[i] = missing
del i, min_hmm_cov, min_target_cov, length, coverage, to_remove, to_add, j, missing

lens_comparison_og = [(len(x), x) for x in final_seqs_for_tree.values()]
#%%
print("Identifying new hits from hmmsearch and assigning them to their orthogroups\n")
# repeats_threshold = {}
threshold_used_new = {}
coverage_stats_new = {} 
for i in results:
    print(f"\n{i}")
    threshold = {}
    covs = {}
    check = [x for x in assigned_ogs_updated.keys() if i in assigned_ogs_updated[x]][0]
    rejections = {x: [] for x in speciesdict}
    rejections1 = pd.DataFrame(columns=["Sequence ID", "Reason"])
    ind1 = len(rejections)
    # input()
    # repeat_threshold = int(stat.mean([coverage_stats[i][x][2] for x in coverage_stats[i] if coverage_stats[i][x]]))
    # if repeat_threshold < 2:
    #     repeats_threshold[i] = (repeat_threshold, 3)
    #     repeat_threshold = 3
    # else:
    #     repeats_threshold[i] = (repeat_threshold, repeat_threshold)
   
    for j in speciesdict:        
        temp = results[i][results[i]["target name"].str.startswith(j)]
        temp1 = list(map(sum, zip(list(temp["full sequence score"]), list(temp["best 1 domain score"]))))
        
        if len(list(dict.fromkeys(temp1))) == 1:
            temp = list(temp["target name"])
            if temp[0] in added_seqs or temp[0] in rejections[j]:
                continue
            added_seqs += [temp[0]]
            final_seqs_for_tree[check] += [temp[0]]
            continue
            
        if len(list(dict.fromkeys(temp1))) > 2:
            jnb.fit(temp1)
            threshold[f"bitscore_threshold_{j}"] = [max(jnb.groups_[0])]
            temp_ind = np.add(temp["full sequence score"], temp["best 1 domain score"]) > max(jnb.groups_[0])
            temp = [temp["target name"][x] for x in temp_ind.index if temp_ind[x]]
            
        else:
            temp = list(temp["target name"])
        # temp = [x for x in temp if x not in added_seqs]
        if not temp:
            rejections[j] = []
            covs[j] = []
            rejections1.loc[ind1] = [[], f"No sequence found in {j}"]
            ind1 += 1
            # threshold[j]
            continue
        temp_cov = coverage_calc(i[:-2], temp, results_dom, results_hmm_len)
        to_remove = [x for x in temp_cov if not temp_cov[x] and x not in added_seqs]
        #CHECK THIS
        rejections[j] = to_remove
        rejections1.loc[ind1] = [to_remove, f"Sequences were not found in the hmmsearch for {j}"]
        ind1 += 1
        coverage = {x: y for x,y in temp_cov.items() if y} 
        
        covs[j] = coverage
        temp = list(coverage.keys())
        for k in temp:
            if k in added_seqs or k in rejections[j]:
                continue
            
            
            redundancy_check = [x for x in results if len(results[x][results[x]["target name"] == k]) == 1 and x != i]
            if not redundancy_check:
                if k in added_seqs or k in rejections[j]:
                    continue
                hmm_cov_check = thresholds_used[i][0][0] >= temp_cov[k][0] >= thresholds_used[i][0][1]
                if not hmm_cov_check:
                    hmm_cov_check = thresholds_used[i][0][0] < temp_cov[k][0]
                
                tgt_cov_check = thresholds_used[i][1][0] >= temp_cov[k][1] >= thresholds_used[i][1][1]
                if not tgt_cov_check:
                    tgt_cov_check = thresholds_used[i][1][0] < temp_cov[k][1]
                
                if hmm_cov_check and tgt_cov_check:
                    # if temp_cov[k][2] > repeat_threshold:
                    #     rejections[j] += [k]
                    #     rejections1.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                    #     ind1 += 1
                    #     continue                    
                    # else:
                        added_seqs += [k]
                        final_seqs_for_tree[check] += [k]
                        continue     
                else:
                    rejections[j] += [k]
                    rejections1.loc[ind1] = [k, f"Either the HMM coverage ({temp_cov[k][0]}) or the target coverage ({temp_cov[k][1]}) fell below the thresholds established by the reciprocal best hit sequences. (hmm coverage threshold range - {thresholds_used[i][0][0]}, {thresholds_used[i][0][1]}, tgt coverage threshold range - {thresholds_used[i][1][0]}, {thresholds_used[i][1][1]})"]
                    ind1 += 1
                    continue
            
            
            paralog_check = [x for x in redundancy_check if x not in assigned_ogs_updated[check]]
            if not paralog_check:
                if k in added_seqs or k in rejections[j]:
                    continue
                # input()
                hmm_cov_check = thresholds_used[i][0][0] >= temp_cov[k][0] >= thresholds_used[i][0][1]
                if not hmm_cov_check:
                    hmm_cov_check = thresholds_used[i][0][0] < temp_cov[k][0]
                
                tgt_cov_check = thresholds_used[i][1][0] >= temp_cov[k][1] >= thresholds_used[i][1][1]
                if not tgt_cov_check:
                    tgt_cov_check = thresholds_used[i][1][0] < temp_cov[k][1]
                
                if hmm_cov_check and tgt_cov_check:
                    # if temp_cov[k][2] > repeat_threshold:
                    #     rejections[j] += [k]
                    #     rejections1.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                    #     ind1 += 1
                    #     continue
                    # else:
                        added_seqs += [k]
                        final_seqs_for_tree[check] += [k]
                        continue 
                else:
                    rejections[j] += [k]
                    rejections1.loc[ind1] = [k, f"Either the HMM coverage ({temp_cov[k][0]}) or the target coverage ({temp_cov[k][1]}) fell below the thresholds established by the reciprocal best hit sequences (hmm coverage threshold range - {thresholds_used[i][0][0]}, {thresholds_used[i][0][1]}, tgt coverage threshold range - {thresholds_used[i][1][0]}, {thresholds_used[i][1][1]})"]
                    ind1 += 1
                    continue
            
            temp_og = results[i][results[i]["target name"] == k]
            dist_og_cov_score = dist((0,0,0),temp_cov[k][:-2]+list(temp_og["best 1 domain score"]))
            # input()
            # dist_og_scores = dist((0,0), list(temp_og["full sequence score"]) + list(temp_og["best 1 domain score"]))
            
            for l in paralog_check:  
                # input()                
                temp_new = results[l][results[l]["target name"] == k]
                temp_cov_new = coverage_calc(l[:-2], list(temp_new["target name"]), results_dom, results_hmm_len)
                dist_new_cov_score = dist((0,0,0),temp_cov_new[k][:-2]+list(temp_new["best 1 domain score"]) )
                # dist_new_scores = dist((0,0), list(temp_new["full sequence score"]) + list(temp_new["best 1 domain score"]))
                
                # if dist_og_cov >= dist_new_cov:
                if dist_og_cov_score >= dist_new_cov_score:
                    if k in added_seqs or k in rejections[j]:
                        continue                    
                    # if temp_cov[k][2] > repeat_threshold:
                    #     rejections[j] += [k]
                    #     rejections1.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                    #     ind1 += 1
                    #     continue      
                    # # if list(temp_og["best 1 domain score"][0]) >= list(temp_new["best 1 domain score"][0]):
                    final_seqs_for_tree[check] += [k]
                    added_seqs.append(k)
                else:
                    if k in added_seqs or k in rejections[j]:
                        continue
                    # if temp_cov_new[k][2] > repeat_threshold:
                    #     rejections[j] += [k]
                    #     rejections1.loc[ind1] = [k, f"{k} likely contains a lot of repeats ({temp_cov[k][2]}) as determined by the code. The repeat threshold for this HMM was {repeat_threshold}. This was determined by checking how many times the same segment of the HMM mapped to multiple regions in the protein"]
                    #     ind1 += 1
                    #     continue   
                    check_new = [x for x in assigned_ogs_updated.keys() if l in assigned_ogs_updated[x]][0]                    
                    added_seqs.append(k)
                    final_seqs_for_tree[check_new] += [k]
        
    coverage_stats_new[i] = covs
    threshold_used_new[i] = threshold
    rejected_seqs[i] = rejections
    reject_seqs[i] = rejections1
del i, j, temp, temp1, threshold, jnb, check, temp_new, temp_og, redundancy_check, paralog_check, k, l, covs, to_remove, rejections, temp_cov, temp_ind, dist_og_cov_score, dist_new_cov_score, temp_cov_new, check_new, coverage, hmm_cov_check, tgt_cov_check, ind1

lens_comparison_new = [(len(x), x) for x in final_seqs_for_tree.values()]

#%%
os.chdir(os.path.dirname(base_path))
os.mkdir("post_hmmsearch_orthogroups_final")
os.chdir("post_hmmsearch_orthogroups_final")
#
path2db = r"/g/dey/Comparativegenomics/Mylan/Opis_Choano_92/opis_choanos.fasta"
if os.name == 'nt':
    path2db = serverpath2winpath(path2db)
    seqdb = [x for x in SeqIO.parse(path2db, 'fasta')]
else:
    seqdb = [x for x in SeqIO.parse(path2db, 'fasta')]
# seqdb_ids = [x.id for x in seqdb]
path2db = r"/g/dey/Comparativegenomics/Mylan/Species_database/new_species_added_db_18032024.fsa"
if os.name == 'nt':
    path2db = serverpath2winpath(path2db)
    seqdb += [x for x in SeqIO.parse(path2db, 'fasta')]
else:
    seqdb += [x for x in SeqIO.parse(path2db, 'fasta')]
seqdb_ids = [x.id for x in seqdb]
#
temp = os.getcwd()
for i in final_seqs_for_tree:
    os.mkdir(i)
    os.chdir(i)
    print(f"Writing sequences for {i}")
    with open(f"{i}.fasta", 'w') as f:
        [SeqIO.write(seqdb[seqdb_ids.index(x)], f, "fasta") for x in final_seqs_for_tree[i]]
        f.close()
    os.chdir(temp)
del temp, i
#%%
os.chdir(os.path.dirname(base_path))
os.chdir([x for x in os.listdir() if x.startswith("log")][0])
with open(f"hmmsearch_log_{datetime.today().strftime('%Y_%m_%d')}.txt", 'w') as f:
    # f.write("For a given list of proteins, HMMs were built with default hmmbuild parameters using the best, bidirectional hits (BBHs) whose orthology was verified through maximum likelihood gene tree construction. These HMMs were then used to search for homologs in the proteome of each species separately.\n\nFollowing HMM searches, paralogous relationships between the originally given list of proteins were inferred using a clustering approach that incorporated three measures:\n\n1. extent of overlap of search results between different HMM searches (how many hits from the search results of protein A can also be found in the search results of protein B),\n2. Average full sequence bitscore of the overlapping search results, and\n3. Average best one domain bitscore of the overlapping search results\n\nThe resulting orthogroups are:\n\n")
    # f.write("----Thresholds used----")
    # f.write(f"\nMinimum overlap percentage between results\t{overlap_minimum_dist}\n")
    # f.write(f"\nMaximum coefficient of variation in domain bitscores\t{overlap_max_score_coev}\n\n")
    f.write("----Orthologous groups----")
    f.write("\nOrthogroup name\tMember proteins\n")
    for i in assigned_ogs_updated:
        f.write(f"{i}\t")
        [f.write(f"{x[:-2]}\t") for x in assigned_ogs_updated[i]]
        f.write("\n")
    f.write("--------\n")    
    f.write(f"\nThe following species (total = {len(speciesdict)}) were used in this analysis:\n")
    f.write("\n----Species used----\n")
    f.write("Species name\tAbbreviation\n")
    for i in speciesdict:
        f.write(f"{speciesdict[i]}\t{i}\n")
    f.write("--------\n\n")
    f.write("The following proteins were recovered and assigned to each orthogroup:\n\n")
    f.write("----Orthogroup assignments----\n\n")
    for i in final_seqs_for_tree:
        f.write(f"Orthogroup {i} - ")
        [f.write(f"{x} ") for x in assigned_ogs_updated[i]]
        f.write(f"\t(Number of proteins = {len(final_seqs_for_tree[i])})\n\n")
        temp = final_seqs_for_tree[i]
        temp.sort()
        [f.write(f"{x}\n") for x in temp]
        f.write("\n")
    f.write("\n----Rejected proteins----\n")
    # f.write("Since the HMMs were custom-built and sequence-specific rather than domain-specific, a t")
    f.write("Protein ID\tReason")
    for i in reject_seqs:
        f.write(f"\n\n{i}\n")
        for j in reject_seqs[i].index:
            if reject_seqs[i]["Sequence ID"][j]:
                f.write(f"{reject_seqs[i]['Sequence ID'][j]}\t{reject_seqs[i]['Reason'][j]}\n")
    f.close()
