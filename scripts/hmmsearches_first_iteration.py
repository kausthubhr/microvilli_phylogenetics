# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 15:35:59 2023

@author: Kausthubh R
"""

import os
import subprocess as sp
import pandas as pd
import shutil
import copy
from concurrent.futures import ThreadPoolExecutor as tpe
#%
def exec_(cmd):
    '''
    single command run
    command must be prepared as subprocess.call
    '''
    sp.getoutput(cmd)

        
def pull_run(obj1, obj2):
    # print(len(obj2))
    if len(obj2) > obj1:
        with tpe(max_workers=obj1) as exe:
            futures = exe.map(exec_, obj2)
            # wait(futures)
    else:
        with tpe(max_workers=len(obj2)) as exe:
            futures = exe.map(exec_, obj2)
        
    return futures

def setpath():  
    choice = 'y'
    print(os.getcwd())
    choice = input("Do you wish to change your current working directory? [y/n]: \n")    
    if (choice == 'n'):
        print("The working directory will not be changed")        
        temp1 = 'n'
        temp1 = input("Do you want to list all the directories? [y/n]: \n")
        if(temp1 == 'y'):
            dir_list = os.listdir()
            for i in range(len(dir_list)):
                print("\n{0} - {1}\n".format(int(i+1), dir_list[i]))
            dir_select = int(input("Please enter the number correspoding to the directory you want to select: \n"))
            dir_select = dir_select - 1
            return dir_list[dir_select]
        else:
            return os.getcwd()
        
    if (choice == 'y'):
        path1 = input("Enter the path to your query file/target database/program: \n" )
        os.chdir(path1)
        
        temp1 = 'n'
        temp1 = input("Do you want to list all the directories? [y/n]: \n")
        if(temp1 == 'y'):
            dir_list = os.listdir()
            for i in range(len(dir_list)):
                print("\n{0} - {1}\n".format(int(i+1), dir_list[i]))
            dir_select = int(input("Please enter the number correspoding to the directory you want to select: \n"))
            dir_select = dir_select - 1
            return dir_list[dir_select]
        else:
            return os.getcwd()

#obj1 is the name of the tblout file
def read_phmmer_tblout(obj1):
    temp = open(obj1, 'r')
    temp_con = temp.read().split("\n")[0:-10]
    temp1 = []
    for i in range(len(temp_con[1:])):
        temp2 = copy.deepcopy(temp_con[i+1])
        temp2 = temp2.split(" ")
        temp2 = list(filter(None, temp2))
        temp1.append(temp2)
    del i, temp2
    for i in temp1:
        try:
            for j in range(0,18):
                try:
                    i[j] = float(i[j])
                except ValueError:
                    continue
        except IndexError:
            continue
    del i, j, temp1[1], temp1[0][0]
    #%
    temp1[0][3] = temp1[0][3] + " " + temp1[0][4]
    temp1[0][0] = temp1[0][0] + " " + temp1[0][1]
    temp1[0][20] = temp1[0][20] + " " + temp1[0][21] + " " + temp1[0][22]
    temp1[0][6] = "full sequence " + temp1[0][6]
    temp1[0][7] = "full sequence " + temp1[0][7]
    temp1[0][8] = "full sequence " + temp1[0][8]
    temp1[0][9] = "best 1 domain " + temp1[0][9]
    temp1[0][10] = "best 1 domain " + temp1[0][10]
    temp1[0][11] = "best 1 domain " + temp1[0][11]
    temp1[0][12] = "domain number estimation " + temp1[0][12]
    temp1[0][13] = "domain number estimation " + temp1[0][13]
    temp1[0][14] = "domain number estimation " + temp1[0][14]
    temp1[0][15] = "domain number estimation " + temp1[0][15]
    temp1[0][16] = "domain number estimation " + temp1[0][16]
    temp1[0][17] = "domain number estimation " + temp1[0][17]
    temp1[0][18] = "domain number estimation " + temp1[0][18]
    temp1[0][19] = "domain number estimation " + temp1[0][19]
    
    del  temp1[0][22], temp1[0][21], temp1[0][4], temp1[0][1], temp1[int(len(temp1)-1)]
    
    for i in temp1[1:]:
        for j in range(19, len(i)):
            i[18] = i[18] + " " + i[j]
            i[j] = []
        
    del i, j
    for i in range(len(temp1)):
        temp1[i] = [ele for ele in temp1[i] if ele != []]
    del i
    
    data = pd.DataFrame(temp1[1:len(temp1)], columns=temp1[0])
    return data

#obj1 is the name of the domtblout file
def read_phmmer_domtblout(obj1):
    temp = open(obj1, 'r')
    temp_con = temp.read().split("\n")[0:-10]
    temp1 = []
    for i in range(len(temp_con[1:])):
        temp2 = copy.deepcopy(temp_con[i+1])
        temp2 = temp2.split(" ")
        temp2 = list(filter(None, temp2))
        temp1.append(temp2)
    del i, temp2
    for i in temp1:
        try:
            for j in range(0,18):
                try:
                    i[j] = float(i[j])
                except ValueError:
                    continue
        except IndexError:
            continue
    del i, j, temp1[1], temp1[0][0]
    #%
    temp1[0][0] += " " + temp1[0][1]
    temp1[0][3] = "target length"  
    temp1[0][4] += " " + temp1[0][5]
    temp1[0][7] = "query length"
    temp1[0][8] = "full sequence " + temp1[0][8]
    temp1[0][9] = "full sequence " + temp1[0][9]
    temp1[0][10] = "full sequence " + temp1[0][10]
    temp1[0][11] = "this domain " + temp1[0][11]
    temp1[0][12] = "this domain " + temp1[0][12]
    temp1[0][13] = "this domain " + temp1[0][13]
    temp1[0][14] = "this domain " + temp1[0][14]
    temp1[0][15] = "this domain " + temp1[0][15]
    temp1[0][16] = "this domain " + temp1[0][16]
    temp1[0][17] = "hmm coord " + temp1[0][17]
    temp1[0][18] = "hmm coord " + temp1[0][18]
    temp1[0][19] = "ali coord " + temp1[0][19]
    temp1[0][20] = "ali coord " + temp1[0][20]
    temp1[0][21] = "envelope coord " + temp1[0][21]
    temp1[0][22] = "envelope coord " + temp1[0][22]
    temp1[0][23] = "mean posterior probability (acc)" # in maximum expected accuracy  (mea) alignment
    temp1[0][24] += " " + temp1[0][25] + " " + temp1[0][26]

    del temp1[0][26], temp1[0][25], temp1[0][5], temp1[0][1], temp1[int(len(temp1)-1)]

    for i in range(1, len(temp1)):
        for j in range(23, len(temp1[i])):
            temp1[i][22] += " " + temp1[i][j]
            temp1[i][j] = []
        temp1[i] = [x for x in temp1[i] if x != []]
    del i, j


    data = pd.DataFrame(temp1[1:len(temp1)], columns=temp1[0])
    return data

def remove_sequence_duplicates(obj1):
    temp = [x.seq for x in obj1]
    #get ids of duplicate sequences
    if len(temp) <= 100:
        res = [idx for idx, val in enumerate(temp) if val in temp[:idx]]
    if len(temp) > 100:
        oc_set = set()
        res = []
        for idx, val in enumerate(temp):
            if val not in oc_set:
                oc_set.add(val)        
            else:
                res.append(idx) 
    new = [obj1[x] for x in range(len(obj1)) if x not in res]
    removed = [obj1[x] for x in res]
    return new, removed

if os.name == 'nt':
    def winpath2linuxpath(obj1):
        temp = obj1
        temp = temp.replace("\\", "/")
        temp = "/mnt/" + temp[0].lower() + "/" + temp[3:]
        return temp
    
    def serverpath2winpath(obj1):
        temp = obj1
        temp = r"Z:" + temp[6:]
        temp = temp.replace("/", r"\\")
        return temp
    
    def winpath2serverpath(obj1):
        temp = obj1
        temp = temp.replace("Z:", "/g/dey")
        temp = temp.replace("\\", "/")
        temp = temp.replace("//", "/")
        return temp
    
#%%
print("Please ensure that this code is run with the folder containing your bbh run as a starting point\n\n")
os.chdir(r"verified_hits_0")
base_list = os.listdir()
base_path = os.getcwd()

print("Please navigate to the directory containing your proteome files (.fasta). The path to the directory is the only thing required for the code\n\n")
path2db = setpath()

if os.name == "nt":
    path2databases = [winpath2linuxpath(os.path.join(serverpath2winpath(path2db), x)) for x in os.listdir(serverpath2winpath(path2db)) if not x.startswith(".")]
else:
    path2databases = [os.path.join(path2db, x) for x in os.listdir(path2db)]

#%
for i in base_list:
    os.chdir(i)
    path2start = os.getcwd()
    hmmsearch_check = [x for x in os.listdir() if x == "HMM_searches"]
    if hmmsearch_check:
        print(f"HMM searches appear to have been completed for {i}. Proceeding to the next protein\n")
        os.chdir(base_path)
        continue
    os.mkdir("HMM_searches")
    os.chdir("HMM_searches")
    
    hmmsearch_check = [x for x in os.listdir() if x.startswith("hmm_searches")]
    if not hmmsearch_check:
        iteration = 0
    
    os.mkdir(f"hmm_searches_{iteration}")
    os.chdir(f"hmm_searches_{iteration}")
    path2hmm = os.getcwd()
    
    os.mkdir("hmmbuild_files")
    path2hmmbuild = os.path.join(os.getcwd(), "hmmbuild_files")
    
    os.chdir(path2start)
    [shutil.copy(x, path2hmmbuild) for x in os.listdir() if x.endswith("_updated.fasta")]
    [shutil.copy(x, path2hmmbuild) for x in os.listdir() if x.endswith("_trim.phy")]
    [shutil.copy(x, path2hmmbuild) for x in os.listdir() if x.endswith("_align.phy")]
    
    os.chdir(path2hmmbuild)
    align_name = [x for x in os.listdir() if x.endswith("_trim.phy")][0]
    untrim_align = [x for x in os.listdir() if x.endswith("_align.phy")][0]    
    print(f"Building HMMs for {i}\n")
    if os.name == 'nt':
        sp.getoutput("wsl hmmbuild -o {0}_log.txt -O {1}_updated_post_hmmbuild.phy {2}.hmm {3}".format(align_name, align_name, align_name[:-4], align_name))
    else:
        sp.getoutput("hmmbuild -o {0}_log.txt -O {1}_updated_post_hmmbuild.phy {2}.hmm {3}".format(align_name, align_name, align_name[:-4], align_name))
    
    check = [x for x in os.listdir() if x.endswith(".hmm")]
    
    if not check:
        print(f"The HMM has not been built in {i}. Please check. Proceeding to the next protein\n")
        os.chdir(base_path)
        continue
    
    print(f"Performing HMM searches for {i}\n")
    
    if os.name == "nt":
        hmm = winpath2linuxpath(os.path.join(path2hmmbuild, [x for x in os.listdir() if x.endswith(".hmm")][0]))
    else:
        hmm = os.path.join(path2hmmbuild, [x for x in os.listdir() if x.endswith(".hmm")][0])
    
    os.chdir(path2hmm)
    for j in path2databases:
        db_name = j[j.rfind("/")+1:][:-6]
    if os.name == 'nt':
        cmd_list = [("wsl hmmsearch --tblout {0}_iteration_{1}_hmmsearch_{2}_log.txt --domtblout {0}_iteration_{1}_hmmsearch_log_{2}_dom.txt -E 1e-10 {3} {4}".format(i, iteration, x[x.rfind("/")+1:], hmm, x)) for x in path2databases]
    else:
        cmd_list = [("hmmsearch --tblout {0}_iteration_{1}_hmmsearch_{2}_log.txt --domtblout {0}_iteration_{1}_hmmsearch_log_{2}_dom.txt -E 1e-10 {3} {4}".format(i, iteration, x[x.rfind("/")+1:], hmm, x)) for x in path2databases]
        
    pull_run(os.cpu_count(), cmd_list)
    print(f"Combining HMM search results for {i} into a single CSV file\n")
    list_hmm = [x for x in os.listdir() if x.endswith(".txt") and "dom" not in x] #fetching all individual results into one list of dataframes
    list_hmm_df = []
    for j in list_hmm:
        try:
            df = read_phmmer_tblout(j)
        except UnboundLocalError:
            df = j
        list_hmm_df.append(df)
    del j, df, list_hmm
    list_hmm_doms = [x for x in os.listdir() if x.endswith(".txt") and "dom" in x]
    list_hmm_doms_df = []
    for j in list_hmm_doms:
        try:
            df = read_phmmer_domtblout(j)
        except UnboundLocalError:
            df = j
        list_hmm_doms_df.append(df)
    del j, df, list_hmm_doms
    
    os.mkdir(f"{i}_individual_searches")
    [os.rename(x, os.path.join(f"{i}_individual_searches", x)) for x in os.listdir() if x.endswith(".txt")]
    
    for j in range(len(list_hmm_df)): #combining all individual results into one large dataframe
        if type(list_hmm_df[j]) is not str:
            combined_hmm_df = list_hmm_df[j].iloc[:0]
            combined_hmm_doms_df = list_hmm_doms_df[j].iloc[:0]
            break
    del j
    for j in range(len(list_hmm_df)):
        if type(list_hmm_df[j]) is str:
            species_name = list_hmm_df[j][:list_hmm_df[j].find(".")][list_hmm_df[j][:list_hmm_df[j].find(".")].rfind("hmmsearch")+10:].replace("_", " ")        
            append_list_temp = [None]*(len(combined_hmm_df.columns)-1)+[f"No {i} homolog found in [{species_name}]"]
            append_list_temp_1 = [None]*(len(combined_hmm_doms_df.columns)-1)+[f"No {i} homolog found in [{species_name}]"]
            combined_hmm_df.loc[-1] = append_list_temp
            combined_hmm_doms_df.loc[-1] = append_list_temp_1
            combined_hmm_df.reset_index(drop=True, inplace=True)
            combined_hmm_doms_df.reset_index(drop=True, inplace=True)
        if type(list_hmm_df[j]) is not str:
            combined_hmm_df = pd.concat([combined_hmm_df, list_hmm_df[j]])
            combined_hmm_doms_df = pd.concat([combined_hmm_doms_df, list_hmm_doms_df[j]])
            combined_hmm_df.reset_index(drop=True, inplace=True)
            combined_hmm_doms_df.reset_index(drop=True, inplace=True)
    

    combined_hmm_df.to_csv(f"{i}_combined_hmmsearch.csv")
    combined_hmm_doms_df.to_csv(f"{i}_combined_hmmsearch_doms.csv")
    print("\n\n\n")
    os.chdir(base_path)
#%%
