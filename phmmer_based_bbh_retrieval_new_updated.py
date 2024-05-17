# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:02:02 2024

@author: Kausthubh R
"""

import os
from Bio import SeqIO
import pandas as pd
import copy
import subprocess as sp
import csv
from tqdm import tqdm
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor as tpe
from datetime import datetime
import tracemalloc


#% function cell
def exec_(cmd):
    '''
    single command run
    command must be prepared as subprocess.call
    '''
    sp.getoutput(cmd)

        
def pull_run(obj1, obj2):
    print(f"Number of searches being performed: {len(obj2)}\n")
    if len(obj2) > obj1:
        with tpe(max_workers=obj1) as exe:
            futures = exe.map(exec_, obj2)
            # wait(futures)
    else:
        with tpe(max_workers=len(obj2)) as exe:
            futures = exe.map(exec_, obj2)
        
    return futures

def search_completenes_check(obj1):
    os.chdir(obj1)
    incomplete_searches = []
    temp_list_reverse_searches = os.listdir()
    for i in temp_list_reverse_searches:
        if i[0] != ".":
            os.chdir(i)
        for j in os.listdir():
            if os.stat(j).st_size == 0:
                if os.name =='nt':
                    incomplete_searches.append(winpath2linuxpath(os.path.join(os.getcwd(), j)))
                else:
                    incomplete_searches.append(os.path.join(os.getcwd(), j))
        os.chdir(obj1)
    return incomplete_searches


#setting path to required file/program working directory
#if you don't change the directory, this function will return the current working directory as a string

#Checking the species available in the target database.
#obj1 should be the path to the target database (string).
#function returns a list of species IDs
#Target database should be in .fasta file format. Each sequence header must end with the species name in square brackets []

def getlistofspeciesindb(obj1):
    listsp = []
    spdict = []
    for record in tqdm(SeqIO.parse(obj1, 'fasta'), unit=" sequences"):
        x = record.description
        y = record.id
        ind1 = x.rfind("]")
        ind2 = x.rfind("[")
        if(ind1<10):
            continue
        else:
            listsp.append(x[(ind2+1):(ind1)])
            spdict.append(y[:y.find("_")])
    print(len(listsp))
    print(len(spdict))
    spdict = OrderedDict.fromkeys(zip(spdict, listsp))
    listsp = list(dict.fromkeys(listsp))
    spdict = dict(OrderedDict(spdict))
    spdict = dict([(i[0], i[1]) for i in spdict])
    return listsp, spdict

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
     
    if len(temp1) == 0:
       return None
     
    for i in temp1[1:]:
        for j in range(19, len(i)):
            i[18] = i[18] + " " + i[j]
            i[j] = []
        
    # del i, j
    for i in range(len(temp1)):
        temp1[i] = [ele for ele in temp1[i] if ele != []]
    del i
     
    data = pd.DataFrame(temp1[1:len(temp1)], columns=temp1[0])
    return data

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
    
# obj1 is the path to the folder containing reverse searches
# obj2 is the path to the isoform dictionary
def retrievebbh(obj1, obj2):
    os.chdir(obj1) # function input
    temp_list = os.listdir()
    list_of_verified_hits = {}
    for i in temp_list:
        if i[0] != ".":
            print(f"Retriving the BBHs for {i}")
            # list_of_verified_hits.append(i)
            temp_verified = []
            with open(obj2) as f:  #function input
                x = csv.reader(f, delimiter='\t')
                for line in x:
                    for word in line:
                        if "|{0}".format(i) in word:
                            query_overall = line
                            break
                f.close()
            query_overall = list(filter(None, query_overall))
            os.chdir(i)
            temp_list1 = [x for x in os.listdir() if not x.startswith("dom_")]
            for j in temp_list1:
                temp_df = read_phmmer_tblout(j)
                if type(temp_df) is None:
                    os.chdir(obj1)
                    continue
                temp_id_list = list(temp_df.loc[:,'target name'])
                temp_bitscores = list(temp_df.loc[:,'full sequence score'])
                temp_domain_scores = list(temp_df.loc[:,'best 1 domain score'])
                # del temp_bitscores[0], temp_domain_scores[0]
                if not temp_bitscores:
                    continue
                # best_seq_score_id = temp_id_list[temp_bitscores.index(max(temp_bitscores))]
                best_dom_score_id = temp_id_list[temp_domain_scores.index(max(temp_domain_scores))]
                # if best_seq_score_id == best_dom_score_id:
                #     best_id = best_seq_score_id
                #     if best_id in query_overall:
                #         temp_verified.append(j[:j.find(".fasta")+6])
                # else:
                best_id = best_dom_score_id
                if best_id in query_overall:
                    temp_verified.append(j[:j.find("reverse")-1] + ".fasta")
            temp_verified = list(dict.fromkeys(temp_verified))
            temp_verified = list(filter(None, temp_verified))
            list_of_verified_hits[i] = temp_verified
            os.chdir(obj1)
            # del i, temp_list1, best_id, best_seq_score_id, best_dom_score_id, temp_list, temp_verified, temp_bitscores, temp_id_list, temp_domain_scores, f, word, x, temp_df, query_overall, j, line

    return list_of_verified_hits
        
#%%
print("\nCreate the directory where you wish to store the results of the run\n\n")
print("\n\nSTEP A: Decide where in your system you want to store the results of your searches. All the results of any one run are contained in a single folder with a name and location determined by the user.\n\nThe script will now tell you what your current location in the system is, and will ask you if you want to change it. Please follow the instructions.\n\nIt is good practice to avoid spaces/gaps in your filenames. Please consider using _ (underscore) instead.\n\n")
print(f"Current folder location = {os.getcwd()}\n\n")
choice = input("Do you want to change your location to a different folder: [y/n] \n")
if choice != 'y':
    pass
if choice =='y':
    temp_path = input("Please enter the path to your new location:\n")
    os.chdir(temp_path)

run_folder = input("\nGive a name for the folder where you will be storing the results of this run:\n\n")


os.mkdir(run_folder)
path2runfolder = os.path.join(os.getcwd(), run_folder)

os.chdir(path2runfolder)
os.mkdir(r"query_sequences")
path2queryseqs = os.path.join(path2runfolder, r"query_sequences")
if os.name == 'nt':
    linuxpath2queryseqs = winpath2linuxpath(path2queryseqs) + '/'

print("\n\n\nFor the subsequent steps, it might be more convenient if you, the user, can have all relevant files in one folder\n\n")

print("STEP B: Select a target database to search. \n\nThis should be one FASTA file containing a set of sequences that can be from one or more species. Usually, this is a collection of proteomes (one proteome per species) that has previously been prepared using the script uniprot_fasta_id_modifier.py. For further details, please consult the Google doc.\n\nIf your target database is in a different folder than the one you are currently at, please navigate to that folder\n\n")
print(f"Current folder location = {os.getcwd()}\n\n")
choice = input("Do you want to change your location to a different folder: [y/n] \n")
if choice != 'y':
    pass
if choice =='y':
    temp_path = input("Please enter the path to your new location:\n")
    os.chdir(temp_path)

print("\n\nTo make the target database selection easier, the script will now print all the available files in the folder and ask you to select one:\n\n")
dir_list = os.listdir()
for i in range(len(dir_list)):
    print("\n{0} - {1}\n".format(int(i+1), dir_list[i]))
dir_select = int(input("\nPlease enter the number corresponding to the target database: \n"))
dir_select = dir_select - 1
targetdbname = dir_list[dir_select]
path2targetdb = os.path.join(os.getcwd(), targetdbname)

if os.name == 'nt':
    linuxpath2targetdb = winpath2linuxpath(path2targetdb)
listofspecies, speciesdict = getlistofspeciesindb(path2targetdb)
del targetdbname, dir_list, dir_select, i 
#%
print("\n\nSTEP C: Select a starting database that will be the source of your query sequences. \n\nThis should be one FASTA file containing a set of sequences that belongs to a single species. Currently, this has to be a reference proteome from UniProt. For further details, please consult the Google doc.\n\nIf your starting database is in a different folder than the one you are currently at, please navigate to that folder\n\n")
print(f"Current folder location = {os.getcwd()}\n\n")
choice = input("Do you want to change your location to a different folder: [y/n] \n")
if choice != 'y':
    pass
if choice =='y':
    temp_path = input("Please enter the path to your new location:\n")
    os.chdir(temp_path)

print("\n\nTo make the source database selection easier, the script will now print all the available files in the folder and ask you to select one:\n\n")
dir_list = os.listdir()
for i in range(len(dir_list)):
    print("\n{0} - {1}\n".format(int(i+1), dir_list[i]))
dir_select = int(input("\nPlease enter the number corresponding to the source database: \n"))
dir_select = dir_select - 1
sourcedbname = dir_list[dir_select]
path2sourcedb = os.path.join(os.getcwd(), sourcedbname)

if os.name == 'nt':
    linuxpath2sourcedb = winpath2linuxpath(path2sourcedb)
del sourcedbname, dir_list, dir_select, i 

print("\n\nSTEP D: Select a CSV file containing a list of identifiers. \n\nThese identifiers should have been taken from the starting database. The file should have been formatted using uniprot_id_extraction.py. For further details, please consult the Google doc.\n\nIf your list is in a different folder than the one you are currently at, please navigate to that folder\n")
print(f"Current folder location = {os.getcwd()}\n\n")
choice = input("Do you want to change your location to a different folder: [y/n] \n")
if choice != 'y':
    pass
if choice =='y':
    temp_path = input("Please enter the path to your new location:\n")
    os.chdir(temp_path)

print("\n\nTo make the list selection easier, the script will now print all the available files in the folder and ask you to select one:\n\n")
dir_list = os.listdir()
for i in range(len(dir_list)):
    print("\n{0} - {1}\n".format(int(i+1), dir_list[i]))
dir_select = int(input("\nPlease enter the number corresponding to the list of queries: \n"))
dir_select = dir_select - 1
queryfilename = dir_list[dir_select]
querylist = os.path.join(os.getcwd(), queryfilename)

path2querylist = os.getcwd()
del queryfilename, dir_list, dir_select, i

print("\n\nSTEP E: Select a file containing an isoform dictionary. \n\nThis file should have been formatted using isoform_grouping.py. For further details, please consult the Google doc.\n\nIf your dictionary is in a different folder than the one you are currently at, please navigate to that folder\n")
print(f"Current folder location = {os.getcwd()}\n\n")
choice = input("Do you want to change your location to a different folder: [y/n] \n")
if choice != 'y':
    pass
if choice =='y':
    temp_path = input("Please enter the path to your new location:\n")
    os.chdir(temp_path)

print("\n\nTo make the dictionary selection easier, the script will now print all the available files in the folder and ask you to select one:\n\n")
dir_list = os.listdir()
for i in range(len(dir_list)):
    print("\n{0} - {1}\n".format(int(i+1), dir_list[i]))
dir_select = int(input("\nPlease enter the number corresponding to the isoform dictionary: \n"))
dir_select = dir_select - 1
isodictname = dir_list[dir_select]
path2isodict = os.path.join(os.getcwd(), isodictname)
del isodictname, dir_select, dir_list, i

 
#%%
os.chdir(r"query_sequences")
queries = open(querylist)
queries = queries.read().split(",")
queries = list(filter(None, queries))
queries = list(dict.fromkeys(queries))

print("\nExtracting query sequences\n")
p = []
query_seqs = []
query_seq_ids = []
for record in SeqIO.parse(open(path2sourcedb), 'fasta'):
    p.append(record)
for i in queries:
    for record1 in p:
        if i == record1.id:
            query_seqs.append(record1)
            query_seq_ids.append(record1.id)
            break
del i, record, record1, p

#checking if all ids in querylist have sequences
missing_seq_ids = []
for i in queries:
    if i not in query_seq_ids:
        missing_seq_ids.append(i)
        print(i)
del i, query_seq_ids

if not missing_seq_ids:
    print("\nSequences for all submitted query IDs have been extracted\n")
    del missing_seq_ids

print("\nWriting query sequences into the created directory\n")
# os.chdir(path2queryseqs)
for i in query_seqs:
    ind1 = i.id.rfind("|")
    ind1 = ind1 + 1
    SeqIO.write(i, "{0}.fasta".format(i.id[ind1:]), 'fasta')
del i, queries, query_seqs, ind1

queryfilelist = os.listdir()

print("\nInitiating forward searches using phmmer from the HMMER suite of tools\n")
os.chdir(path2runfolder)
iteration = 0
forward_searches = "forward_searches_{0}_".format(iteration) + run_folder
os.mkdir(forward_searches)
path2forwardsearches = os.path.join(path2runfolder, forward_searches)
os.chdir(path2forwardsearches)

for i in range(len(queryfilelist)):
    try:
        open('{0}.txt'.format(str(queryfilelist[i][:-6] +"_forward_searches")), 'a').close()
    except OSError:
        print('Failed creating the file. Possible reasons include lack of memory.')
    # queryfilelist[i] = os.path.join()
del i
forward_searches_files_list = [x for x in os.listdir() if not x.startswith("dom")] #check if needed after reverse searches cell
#%
if os.name == 'nt':
    cmd_list = ['wsl phmmer --tblout {0} --domtblout dom_{0} {1} {2}'.format(forward_searches_files_list[i], (os.path.join(path2queryseqs, queryfilelist[i])), winpath2linuxpath(path2targetdb)) for i in range(len(forward_searches_files_list))]
else:
    cmd_list = ['phmmer --tblout {0} --domtblout dom_{0} {1} {2}'.format(forward_searches_files_list[i], (os.path.join(path2queryseqs, queryfilelist[i])), (path2targetdb)) for i in range(len(forward_searches_files_list))]

with open(os.path.join(path2runfolder, f"{run_folder}_forward_searches_memory_usage.txt"), 'w') as f:
    f.write("Tracking memory usage for forward searches\n\n")
    f.write(f"{datetime.now()} - Current,Peak\n")
    f.close()
tracemalloc.start()
no_of_workers = 8
# Parallel(n_jobs=16)(delayed(sp.getoutput)(i) for i in cmd_list)
cmd_list = [cmd_list[i:i + 100] for i in range(0, len(cmd_list), 100)]
for i in cmd_list:
    pull_run(no_of_workers, i)    
    with open(os.path.join(path2runfolder, f"{run_folder}_forward_searches_memory_usage.txt"), 'a') as f:
        # old = f.read()
        f.write(f"\n{datetime.now()} - {round(tracemalloc.get_traced_memory()[0]/1e3, 4)} KB,{round(tracemalloc.get_traced_memory()[1]/1e3, 4)} KB\n")
        f.close()
del i, f

tracemalloc.stop()

print("\nForward searches completed\n")
#%%
incomplete_forward_searches = []
for i in os.listdir():
    if os.stat(i).st_size == 0:
        print(i)
        incomplete_forward_searches.append(i)
del i

if incomplete_forward_searches:
    os.chdir(path2runfolder)
    with open(f"forward_searches_log_{run_folder}.txt", 'w') as f:
        f.write("The following forward searches were not completed. Please check before proceeding\n\n")
        [f.write(f"{x}\n") for x in incomplete_forward_searches]
        f.close()
    del f

    if os.name == 'nt':
        cmd_list = ['wsl phmmer --tblout {0} --domtblout dom_{0} {1} {2}'.format(forward_searches_files_list[i], (os.path.join(path2queryseqs, queryfilelist[i])), (path2targetdb)) for i in range(len(forward_searches_files_list))]
    else:
        cmd_list = ['phmmer --tblout {0} --domtblout dom_{0} {1} {2}'.format(forward_searches_files_list[i], (os.path.join(path2queryseqs, queryfilelist[i])), (path2targetdb)) for i in range(len(forward_searches_files_list))]
    
    for i in incomplete_forward_searches:
        for j in range(len(cmd_list)):
            if i not in cmd_list[j]:
                cmd_list[j] = []
    del i, j
    
        
    cmd_list = list(filter(None, cmd_list))
    if cmd_list:
        for i in cmd_list:
            sp.getoutput(i)
        del i
    try:
        del linuxpath2queryseqs, linuxpath2targetdb, cmd_list
    except NameError:
        pass
else:
    del incomplete_forward_searches
#%%
print("Retrieving top hit for each species\n")    
os.chdir(path2forwardsearches)
best_hit_per_species = {}
hits_profile = pd.DataFrame(index=listofspecies, columns=[x[:-6] for x in queryfilelist])
for i in queryfilelist:
# i = forward_searches_files_list[0]
    temp = read_phmmer_tblout([x for x in forward_searches_files_list if i[:-6] in x][0])
    temp_best_hits = []
    for j in listofspecies:
        
        temp1 = temp[temp["description of target"].str.endswith(f"[{j}]")] 
        hits_profile[i[:-6]][j] = len(temp1)
        if len(temp1) == 0:
            # hits_profile[i[:-6]][j] = 0
            continue
        temp1 = temp1[temp1["best 1 domain score"] == max(list(temp1["best 1 domain score"]))]
        
        temp_best_hits.append(list(temp1["target name"])[0])
    best_hit_per_species[i[:-6]] = temp_best_hits
del i, j, temp, temp1, temp_best_hits, cmd_list
#%%
temp_seq = []
temp_seq_id = []

for i in SeqIO.parse(open(path2targetdb), 'fasta'):
    temp_seq.append(i)
    temp_seq_id.append(i.id)
del i
#% writing hits to new folder
os.chdir(path2runfolder)
#%%
forward_searches_hits = "top_hits_" + forward_searches
os.mkdir(forward_searches_hits)
os.chdir(forward_searches_hits)
path2forwardsearches_hits = os.getcwd()

for i in best_hit_per_species:
    os.mkdir(i)
    os.chdir(i)
    print(f"Writing best hits for {i}")
    # for j in best_hit_per_species[i]:
    inds = [temp_seq_id.index(x) for x in best_hit_per_species[i]]
    records_to_write = [temp_seq[x] for x in inds]
    p = 0
    for j in records_to_write:
        name = str(f"{p}_{j.id[:j.id.find('_')]}.fasta")
        SeqIO.write(j, name, 'fasta')
        p += 1
    os.chdir(path2forwardsearches_hits)
del i, p, j, name, inds, temp_seq, temp_seq_id
#%%
# os.chdir(os.path.dirname(os.getcwd()))
best_hit_per_species_paths = {}
for i in best_hit_per_species:
    os.chdir(i)
    if os.name == 'nt':
        best_hit_per_species_paths[i] = [winpath2linuxpath(os.path.join(os.getcwd(),x)) for x in os.listdir()]
    else:
        best_hit_per_species_paths[i] = [os.path.join(os.getcwd(),x) for x in os.listdir()]
    os.chdir(path2forwardsearches_hits)
del i

os.chdir(path2runfolder)
reverse_searches = "reverse_searches_{0}_{1}".format(iteration, run_folder)
os.mkdir(reverse_searches)
path2reversesearches = os.path.join(path2runfolder, reverse_searches)
os.chdir(path2reversesearches)

file_paths = {}
for i in best_hit_per_species:
    os.mkdir(i)
    if os.name == 'nt':
        file_paths[i] = winpath2linuxpath(os.path.join(os.getcwd(), i))
    else:
        file_paths[i] = os.path.join(os.getcwd(), i)
del i


# dill.dump_session(session_pickle)
#%
prefix = input("\nPlease enter a suitable prefix to add to the reverse search result files. This is to make it more straightforward to understand which database was the search done\n\n")
cmd_list = []
for i in best_hit_per_species_paths:
    if os.name == 'nt':
        # path_name = 
        name = {x[x.rfind("/")+1:-6] : x for x in best_hit_per_species_paths[i]}
        [cmd_list.append(f"wsl phmmer --tblout {file_paths[i]}/{x}_reverse_searches_{prefix}.txt --domtblout {file_paths[i]}/dom_{x}_reverse_searches_{prefix}.txt {name[x]} {winpath2linuxpath(path2sourcedb)}") for x in name]
    else:
        name = {x[x.rfind("/")+1:-6] : x for x in best_hit_per_species_paths[i]}
        [cmd_list.append(f"phmmer --tblout {file_paths[i]}/{x}_reverse_searches_{prefix}.txt --domtblout {file_paths[i]}/dom_{x}_reverse_searches_{prefix}.txt {name[x]} {path2sourcedb}") for x in name]
del i, name

tracemalloc.start()
with open(os.path.join(path2runfolder, f"{run_folder}_reverse_searches_memory_usage.txt"), 'w') as f:
    f.write("Tracking memory usage for reverse searches\n\n")
    f.write(f"{datetime.now()} - Current,Peak\n")
    f.close()
cmd_list = [cmd_list[i:i + 100] for i in range(0, len(cmd_list), 100)]
for i in cmd_list:
    pull_run(no_of_workers, i)    
    with open(os.path.join(path2runfolder, f"{run_folder}_reverse_searches_memory_usage.txt"), 'a') as f:
        # old = f.read()
        f.write(f"\n{datetime.now()} - {round(tracemalloc.get_traced_memory()[0]/1e3, 4)} KB,{round(tracemalloc.get_traced_memory()[1]/1e3, 4)} KB\n")
        f.close()
del i, f
tracemalloc.stop()
# pull_run(os.cpu_count(), cmd_list)
#%%
incomplete_searches = search_completenes_check(path2reversesearches)
iter1 = 0
while incomplete_searches:
    
    os.chdir(path2runfolder)
    
    for i in range(len(incomplete_searches)):
        incomplete_searches[i] = [x for x in cmd_list if incomplete_searches[i] in x][0]
    del i
    with open(f"completeness_check_{iter1}_log.txt", 'w') as f:
        f.write("The following commands failed to execute:\n\n")
        for i in incomplete_searches:
            f.write(f"{i}\n\n")
        f.write("The script has attempted to complete them. If you see the same command occuring in the log files. Please check\n")
        f.close()
    del i, f
#% reexecuting the failed commands
    for i in incomplete_searches:
        print(f"{i}\n\n")
        sp.getoutput(i)
    del i
    iter1 += 1
    incomplete_searches = search_completenes_check(path2reversesearches)
del iter1, cmd_list
#%%
print("\nSelecting best, bidirectional hits (BBHs)\n")
if os.name == 'nt':
    list_of_verified_hits = retrievebbh(path2reversesearches, serverpath2winpath(path2isodict))
else:
    list_of_verified_hits = retrievebbh(path2reversesearches, path2isodict)
#%%

verified_hits_src_path = {}
os.chdir(path2forwardsearches_hits)
for i in list_of_verified_hits:
    os.chdir(i)
    # if os.name == 'nt':
    #     verified_hits_src_path[i] = [winpath2linuxpath(os.path.join(os.getcwd(), x)) for x in list_of_verified_hits[i]]
    # else:
    verified_hits_src_path[i] = [os.path.join(os.getcwd(), x) for x in list_of_verified_hits[i]]
    os.chdir(path2forwardsearches_hits)
del i

os.chdir(path2runfolder)
verified_hits = "verified_hits_{0}".format(iteration)
os.mkdir(verified_hits)
path2verifiedhits_0 = os.path.join(path2runfolder, verified_hits)
os.chdir(path2verifiedhits_0)
verified_hits_ids = {}    
for i in verified_hits_src_path:
    os.mkdir(i)
    os.chdir(i)
    temp = []
    for j in verified_hits_src_path[i]:
        [temp.append(x) for x in SeqIO.parse(j, 'fasta')]
    with open(f"{i}_combined_hits_{run_folder}.fasta", 'w') as g:
        [SeqIO.write(x, g, 'fasta') for x in temp]
        g.close()
    verified_hits_ids[i] = temp
    os.chdir(path2verifiedhits_0)
del i, j, temp, g

#%%
for i in verified_hits_ids:
    verified_hits_ids[i] = [x.id for x in verified_hits_ids[i]]
del i
#%
print("\nWriting log files for the run\n")
os.chdir(path2runfolder)
os.mkdir("log_of_runs_for_{0}".format(run_folder))
os.chdir("log_of_runs_for_{0}".format(run_folder))
with open("list_of_species.txt", 'w') as f:
    f.write("The species used in this search are:\n")
    for i in listofspecies:
        f.write("{0}\n".format(i))
    f.close()
del f


with open("Map_of_created_files.txt", 'w') as f:
    f.write("All the commands used for phmmer in this run can be found in the results files\n")
    f.write("\nThe run folder itself can found at: {0} \n".format(path2runfolder))
    f.write("\nThe queries can found at: {0} \n".format(querylist))
    f.write("\nThe target database can be found at: {0} \n".format(path2targetdb))
    f.write("\nThe source database can be found at: {0} \n".format(path2sourcedb))
    f.write("\nThe isoform dictionary used can be found at: {0} \n".format(path2isodict))
    f.write("\nThe search results for the first forward searches can be found at: {0} \n".format(path2forwardsearches))
    f.write("\nThe top hits for each species from the first forward searches can be found at: {0} \n".format(path2forwardsearches_hits))
    f.write("\nThe search results for the reverse searches using the selected top hits can be found at: {0} \n".format(path2reversesearches))
    f.write("\nThe best, bidirectional hits (BBHs) selected from the top hits can be found at: {0} \n".format(path2verifiedhits_0))
    f.write("\nThe set of hits selected for the next iteration of searches can be found at: {0} \n".format(path2forwardsearches_hits))
del f

# if ids_flagged:
#     with open("Flagged_IDs.txt", 'w') as f:
#         f.write("These IDs were flagged because they returned no best, bidirectional hits. There can be a number of reasons for this. Examining the proteins referenced by these IDs along with closely related proteins will likely provide information for further analysis\n\n")
#         for i in ids_flagged:
#             f.write("{0}\n".format(i))
#         f.close()
#     del i, f

df_0 = pd.DataFrame(0, index=listofspecies, columns=[x[:-6] for x in queryfilelist])

for i in list_of_verified_hits:
    
    for j in list_of_verified_hits[i]:
        temp2 = j[j.find("_")+1:-6]
        df_0[i][speciesdict[temp2]] = df_0[i][speciesdict[temp2]] + 1
del i, j, temp2

df_0.to_csv("hits_per_species_0.csv")

with open("species_with_no_bbhs.txt", 'w') as f:
    for i in df_0.columns:        
        f.write("{0} - ".format(i))
        x = df_0.index[df_0[i] == 0].tolist()
        for j in x:
            f.write("{0},".format(j))
            
        f.write("\n")
    f.close()
del f, i, j, x

df_1 = pd.DataFrame(None, index=listofspecies, columns=[x[:-6] for x in queryfilelist])
for i in list_of_verified_hits:    
    for j in verified_hits_ids[i]:
        temp2 = j[:j.find("_")]
        if not type(df_1[i][speciesdict[temp2]]) == float:
            df_1[i][speciesdict[temp2]] += f",{j}"
        else:
            df_1[i][speciesdict[temp2]] = j
del i, temp2, j
df_1.to_csv("list_of_verified_hits.csv")
