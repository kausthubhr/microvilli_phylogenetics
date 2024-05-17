# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:06:12 2024

@author: Kausthubh R
"""


import os
import subprocess as sp
from Bio import SeqIO

trimv = sp.getoutput("trimal --version")
folder_path = os.getcwd()

align_check = [x for x in os.listdir() if x.endswith("_trim.phy")]
if align_check:
    with open("error.txt", 'w') as f:
        f.write("this folder already seems to contain alignments. Please check")
        f.close()
    exit()
    
name_list = [f for f in os.listdir() if f.endswith('interpro.fasta')]

if not name_list:
    with open("error.txt", 'w') as f:
        f.write("this folder does not seem to have any fasta files. Please check")
        f.close()
    exit()

name = name_list[0]
if not name:
    exit()
name_out = name[:-6] + "_align.phy"
print(name)
temp_record = []
for record in SeqIO.parse(open(name), 'fasta'):
    temp_record.append(record)
del record

# =============================================================================
# Exiting if the number of sequences available is less than 4
# =============================================================================

if len(temp_record) <= 4:
    with open("{0}_record.txt".format(name), 'w') as f:
        f.write("\nThe number of available sequences is less than 4. It is therefore meaningless to build alignments / trees with bootstrap\n")
    exit()


# =============================================================================
# Checking for human isoforms and fragments as annotated in the UniProt 
# reference proteome
# =============================================================================
human_check = []
for i in range(len(temp_record)):
    if temp_record[i]:
        if temp_record[i].id.startswith("Hsap"):
            if "Isoform" in temp_record[i].description:
                human_check.append(temp_record[i])
                temp_record[i] = []
                continue
            if "isoform" in temp_record[i].description:
                human_check.append(temp_record[i])
                temp_record[i] = []
                continue
            if "Fragment" in temp_record[i].description:
                human_check.append(temp_record[i])
                temp_record[i] = []
                continue
            if "fragment" in temp_record[i].description:
                human_check.append(temp_record[i])
                temp_record[i] = []
                continue           
del i

if human_check:
    with open(f"{name[:name.rfind('.')]}_human_isoforms_fragments.fasta", 'w') as f:
        [SeqIO.write(x, f, 'fasta') for x in human_check]
        f.close()

temp_record = list(filter(None, temp_record))
#%
# identical = []
# abbrevs = []
# for i in temp_record:
#     abbrevs.append(i.id[:i.id.find("_")])
# del i

# abbrevs = list(set(abbrevs))
# for i in temp_record:
#     other_temp = [x for x in temp_record if x.id != i.id and x.id[:x.id.find("_")]]
#     if any([x for x in other_temp if x.seq == i.seq]):
#         identical.append(i)
# del i, other_temp
#%

temp_record_len = []
for i in temp_record:
    temp_record_len.append(len(i))
del i


# =============================================================================
# Removing sequences that are smaller than 10% of the average length of the 
# sequences in the dataset
# =============================================================================
avg = sum(temp_record_len)/(len(temp_record_len))
print(avg)
tiny_seqs = []
y = 0.1
for i in range(len(temp_record)):
    if len(temp_record[i]) < y*avg:
        tiny_seqs.append(temp_record[i])
        temp_record[i] = []
del i
temp_record = list(filter(None, temp_record))

if tiny_seqs:
    with open(f"{name[:name.rfind('.')]}_small_seqs.fasta", "w") as f:
        [SeqIO.write(x, f, 'fasta') for x in tiny_seqs]
        f.close()


name_updated = f"{name[:name.rfind('.')]}_updated.fasta"
# name_updated = name


with open((name_updated), 'w') as f:
    for i in temp_record:
        SeqIO.write(i, f, 'fasta')
    f. close()
del f, i

updated_len = len(temp_record)
print("Aligning the updated {0} sequences\n".format(name))

name_out = name_updated[:name_updated.rfind(".")] + "_align.phy"

# =============================================================================
# Exiting if the number of sequences available is less than 4
# =============================================================================

if len(temp_record) <= 4:
    with open("{0}_record.txt".format(name), 'w') as f:
        f.write("\nThe number of available sequences is less than 4. It is therefore meaningless to build alignments / trees with bootstrap\n")
    exit()


#mafft
if len(temp_record) <= 500:
    mafft_cmd = ("mafft --maxiterate 1000 --bl 45 --genafpair --reorder {0} > {1}").format(name_updated, name_out)
if len(temp_record) >= 500:
    mafft_cmd = ("mafft --retree 2 --maxiterate 1000 --bl 45 --reorder {0} > {1}").format(name_updated, name_out)

x = sp.getoutput(mafft_cmd) 

with open("mafft_record_{0}.txt".format(name_out[:-4]), 'w') as f:
    f.write(x)
    f.close()
del f, x

print("Trimming the genafpair aligned {0} sequences\n".format(name))
name_trim = name_out[:-4]+"_trim.phy"
# trim_cmd = "bmge -i {0} -opp {1} -m BLOSUM45 -t AA -oh {1}_log.html".format(name_out, name_trim)
trim_cmd = "trimal -in {0} -out {1} -gappyout -htmlout {1}_log.html".format(name_out, name_trim)
sp.getoutput(trim_cmd)

# =============================================================================
# Exiting if the number of sequences available is less than 4
# =============================================================================
if len(temp_record) < 4:
    with open("{0}_record.txt".format(name), 'w') as f:
        f.write("\nThe number of available sequences is less than 4. It is therefore meaningless to try and build trees with bootstrap\n")
    exit()


print("Building trees using gappyout trimmed alignments of the {0} sequences\n\n\n".format(name))
# tree_cmd = ('iqtree2 -m MFP -madd LG+C20,LG+C40,LG+C60 --seqtype AA -B 1000 --boot-trees --wbtl --score-diff ALL -mwopt -s {0}'.format(name_trim))
tree_cmd = ('FastTree -spr 4 -mlacc 2 -slownni -n 1000 -gamma -log {0}_fasttree_log.txt {0} > {0}_fasttree.tree'.format(name_trim))
# tree_cmd = ('iqtree2 -m LG+G --seqtype AA -B 1000 --boot-trees --wbtl --score-diff ALL -s {0}'.format(name_trim))
sp.getoutput(tree_cmd)

with open("{0}_record.txt".format(name), 'w') as f:
    f.write("Record of commands used for aligning, trimming, and building trees for {0}\n\n\n".format(name))
    
    f.write("Alignment tool used: MAFFT (version specified in tool's log file mafft_record_{0}.txt)\n".format(name_out[:-4]))
    f.write("Alignment tool command used: {0}\n".format(mafft_cmd))
    
    f.write("Aligment trimming tool used: TrimAl (version specified in tool's HTML output) \n")
    f.write("Aligment trimming tool used: {0}) \n".format(trimv))
    f.write("Alignment trimming tool command used: {0}\n".format(trim_cmd))
    f.write("Tree building tool used: FastTree (version specified in tool's log file)\n")
    f.write("Tree building tool command used: {0}\n".format(tree_cmd))
    f.write("Average protein length {0}\n\n\n".format(avg))
    f.write('Proteins that were smaller than {0} of the average protein length were removed from the alignment\n\n\n'.format(y))
    
    if human_check:
        f.write("The following sequences were removed prior to building alignments as the proteins were annotated to be either an isoform or a fragment in the Uniprot Homo sapiens reference proteome: \n\n")
        [f.write(f"{x.id}\n") for x in human_check]
    
    
    if tiny_seqs:
        f.write(f"The following sequences were removed prior to building alignments as the proteins were smaller than {y} times the average length ({avg} amino acids) of the proteins in the set: \n\n")
        [f.write(f"{x.id}\n") for x in tiny_seqs]
    
    f.close()
