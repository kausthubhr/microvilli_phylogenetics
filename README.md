The code in this repository was used to generate the data used for the review - Origin and Evolution of Microvilli. 
This README is a quick introduction to the code. For information on how exactly to use it, please reach out!

## Data availability
The proteomes used for the analysis, and the results can found on [Figshare](https://figshare.com/articles/online_resource/Origin_and_evolution_of_microvilli/25828882/2)

## Pre-analysis work
The header for every sequence in the proteomes were formatted - 
```
>abcd_protein_name [Species name]
```
where "abcd" is an abbreviation unique to each species. For example, "Hsap" would be the abbreviation for "Homo sapiens". Please take care to ensure that the abbreviations are unique to each species.

Since the UniProt human proteome was used as the reference for performing the best, bidirectional hits searches, an isoform dictionary was built. The proteome file and the dictionary can be found in the folder labelled "human_proteome" in this repo. The Uniprot IDs for the query proteins were organized into a CSV file which can also be found in the same folder.

## Setting up the conda environments
Two conda environments need to be set up. To install conda (Miniconda was used for this project), please refer to the [documentation](https://docs.anaconda.com/free/miniconda/index.html)
- For performing the BUSCO v5.7.1 analysis 
  ```bash
  conda env create -f busco_scoring_specs.yml
  ```
- For performing the remaining analysis
  ```bash
  conda env create -f microvilli_phylogenetics_specs.yml
  ```

## Running the scripts

### 1. Performing the BUSCO v5.7.1 analysis
1. Create the first conda environment "busco_scoring" as described above
2. Organize all your proteomes into a single directory
3. Activate the conda environment
   ```bash
   conda activate busco_scoring
   ```
4. Run BUSCO v5.7.1 using this command. This will perform the analysis using the eukaryota_odb10 lineage. For more information on BUSCO and lineages, please refer to the [BUSCO user guide](https://busco.ezlab.org/busco_userguide.html)
   ```bash
   busco -m prot -l eukaryota -i </path/to/directory/containing/proteomes> -o <name_of_output_directory> 
   ```
  >This analysis can take some time, depending on the number of proteomes used for the analysis. It took around ~5 hours for analysing 105 proteomes
5. The batch_summary.txt is the ouput (along with other stuff) 

For the remaining analysis, please ensure that you have set up the second conda environment and that it is active. The environment can be activated using
```bash
conda activate microvilli_phylogenetics
```

### 2. Identifying best, bidirectional hits for the query sequences
Please run the script for [retrieving hits](https://github.com/kausthubhr/microvilli_phylogenetics/blob/e19201428b6425a960b6fc9c519dc54d5e71df83/scripts/phmmer_based_bbh_retrieval_new_updated.py). All the required steps is detailed in the script itself.

This script uses phmmer from the HMMER suite of tools in the environment to identify the best hit in each proteome that maps back to the original query protein in the human proteome. The searches are parallelized and can use a significant amount of memory

### 3. Building alignments and gene trees
Please set your working directory to the "verified_hits_0" folder that should have been created after running the script in step 2. Within this directory, there will be one directory for each query protein. You will need to navigate to each directory and run the [tree building script](https://github.com/kausthubhr/microvilli_phylogenetics/blob/e19201428b6425a960b6fc9c519dc54d5e71df83/scripts/align_trim_tree_builder_new.py) there. 

This script uses 
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) for building alignments
- [TrimAl](https://vicfero.github.io/trimal/) in the "gappyout" mode for trimming the alignments
- [FastTree](http://www.microbesonline.org/fasttree/) or [IQ-TREE](http://www.iqtree.org/) for building maximum-likelihood phylogenetic trees
  - IQ-TREE is very much slower and computationall)y intensive when compared to FastTree. You can decide which algorithm to use by checking line number 179/180 in the tree building script. Make sure that the algorithm you want to use is not commented i.e. starts with a "#"
 
The phylogenetic trees and alignments were manually analysed prior to continuing the analysis. The trees were visualized using either [FigTree](https://github.com/rambaut/figtree) or [iTOL](https://itol.embl.de/).

### 4. Building Hidden Markov Models (HMMs)
Please set your working directory to the "verified_hits_0" folder that should have been created after running the script in step 2. Please run the [hmmsearch script](https://github.com/kausthubhr/microvilli_phylogenetics/blob/e19201428b6425a960b6fc9c519dc54d5e71df83/scripts/hmmsearches_first_iteration.py). Further instructions are available in the script itself.

### 5. Grouping the HMM search results into orthogroups
>Orthogroups for the proteins we used were taken from the literature. 

Please set your working directory to the "verified_hits_0" folder that should have been created after running the script in step 2. Please run the [hmmsearch sequence retrieval script](https://github.com/kausthubhr/microvilli_phylogenetics/blob/e19201428b6425a960b6fc9c519dc54d5e71df83/scripts/microvili_phylo_hmmsearch_seq_selection.py)

### 6. Performing the Interproscan analysis
#### Installing and running Interproscan
1. Please download and install Interproscan as described in the [documentation](https://interproscan-docs.readthedocs.io/en/latest/index.html)
2. Please organize the sequences you wish to analyse using Interproscan into a separate folder, with one folder per orthogroup
   >For this analysis, a folder was created within each orthogroup folder labelled "interpro". It is recommended that this folder organization be followed
4. Navigate to each folder, and run Interproscan using the following command
   ```bash
   /path/to/interproscan-5.67-99.0/interproscan.sh -i <fasta_file> -appl SMART,PANTHER,SUPERFAMILY,AntiFam -f tsv -vtsv
   ```
5. Once you have the Interproscan results, please run the [Interproscan analysis script](https://github.com/kausthubhr/microvilli_phylogenetics/blob/e19201428b6425a960b6fc9c519dc54d5e71df83/scripts/interproscan_analysis.py)

Following this step, there should be a final set of sequences for each orthogroup. Analyse the sequences in each orthogroup by repeating step 3. Further details on the final analysis are available in the manuscript.


### OPTIONAL: Adding new species to this analysis
If you wish to add new species to this analysis without redoing the analysis from scratch, please follow the following steps **prior to running the Interproscan analysis**
1. Organize the proteomes for the new species into a separate directory. Please include a few species from your original analysis in this directory.
   >NOTE - this directory should be within the run folder created in the initial analysis
2. Repeat steps 1 and 2
3. Run the [script for adding new species](https://github.com/kausthubhr/microvilli_phylogenetics/blob/e19201428b6425a960b6fc9c519dc54d5e71df83/scripts/microvili_phylo_hmmsearch_adding_new_species.py)

This script is still being developed. Please reach out if you are running into issues while running this script.
