The code in this repository was used to generate the data used for the review - Origin and Evolution of Microvilli. 
This README is a quick introduction to the code. For information on how exactly to use it, please reach out!

## Data availability
The proteomes used for the analysis, and the results can found on - https://figshare.com/articles/online_resource/Origin_and_evolution_of_microvilli/25828882/2

## Pre-analysis work
The header for every sequence in the proteomes were formatted - >abcd_protein_name [Species name] where "abcd" is an abbreviation unique to each species

Since the UniProt human proteome was used as the reference for performing the best, bidirectional hits searches, an isoform dictionary was built. The proteome file and the dictionary can be found in the folder labelled "human_proteome" in this repo. The Uniprot IDs for the query proteins were organized into a CSV file which can also be found in the same folder
