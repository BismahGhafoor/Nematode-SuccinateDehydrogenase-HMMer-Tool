# This program is a specialized tool designed for bioinformatics research, specifically focusing on studying Succinate Dehydrogenase in nematodes. Firstly, it downloads and unzips FASTA data for three chosen nematode species from a list given to the user taken from the latest version of WormBase: https://parasite.wormbase.org/ftp.html . And then it extracts from a specified TSV file to PFAM identifiers related to Succinate Dehydrogenase (sourced from the InterPro database) and downloads and unzips them. After that, it generates a shell submission script for the SLURM scheduler on ALICE (the University of Leicester HPC) to run HMMer searches using the downloaded HMMs against the downloaded FASTA files. Finally, it parses the HMMer results, generating a detailed table of all of the hits, a score heatmap of all of the hits, and a bar chart of the top 10 score hits.
# Last Modification: 14 Dec 2023

import pandas as pd
import subprocess
from urllib import request
from Bio import SearchIO # BioPython used tutorial for this program: https://biopython.org/DIST/docs/tutorial/Tutorial.html
import re
import requests
from bs4 import BeautifulSoup   # BeautifulSoup used tutorial for this code: https://www.geeksforgeeks.org/how-to-scrape-websites-with-beautifulsoup-and-python/
import matplotlib.pyplot as plt # Learned from matplotlib tutorial: https://www.geeksforgeeks.org/matplotlib-tutorial/
import seaborn as sns   # Learned from seaborn tutorial: https://www.geeksforgeeks.org/python-seaborn-tutorial/



# Function 1: It downloads and unzippes a .fa.gz FASTA file of a given url that is from the wormbase website.
def download_and_unzip_fasta (url):
    #file name will be the same name given in the URL after the last "/"
    local_file = url.split('/')[-1]

    # To download the file.
    request.urlretrieve(url, local_file)

    # To unzip the file using a bash command.
    subprocess.call('gunzip ' + local_file, shell=True)
    
    # fucntion reutrns generated file name, since it is unzipped it will return it without ".gz" the last 3 letters
    return(local_file[:-3])



# Function 2: Extract Pfam HMM profiles from Accession column in a tsv file, downloads them, and unzippes them. After that it returns the created hmm file names.
def extract_and_download_Pfam (tsv):
    # Read the tsv file into a DataFrame
    df = pd.read_csv(tsv, sep='\t')

    # Extract Pfam identifiers that start with 'PF' from Accession column
    pfam_identifiers = [identifier for identifier in df['Accession'] if identifier.startswith('PF')]

    for pf_id in pfam_identifiers :
        # Download and unzips the HMM files using bash commands
        subprocess.call('wget https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/' + pf_id + '?annotation=hmm -O ' + pf_id + '.hmm.gz', shell=True)
        subprocess.call('gunzip ' + pf_id + '.hmm.gz', shell=True)

    # fucntion reutrns a list of the generated hmm file names
    return pfam_identifiers



# Function 3: generates a shell submission script for the SLURM scheduler on ALICE (the University of Leicester HPC) to run HMMer searches
def create_hmmer_alice_script(email, hmms, fastas) :    
    with open("HMMsearch.sh", 'w') as file:
        
        # SLURM directives
        file.write("#!/bin/bash\n")
        file.write("#SBATCH --job-name=HMMer_Nematodes\n")
        file.write("#SBATCH --nodes=1\n")
        file.write("#SBATCH --tasks-per-node=1\n")
        file.write("#SBATCH --mem=8gb\n")
        file.write("#SBATCH --time=02:00:00\n")
        file.write("#SBATCH --mail-type=BEGIN,END,FAIL\n")
        file.write("#SBATCH --mail-user=" + email + "\n\n")
        
        # Script to run the HMMer searches
        file.write("executable from ALICE")
        file.write("hmmsearch=/cm/shared/spack/opt/spack/linux-rocky9-x86_64_v3/gcc-12.3.0/hmmer-3.3.2-ipmjfm2vvzhroirpnpn5i4rw5wptqf7r/bin/hmmsearch\n\n")
        file.write("# Module needed for using HPC installed software\n")
        file.write("module load gcc/12.3.0-yxgv2bl\n")
        file.write("module load openmpi/4.1.5-fzc7xdf\n")
        file.write("module load hmmer/3.3.2-ipmjfm2\n\n")

        # All files are assumed to be in the current directory of the script
        file.write("# HMM and FASTA files (assumed to be the current directory)\n")
        file.write("hmm_dir=$(pwd)\n")
        file.write("fasta_dir=$(pwd)\n")
        file.write("output_dir=$(pwd)\n\n")

        # To write the HMMer commands
        for hmm in hmms:
            for fasta in fastas:
                file.write("hmmsearch --tblout ${output_dir}/" + hmm + "_" + fasta + ".out -E 0.1 --noali ${hmm_dir}/" + hmm + ".hmm ${fasta_dir}/" + fasta + "\n")
    print("\nALICE script (HMMsearch.sh) created successfully in the current directory\n")



#Function 4: it parses the HMMer results, generating a detailed table of all of the hits, a score heatmap of all of the hits, and a bar chart of the top 10 score hits.
def parse_analyse_hmmer_outputs(hmms, fastas):
    # Parse the HMMer output files
    extracted_data = []
    for hmm in hmms:
        for fasta in fastas:
            file_path = hmm + "_" + fasta + ".out" # This is how the output files were named as in the generated shell script (Function 3)
            hmmer_results = SearchIO.parse(file_path, "hmmer3-tab") # All output files are assumed to be in the current directory
            
            for query_result in hmmer_results:
                for hit in query_result.hits:
                    data = {
                        'target_name': hit.id,
                        'query_name' : hit.query_id,
                        'e_value': hit.evalue,
                        'score': hit.bitscore
                    }
                    extracted_data.append(data)
    
    # Convert the extracted data list to a pandas DataFrame
    df = pd.DataFrame(extracted_data)
    print(df) # show table to the user in the terminal
    
    # To export the DataFrame to a CSV file, generating a detailed table of all of the hits
    csv_file = "hmmer_output_summary.csv"
    df.to_csv(csv_file, sep='\t', index=False)

    # Create a pivot table for the heatmap
    pivot_table = df.pivot(index='target_name', columns='query_name', values='score')

    # Plotting the heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(pivot_table, annot=False, fmt=".1f", cmap="YlGnBu") # To show values change annot to True, but it would be messy with a lot of hits
    plt.title("HMMer Output Heatmap")
    plt.ylabel("Target Name")
    plt.xlabel("Query Name")
    plt.tight_layout()

    # Save the heatmap and show it to the user
    plt.savefig("hmmer_output_heatmap.png")
    plt.show()

    # Extract top 10 hits based on score
    top_hits = df.nlargest(10, 'score')

    # Plotting a bar chart for the top 10 score hits
    plt.figure(figsize=(10, 6))
    sns.barplot(x='score', y='target_name', data=top_hits, palette="viridis")
    plt.title("Top 10 Hits by Score")
    plt.xlabel("Score")
    plt.ylabel("Target Name")
    plt.tight_layout()
    
    # Save the bar chart and show it to the user
    plt.savefig("hmmer_top_hits_bar_chart.png")
    plt.show()





# Firstly, we need to extract all URLs that ends with ".protein.fa.gz" , from the latest version of WormBase:https://parasite.wormbase.org/ftp.html .
response = requests.get("https://parasite.wormbase.org/ftp.html")
soup = BeautifulSoup(response.text, 'html.parser')

urls = []
for link in soup.find_all('a'):
    href = link.get('href')
    if href and href.endswith('.protein.fa.gz'):
        urls.append(href)



# input 1: To list all speices names and realted BioProject IDs and to prompt the user to input 3 index numbers from the list, and then using Function 1 to download and unzip the 3 fasta files for them.
fastas = [] #to save the fasta file names
while True:
    # First it will list all of the Species Names with its BioProject IDs for the User.
    i=0 # 'i' will be the total number of URLs after the loop, it will be used later.
    for url in urls:
        BioProject = url.split('/')[-2] # The bioproject ID is written inside the WormBase URL
        Specie_name = url.split('/')[-3].replace("_", " ").capitalize() # The specie name is written inside the WormBase URL
        print(str(i) + "\t\t" + Specie_name + "\t\t" + BioProject)
        i=i+1

    # Ask the user to choose 3 speices index numbers from the list
    user_input = input("Enter three index numbers (from the list) separated by space to download their corresponding FASTA files: \n")
    index_ids = user_input.split()

    # this is to check if user did input three correct numbers sperated by space and to check if the number is within the range of urls list ('i' is the max in the range)
    if len(index_ids) == 3 and all(0 <= int(index_id) < i for index_id in index_ids):
        for index_id in index_ids:
            # To download and unzip fasta file using function 1, also it will save the generated fasta file name in fastas list.
            fastas.append(download_and_unzip_fasta(urls[int(index_id)]))
        print("\nFASTA files downloaded successfully in the current directory\n")
        break
    else:
        print("Invalid input. Please enter three valid index numbers.")



# input 2: To extract the  all PFAM identifiers (Accessions starting with  PF) from the tsv file, and to save and extract the HMM files. User will be asked to confirm tsv file name first.
while True:
    user_input2 = input("If you have the file 'SearchResults-succinatedehydrogenase.tsv' in the current directory, type 'y'. Otherwise, type 'change' to change the name of the TSV file: \n").lower() # .lower to make the user input here not case sensitive
    if user_input2 == 'y': 
        hmms = extract_and_download_Pfam('SearchResults-succinatedehydrogenase.tsv')
        break
    elif user_input2 == 'change':
        user_input2B = input("Input tsv file name: ")
        hmms = extract_and_download_Pfam(user_input2B)
        break
    else:
        print("Invalid input. Please type 'y' or 'change'")
print("\nHMM files downloaded successfully in the current directory\n")



# input 3: Asks the user for the email for the shell submission script for the ALICE HPC system.
user_email = input("To generate a shell submission script for the ALICE HPC system. Please enter your email: \n")
create_hmmer_alice_script(user_email, hmms, fastas)



# input 4: It parses the HMMer results, generating a detailed table of all of the hits, a score heatmap of all of the hits, and a bar chart of the top 10 score hits. But first it will ask the user to confirm getting the output files in the current directory.
while True:
    print("\nRun the script on ALICE to get the output files. The HMMer outputs must be in the current directory")
    user_input4 = input("Please type 'y' once you have generated the HMMer outputs: \n").lower()
    if user_input4 == 'y':
        parse_analyse_hmmer_outputs(hmms, fastas)
        break
