# Nematode-SuccinateDehydrogenase-HMMer-Tool
## Succinate Dehydrogenase Analysis in Nematodes

## Overview
This project is a bioinformatics tool designed to automate the analysis of **Succinate Dehydrogenase** in nematodes. The tool leverages data from **WormBase** and the **InterPro** database to perform HMMer searches on the **University of Leicester's HPC (ALICE)**. It streamlines the process of downloading relevant data, running computational searches, and generating informative visual outputs.

## Features
- **Downloads and unzips FASTA files** for three chosen nematode species from WormBase.
- **Extracts PFAM identifiers** related to Succinate Dehydrogenase from a TSV file.
- **Downloads and unzips** the corresponding HMM profiles.
- **Generates a shell script** to run HMMer searches on the ALICE HPC system.
- **Parses HMMer outputs** to generate:
  - A detailed **summary table of hits**.
  - A **heatmap of hit scores**.
  - A **bar chart of the top 10 hits**.

## Requirements

### Python Libraries:
- pandas
- subprocess
- urllib
- BioPython
- requests
- BeautifulSoup
- matplotlib
- seaborn

### Internet Connection:
- Required for downloading data.

### University of Leicester HPC Access:
- Required to run the generated HMMer search script.

### TSV File:
- **SearchResults-succinatedehydrogenase.tsv**: This file should be in the working directory (optional file name can be changed).

## Usage Instructions

### Step 1: Downloading and Unzipping FASTA Files
1. The program will list all available species from WormBase.
2. You will be prompted to select **three species** by their index numbers, separated by spaces.
3. The selected **FASTA files** will be downloaded and unzipped automatically.

### Step 2: Extracting and Downloading Pfam HMM Profiles
1. Ensure that the **SearchResults-succinatedehydrogenase.tsv** file is in your working directory.
2. You will be prompted to type 'y' to confirm the file's presence, or type 'change' to specify a different file name or path.
3. The program will extract **PFAM identifiers** (starting with 'PF') from the TSV file and download the corresponding **HMM profiles**.

### Step 3: Generating HPC (ALICE) HMMer Run Script
1. You will be prompted to enter your **email address**.
2. The program will generate a **shell submission script (HMMsearch.sh)** for the **SLURM scheduler** on ALICE. This script contains the HMMer commands to run the searches.

### Step 4: Parsing and Analyzing HMMer Outputs
1. Run the generated script on **ALICE** to obtain the **HMMer output files**.
2. After generating the output files, you will be prompted to type 'y' to confirm their presence in the current directory.
3. The program will parse the results and produce:
   - A detailed **table of all hits** (**hmmer_output_summary.csv**).
   - A **heatmap of the scores** (**hmmer_output_heatmap.png**).
   - A **bar chart of the top 10 hits** (**hmmer_top_hits_bar_chart.png**).

## Output Files
- **hmmer_output_summary.csv**: A detailed table of HMMer hits.
- **hmmer_output_heatmap.png**: A heatmap visualizing the scores.
- **hmmer_top_hits_bar_chart.png**: A bar chart showing the top 10 hits.

## Important Notes
- The program is **not case sensitive**.
- Ensure that the program and the **HPC script** are run in the **same directory**.
- All generated files will be saved in the **current directory**.
- Ensure all required libraries are installed in your **Python environment**.
- The program requires a stable **Internet connection**.

## Contributing
Contributions are welcome! Please **fork** the repository and **submit a pull request**.

## Contact
For any questions or issues, please contact **Bismah Ghafoor** at **Bismahghafoor470@gmail.com**.
