# Master's Thesis: Computational Exploration Towards Universal Vaccines for Influenza and Coronaviruses
Menno Van Damme

This repository contains the code, data and results obtained during my thesis. Note that not all data is available here, the raw NCBI downloads, the SARS-CoV-2 fasta files and metadata json were to big to put on GitHub. The complete data is available on Zenodo: https://doi.org/10.5281/zenodo.6621827. If you only want to look at the results you don't really need these files.

If you want to run the entire pipelie completely from scratch, the python scripts for Influenza should be run in this order:
- Data.py
- Alignment.py
- Conservation.py
- Phylogentics.py
- Structure.py

Some of the SARS-CoV-2 steps had to be performed on the Ghent Univeristy supercomputer because of the large number of genome sequences and the length of the genome. These were the alignment, clustering and RNA secondary structure prediction. These are the .sh files in the SARSCoV2 directory.

The SARS-CoV-2 scripts should be run in this order:
- Data.py
- MSA.sh
- Conservation.py
- Cluster.sh
- Phylogenetics.py
- Structure.py & RNAfold.sh

The data direcotories are structured like this:

Influenza A & B:
Data
-  NCBI (not on GitHub, see Zenodo)
-  Sorted
-  Clean
-  Metadata
-  Aligned
-  Dataframes
-  Genome
-  Clusters
-  Tree
-  Structure

SARS-CoV-2:
Data
-  NCBI (not on GitHub, see Zenodo)
-  Clean (not on GitHub, see Zenodo)
-  Metadata
-  Alignment (not on GitHub, see Zenodo)
-  Dataframes
-  Genome
-  Clusters
-  Tree
-  Structure

<pre>
├── Aligned
│   ├── A_1.fasta
│   ├── A_2.fasta
│   ├── A_3.fasta
│   ├── A_4.fasta
│   ├── A_5.fasta
│   ├── A_6.fasta
│   ├── A_7.fasta
│   ├── A_8.fasta
│   ├── B_1.fasta
│   ├── B_2.fasta
│   ├── B_3.fasta
│   ├── B_4.fasta
│   ├── B_5.fasta
│   ├── B_6.fasta
│   ├── B_7.fasta
│   └── B_8.fasta
├── Clean
│   ├── A_1.fasta
│   ├── A_2.fasta
│   ├── A_3.fasta
│   ├── A_4.fasta
│   ├── A_5.fasta
│   ├── A_6.fasta
│   ├── A_7.fasta
│   ├── A_8.fasta
│   ├── B_1.fasta
│   ├── B_2.fasta
│   ├── B_3.fasta
│   ├── B_4.fasta
│   ├── B_5.fasta
│   ├── B_6.fasta
│   ├── B_7.fasta
│   └── B_8.fasta
├── Clusters
│   ├── A_1.fasta
│   ├── A_1.fasta.clstr
│   ├── A_2.fasta
│   ├── A_2.fasta.clstr
│   ├── A_3.fasta
│   ├── A_3.fasta.clstr
│   ├── A_4.fasta
│   ├── A_4.fasta.clstr
│   ├── A_5.fasta
│   ├── A_5.fasta.clstr
│   ├── A_6.fasta
│   ├── A_6.fasta.clstr
│   ├── A_7.fasta
│   ├── A_7.fasta.clstr
│   ├── A_8.fasta
│   ├── A_8.fasta.clstr
│   ├── B_1.fasta
│   ├── B_1.fasta.clstr
│   ├── B_2.fasta
│   ├── B_2.fasta.clstr
│   ├── B_3.fasta
│   ├── B_3.fasta.clstr
│   ├── B_4.fasta
│   ├── B_4.fasta.clstr
│   ├── B_5.fasta
│   ├── B_5.fasta.clstr
│   ├── B_6.fasta
│   ├── B_6.fasta.clstr
│   ├── B_7.fasta
│   ├── B_7.fasta.clstr
│   ├── B_8.fasta
│   ├── B_8.fasta.clstr
│   └── cluster_info.csv
├── Dataframes
│   ├── Complete
│   │   ├── df_A_1.csv
│   │   ├── df_A_2.csv
│   │   ├── df_A_3.csv
│   │   ├── df_A_4.csv
│   │   ├── df_A_5.csv
│   │   ├── df_A_6.csv
│   │   ├── df_A_7.csv
│   │   ├── df_A_8.csv
│   │   ├── df_B_1.csv
│   │   ├── df_B_2.csv
│   │   ├── df_B_3.csv
│   │   ├── df_B_4.csv
│   │   ├── df_B_5.csv
│   │   ├── df_B_6.csv
│   │   ├── df_B_7.csv
│   │   └── df_B_8.csv
│   ├── Scores
│   │   ├── df_A.csv
│   │   ├── df_A_1.csv
│   │   ├── df_A_2.csv
│   │   ├── df_A_3.csv
│   │   ├── df_A_4.csv
│   │   ├── df_A_5.csv
│   │   ├── df_A_6.csv
│   │   ├── df_A_7.csv
│   │   ├── df_A_8.csv
│   │   ├── df_B.csv
│   │   ├── df_B_1.csv
│   │   ├── df_B_2.csv
│   │   ├── df_B_3.csv
│   │   ├── df_B_4.csv
│   │   ├── df_B_5.csv
│   │   ├── df_B_6.csv
│   │   ├── df_B_7.csv
│   │   └── df_B_8.csv
│   ├── Ungapped
│   │   ├── df_A.csv
│   │   ├── df_A_1.csv
│   │   ├── df_A_2.csv
│   │   ├── df_A_3.csv
│   │   ├── df_A_4.csv
│   │   ├── df_A_5.csv
│   │   ├── df_A_6.csv
│   │   ├── df_A_7.csv
│   │   ├── df_A_8.csv
│   │   ├── df_B.csv
│   │   ├── df_B_1.csv
│   │   ├── df_B_2.csv
│   │   ├── df_B_3.csv
│   │   ├── df_B_4.csv
│   │   ├── df_B_5.csv
│   │   ├── df_B_6.csv
│   │   ├── df_B_7.csv
│   │   └── df_B_8.csv
│   ├── complete_score_df_A.csv
│   ├── complete_score_df_B.csv
│   ├── df_A.csv
│   ├── df_B.csv
│   ├── score_df_A.csv
│   ├── score_df_B.csv
│   ├── small_score_df_A_Mutability.csv
│   ├── small_score_df_A_Shannon_entropy.csv
│   ├── small_score_df_B_Mutability.csv
│   └── small_score_df_B_Shannon_entropy.csv
├── Genome
│   ├── Aligned
│   │   ├── A.fasta
│   │   ├── A_1.fasta
│   │   ├── A_2.fasta
│   │   ├── A_3.fasta
│   │   ├── A_4.fasta
│   │   ├── A_5.fasta
│   │   ├── A_6.fasta
│   │   ├── A_7.fasta
│   │   ├── A_8.fasta
│   │   ├── B.fasta
│   │   ├── B_1.fasta
│   │   ├── B_2.fasta
│   │   ├── B_3.fasta
│   │   ├── B_4.fasta
│   │   ├── B_5.fasta
│   │   ├── B_6.fasta
│   │   ├── B_7.fasta
│   │   └── B_8.fasta
│   ├── CDS
│   │   ├── A_1.fasta
│   │   ├── A_2.fasta
│   │   ├── A_3.fasta
│   │   ├── A_4.fasta
│   │   ├── A_5.fasta
│   │   ├── A_6.fasta
│   │   ├── A_7.fasta
│   │   ├── A_8.fasta
│   │   ├── A_CDS.fasta
│   │   ├── B_1.fasta
│   │   ├── B_2.fasta
│   │   ├── B_3.fasta
│   │   ├── B_4.fasta
│   │   ├── B_5.fasta
│   │   ├── B_6.fasta
│   │   ├── B_7.fasta
│   │   ├── B_8.fasta
│   │   └── B_CDS.fasta
│   ├── CDS_df_A.csv
│   ├── CDS_df_B.csv
│   └── Consensus
│       ├── A_1.fasta
│       ├── A_2.fasta
│       ├── A_3.fasta
│       ├── A_4.fasta
│       ├── A_5.fasta
│       ├── A_6.fasta
│       ├── A_7.fasta
│       ├── A_8.fasta
│       ├── A_consensus.fasta
│       ├── B_1.fasta
│       ├── B_2.fasta
│       ├── B_3.fasta
│       ├── B_4.fasta
│       ├── B_5.fasta
│       ├── B_6.fasta
│       ├── B_7.fasta
│       ├── B_8.fasta
│       └── B_consensus.fasta
├── Metadata
│   ├── metadata.json
│   ├── metadata_clean.json
│   ├── metadata_cluster.json
│   └── metadata_complete.json
├── NCBI
│   ├── influenza.fna
│   └── influenza_na.dat
├── Sorted
│   ├── A_1.fasta
│   ├── A_2.fasta
│   ├── A_3.fasta
│   ├── A_4.fasta
│   ├── A_5.fasta
│   ├── A_6.fasta
│   ├── A_7.fasta
│   ├── A_8.fasta
│   ├── B_1.fasta
│   ├── B_2.fasta
│   ├── B_3.fasta
│   ├── B_4.fasta
│   ├── B_5.fasta
│   ├── B_6.fasta
│   ├── B_7.fasta
│   └── B_8.fasta
├── Structure
│   ├── Colors
│   │   ├── A.txt
│   │   ├── A_1.txt
│   │   ├── A_2.txt
│   │   ├── A_3.txt
│   │   ├── A_4.txt
│   │   ├── A_5.txt
│   │   ├── A_6.txt
│   │   ├── A_7.txt
│   │   ├── A_8.txt
│   │   ├── B.txt
│   │   ├── B_1.txt
│   │   ├── B_2.txt
│   │   ├── B_3.txt
│   │   ├── B_4.txt
│   │   ├── B_5.txt
│   │   ├── B_6.txt
│   │   ├── B_7.txt
│   │   └── B_8.txt
│   ├── Consensus
│   │   ├── A.fasta
│   │   ├── A_1.fasta
│   │   ├── A_2.fasta
│   │   ├── A_3.fasta
│   │   ├── A_4.fasta
│   │   ├── A_5.fasta
│   │   ├── A_6.fasta
│   │   ├── A_7.fasta
│   │   ├── A_8.fasta
│   │   ├── B.fasta
│   │   ├── B_1.fasta
│   │   ├── B_2.fasta
│   │   ├── B_3.fasta
│   │   ├── B_4.fasta
│   │   ├── B_5.fasta
│   │   ├── B_6.fasta
│   │   ├── B_7.fasta
│   │   └── B_8.fasta
│   ├── Dotbracket
│   │   ├── A.txt
│   │   ├── A_1.ifold
│   │   ├── A_1.txt
│   │   ├── A_2.ifold
│   │   ├── A_2.txt
│   │   ├── A_3.ifold
│   │   ├── A_3.txt
│   │   ├── A_4.ifold
│   │   ├── A_4.txt
│   │   ├── A_5.ifold
│   │   ├── A_5.txt
│   │   ├── A_6.ifold
│   │   ├── A_6.txt
│   │   ├── A_7.ifold
│   │   ├── A_7.txt
│   │   ├── A_8.ifold
│   │   ├── A_8.txt
│   │   ├── B.txt
│   │   ├── B_1.ifold
│   │   ├── B_1.txt
│   │   ├── B_2.ifold
│   │   ├── B_2.txt
│   │   ├── B_3.ifold
│   │   ├── B_3.txt
│   │   ├── B_4.ifold
│   │   ├── B_4.txt
│   │   ├── B_5.ifold
│   │   ├── B_5.txt
│   │   ├── B_6.ifold
│   │   ├── B_6.txt
│   │   ├── B_7.ifold
│   │   ├── B_7.txt
│   │   ├── B_8.ifold
│   │   └── B_8.txt
│   └── Element
│       ├── A_1_001.element_string
│       ├── A_2_001.element_string
│       ├── A_3_001.element_string
│       ├── A_4_001.element_string
│       ├── A_5_001.element_string
│       ├── A_6_001.element_string
│       ├── A_7_001.element_string
│       ├── A_8_001.element_string
│       ├── B_1_001.element_string
│       ├── B_2_001.element_string
│       ├── B_3_001.element_string
│       ├── B_4_001.element_string
│       ├── B_5_001.element_string
│       ├── B_6_001.element_string
│       ├── B_7_001.element_string
│       └── B_8_001.element_string
└── Tree
    ├── A_1.fasta
    ├── A_1.tree
    ├── A_2.fasta
    ├── A_2.tree
    ├── A_3.fasta
    ├── A_3.tree
    ├── A_4.fasta
    ├── A_4.tree
    ├── A_5.fasta
    ├── A_5.tree
    ├── A_6.fasta
    ├── A_6.tree
    ├── A_7.fasta
    ├── A_7.tree
    ├── A_8.fasta
    ├── A_8.tree
    ├── B_1.fasta
    ├── B_1.tree
    ├── B_2.fasta
    ├── B_2.tree
    ├── B_3.fasta
    ├── B_3.tree
    ├── B_4.fasta
    ├── B_4.tree
    ├── B_5.fasta
    ├── B_5.tree
    ├── B_6.fasta
    ├── B_6.tree
    ├── B_7.fasta
    ├── B_7.tree
    ├── B_8.fasta
    └── B_8.tree
</pre>




