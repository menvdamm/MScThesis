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

Influenza A & B:
For each type (A or B) and each segment (1, 2, 3, 4, 5, 6, 7, 8) there is often a seperate file, shortened tot type_segment.extension.
<pre>
.
├── Aligned
│   ├── type_segment.fasta
├── Clean
│   ├── type_segment.fasta
├── Clusters
│   ├── type_segment.fasta
│   ├── type_segment.fasta.clstr
│   └── cluster_info.csv
├── Dataframes
│   ├── Complete
│   │   ├── df_type_segment.csv
│   ├── Scores
│   │   ├── df_type.csv
│   │   ├── df_type_segment.csv
│   ├── Ungapped
│   │   ├── df_type.csv
│   │   ├── df_type_segment.csv
│   ├── complete_score_df_type.csv
│   ├── df_type.csv
│   ├── score_df_type.csv
│   ├── small_score_df_type_Mutability.csv
│   ├── small_score_df_type_Shannon_entropy.csv
├── Genome
│   ├── Aligned
│   │   ├── type.fasta
│   │   └── type_segment.fasta
│   ├── CDS
│   │   ├── type_segment.fasta
│   │   └── type_CDS.fasta
│   ├── CDS_df_type.csv
│   └── Consensus
│       ├── type_segment.fasta
│       └── type_consensus.fasta
├── Metadata
│   ├── metadata.json
│   ├── metadata_clean.json
│   ├── metadata_cluster.json
│   └── metadata_complete.json
├── NCBI (not on GitHub, see Zenodo)
│   ├── influenza.fna
│   └── influenza_na.dat
├── Sorted
│   ├── type_segment.fasta
├── Structure
│   ├── Colors
│   │   ├── type.txt
│   │   └── type_segment.txt
│   ├── Consensus
│   │   ├── type.fasta
│   │   └── type_segment.fasta
│   ├── Dotbracket
│   │   ├── type.txt
│   │   ├── type_segment.ifold
│   │   └── type_segment.txt
│   └── Element
│       └── type_segment_001.element_string
└── Tree
    ├── type_segment.fasta
    └── type_segment.tree
</pre>




