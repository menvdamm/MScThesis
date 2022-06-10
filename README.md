# Master's Thesis: Computational Exploration Towards Universal Vaccines for Influenza and Coronaviruses
Menno Van Damme

This repository contains the code, data and results obtained during my thesis. Note that not all data is available here, the raw NCBI downloads, the SARS-CoV-2 fasta files and metadata json were to big to put on GitHub. The complete data is available on: https://doi.org/10.5281/zenodo.6621827.

If you want to run the entire pipleine completely form scratch, the python scripts should be run in this order:
- Data.py
- Alignment.py
- Conservation.py
- Phylogentics.py
- Structure.py

Some of the SARS-CoV-2 steps had to be performed on the Ghent Univeristy supercomputer because of the large number of genome sequences and the length of the genome. These were the alignment, clustering and RNA secondary structure prediction.
