Here we present a tool to analyse enterovirus sequencing data. It was desinged focusing on wastewater samples and Oxford Nanopore data. 

Briefly, VSEARCH tool is used to filter and cluster the generated raw reads and then perform a BLAST search against a custom reference database. A final Excel file per sample is generated containing information about the Enterovirus types found as well as its abundance. 


ev_typing_environment.yml creates a Conda environment that contains all the required packages to run the pipeline

ev_typing_nix.py is the script to analyse those samples that were amplified following the Nix et.al., 2006 protocol targeting all Enterovirus. 

ev_typing_shaw.py is the script to analyse those samples that were amplified following the Shaw et.al., 2020 protocol targeting Cluster C Enterovirus. 

ev_reference_sequences.zip contains a FASTA file with reference sequences from all Enterovirus types. There are 6303 references making it a big file. If preferred, a custom reference FASTA file can be used.
