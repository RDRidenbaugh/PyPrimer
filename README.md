
```
 ███████████             ███████████             ███                                    
░░███░░░░░███           ░░███░░░░░███           ░░░                                     
 ░███    ░███ █████ ████ ░███    ░███ ████████  ████  █████████████    ██████  ████████ 
 ░██████████ ░░███ ░███  ░██████████ ░░███░░███░░███ ░░███░░███░░███  ███░░███░░███░░███
 ░███░░░░░░   ░███ ░███  ░███░░░░░░   ░███ ░░░  ░███  ░███ ░███ ░███ ░███████  ░███ ░░░ 
 ░███         ░███ ░███  ░███         ░███      ░███  ░███ ░███ ░███ ░███░░░   ░███     
 █████        ░░███████  █████        █████     █████ █████░███ █████░░██████  █████    
░░░░░          ░░░░░███ ░░░░░        ░░░░░     ░░░░░ ░░░░░ ░░░ ░░░░░  ░░░░░░  ░░░░░     
               ███ ░███                                                                 
              ░░██████                                                                  
               ░░░░░░                                                                   
```

PyPrimer is bioinformatic pipeline designed to identify intron length variation between multiple species for the purpose of designing diagnostic multiplex PCR primers. Originally implemented in *Neodiprion* sawflies, the methodology can be applied to any system with a NCBI annotated reference genome.

PyPrimer was written in Python3 and requires the following input files to function:
- table_OGs_protein_names.txt - Orthologous groups inferred using [Broccoli](https://github.com/rderelle/Broccoli) (Derelle et al., 2020)
- GCF_021901455.1_iyNeoLeco1.1_feature_table.txt - NCBI feature tables in tab delimited format
- GCF_021901455.1_iyNeoLeco1.1_genomic.gff - NCBI GFF files

PyPrimer can optionally filter orthogroup candidates using genome wide FST population data for focal species. The original implementation focused on intron variation between  *N. lecontei* and *N. pinetum* using data from Glover et al. (2024).

PyPrimer can be run using file in the example_file directory. Feature tables and GFF files will need to be downloaded from the [NCBI ftp](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/901/455/GCF_021901455.1_iyNeoLeco1.1/) servers for each species. 

While the methodology is applicable across systems, the code itself has been written specifically to five Diprionidae reference genomes. Analyzing an equal number of species can be done simply through modification of file paths, but less than or greater than five will require editing control flow in the following functions: longest, OG_isoform_filter, OG_fst_filter, pairwise_intron, OG_intron, and intron_output. See annotation in PyPrimer.py for more information. 

**References**
- Derelle, Romain, Hervé Philippe, and John K. Colbourne. "Broccoli: combining phylogenetic and network analyses for orthology assignment." Molecular Biology and Evolution 37.11 (2020): 3389-3396.

- Glover, Ashleigh N., et al. "Recurrent selection shapes the genomic landscape of differentiation between a pair of host‐specialized haplodiploids that diverged with gene flow" Molecular Ecology 33.18 (2024): e17509.
