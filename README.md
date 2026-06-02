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

PyPrimer was written in Python3.8+ and requires the following input files to function:
* Orthologous groups inferred for a set of species using [Broccoli](https://github.com/rderelle/Broccoli) (Derelle et al., 2020)
* NCBI feature tables in tab delimited format for all species within the orthologous groups.
* NCBI GFF files.

PyPrimer can optionally filter orthogroup candidates using genome wide FST population data for focal species. The original implementation focused on intron variation between  *N. lecontei* and *N. pinetum* from Glover et al. (2024) using VCFtools.

PyPrimer can be run as a tutorial using files from the example_file directory on GitHub. Feature tables and GFF files will need to be downloaded from the [NCBI ftp](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/901/455/GCF_021901455.1_iyNeoLeco1.1/) servers for each species.

## Installation

```
pip install pyprimer-ilps
```

## Documentation

```
# All genomes, no FST filter
pyprimer --genomes example_files/genomes.tsv \
         --broccoli-table example_files/broccoli_gene_orthology/run1/dir_step3/table_OGs_protein_names.txt \
         --annotation-type CDS \
         --output-dir example_files/out/focal_pair/ \
         --min-intron-diff 100 --max-intron-len 1100

# Focal pair with FST filtering
pyprimer --genomes example_files/genomes.tsv \
         --broccoli-table example_files/broccoli_gene_orthology/run1/dir_step3/table_OGs_protein_names.txt \
         --annotation-type CDS \
         --output-dir example_files/out/focal_pair/ \
         --use-fst --fst-table example_files/HighSites_FILTERED_fst_GD_GC_pi_TD_RR_dxy_50kbp.txt --fst-cutoff 0.75 \
         --focal-pair pine leco \
         --min-intron-diff 100 --max-intron-len 1100
```

**Mandatory arguments:**\
`--genomes`: path to the tab-seperate value file idetifying the number of species being considered, and the paths to their respecitive feature tables and .gff files.\
`--broccoli-table`: path to the "table_OGs_protein_names.txt" file output by broccoli in the "dir_step3" directory.\
`--annotation-type`: CDS for designing PCR primers, mRNA for RT-qPCR primers. \
`--output-dir out`: path to your desired output directory.\
`--min-intron-diff` and `--max-intron-len`: The minumum required difference between introns and maximum intron length desired between any two species.\

**Optional Arguments:**\
`--use-fst`: Arugment to enable FST filtering.\
`--fst-table` and `--fst-cutoff`: path to the file containing windowed FST values, along with the minimum desired FST threshold. Used in tandem with `--use-fst`.\
`--focal-pair`: Argument to only search between two focal species.\


## References
* Derelle, Romain, Hervé Philippe, and John K. Colbourne. "Broccoli: combining phylogenetic and network analyses for orthology assignment." Molecular Biology and Evolution 37.11 (2020): 3389-3396.

*-* Glover, Ashleigh N., et al. "Recurrent selection shapes the genomic landscape of differentiation between a pair of host‐specialized haplodiploids that diverged with gene flow" Molecular Ecology 33.18 (2024): e17509.
