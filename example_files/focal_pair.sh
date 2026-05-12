#!/bin/bash

python3 PyPrimer.py --genomes example_files/genomes.tsv \
                    --broccoli-table example_files/broccoli_gene_orthology/run1/dir_step3/table_OGs_protein_names.txt \
                    --output-dir example_files/out/focal_pair \
                    --use-fst --fst-table example_files/HighSites_FILTERED_fst_GD_GC_pi_TD_RR_dxy_50kbp.txt --fst-cutoff 0.75 \
                    --focal-pair pine leco \
                    --min-intron-diff 100 --max-intron-len 1100