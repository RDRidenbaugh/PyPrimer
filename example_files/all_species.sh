#!/bin/bash

python3 PyPrimer.py --genomes example_files/genomes.tsv \
                    --broccoli-table example_files/broccoli_gene_orthology/run1/dir_step3/table_OGs_protein_names.txt \
                    --output-dir example_files/out/all_species \
                    --min-intron-diff 100 --max-intron-len 1100