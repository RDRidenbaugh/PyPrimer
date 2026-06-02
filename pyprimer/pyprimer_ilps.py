import argparse
import sys
import os
from re import split, findall
from itertools import combinations


# CLI and genomes TSV parser
# =============================================================================

def parse_args():
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="PyPrimer",
        description=(
            "Identify orthogroups with diagnostic intron length differences "
            "across multiple genomes for PCR primer design."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # -- Required ---------------------------------------------------------------
    parser.add_argument(
        "--genomes",
        required=True,
        metavar="TSV",
        help=(
            "Tab-separated file listing genomes to examine. "
            "Columns (no header): label, feature_table_path, gff_path. "
            "Row order determines genome index used throughout."
        ),
    )
    parser.add_argument(
        "--broccoli-table",
        required=True,
        metavar="FILE",
        help="Broccoli orthogroup protein-name table "
             "(table_OGs_protein_names.txt).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        metavar="DIR",
        help="Directory for output files (created if absent).",
    )

    parser.add_argument(
        "--min-intron-diff",
        required=True,
        type=int,
        default=100,
        metavar="INT",
        help="Minimum intron length difference (bp) required to flag an "
             "orthogroup.",
    )
    parser.add_argument(
        "--max-intron-len",
        required=True,
        type=int,
        default=1100,
        metavar="INT",
        help="Maximum intron length (bp); both introns in a pair must be "
             "shorter than this value.",
    )

    # -- Annotation type --------------------------------------------------------
    parser.add_argument(
        "--annotation-type",
        default="CDS",
        choices=["CDS", "mRNA"],
        help="Annotation feature type to use when building databases.",
    )

    # -- FST filtering ----------------------------------------------------------
    fst_group = parser.add_argument_group("FST filtering (optional)")
    fst_group.add_argument(
        "--use-fst",
        action="store_true",
        help=(
            "Enable FST-based filtering of orthogroups. "
            "Requires --fst-table and --fst-cutoff. "
            "Uses the second genome in --genomes as the FST reference."
        ),
    )
    fst_group.add_argument(
        "--fst-table",
        metavar="FILE",
        help="Population FST metadata file.",
    )
    fst_group.add_argument(
        "--fst-cutoff",
        type=float,
        default=0.75,
        metavar="FLOAT",
        help="Minimum FST value required to retain an orthogroup.",
    )

    # -- Pairwise comparison mode -----------------------------------------------
    mode_group = parser.add_argument_group("Pairwise comparison mode (optional)")
    mode_ex = mode_group.add_mutually_exclusive_group()
    mode_ex.add_argument(
        "--all-genomes",
        action="store_true",
        default=True,
        help="Compare intron lengths across all pairwise genome combinations "
             "(default).",
    )
    mode_ex.add_argument(
        "--focal-pair",
        nargs=2,
        metavar=("LABEL_A", "LABEL_B"),
        help="Compare intron lengths only between two named genomes "
             "(e.g. --focal-pair pine leco). Overrides --all-genomes.",
    )

    args = parser.parse_args()

    # -- Cross-argument validation ----------------------------------------------
    if args.use_fst and not args.fst_table:
        parser.error("--use-fst requires --fst-table.")

    # --focal-pair overrides the --all-genomes default.
    if args.focal_pair:
        args.all_genomes = False

    return args


def parse_genomes_tsv(path):
    """Parse the genomes TSV into an ordered list of genome dicts."""
    genomes = []
    with open(path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 3:
                raise ValueError(
                    f"Genomes TSV line {line_num} has {len(cols)} column(s); "
                    "expected exactly 3 (label, feature_table, gff)."
                )
            genomes.append(
                {"label": cols[0].strip(), "ft_path": cols[1].strip(), "gff_path": cols[2].strip()}
            )
    if not genomes:
        raise ValueError(f"No genome entries found in '{path}'.")
    return genomes


# Database builders
# =============================================================================

def name(protein_names):
    """Load Broccoli orthogroup protein-name table into a dict."""
    temp_dict = {}
    with open(protein_names, "r") as pn:
        for line in pn:
            line = line.strip("\n")
            if line.startswith("#"):
                continue
            temp_split = line.split("\t")
            temp_value = []
            for col in temp_split[1:]:
                if col == '':
                    temp_value.append(["NA"])
                elif len(col) == 14:
                    temp_value.append([col])
                else:
                    temp_value.append(col.split(" "))
            temp_dict[temp_split[0]] = temp_value
    return temp_dict


def build_db_ft(path, feature):
    """Build a feature-table database from an NCBI feature table file."""
    temp_dict = {}
    with open(path, "r") as f:
        for line in f:
            if line.startswith(feature):
                temp_split = line.strip().split("\t")
                if feature == "mRNA":
                    temp_dict[temp_split[10]] = (
                        temp_split[14], temp_split[10], temp_split[12],
                        int(temp_split[18]), int(temp_split[17]),
                        temp_split[5], temp_split[13],
                    )
                elif feature == "CDS":
                    temp_dict[temp_split[10]] = (
                        temp_split[14], temp_split[12], temp_split[10],
                        int(temp_split[18]), int(temp_split[17]),
                        temp_split[5], temp_split[13],
                        temp_split[7], temp_split[8], temp_split[9],
                    )
    return temp_dict


def build_db_exon(path, feature):
    """Build an exon/intron length database from an NCBI GFF file."""
    xp_dict = {}
    with open(path, "r") as f:
        for line in f:
            if findall(feature, line) and not findall("partial=true", line):
                temp_split = line.strip("\n").split("\t")
                meta_split = split(r'[;=]', temp_split[-1])
                xp_id = meta_split[-1]
                if xp_id not in xp_dict:
                    xp_dict[xp_id] = [temp_split[6], [[temp_split[3], temp_split[4]]]]
                else:
                    if xp_dict[xp_id][0] == "+":
                        xp_dict[xp_id][1].append([temp_split[3], temp_split[4]])
                    else:
                        xp_dict[xp_id][1].insert(0, [temp_split[3], temp_split[4]])

    for key in xp_dict:
        exon_list = xp_dict[key][-1]
        exon_len = len(exon_list)
        strand = xp_dict[key][0]
        if exon_len == 1:
            xp_dict[key] = [abs(int(exon_list[0][1]) - int(exon_list[0][0])) + 1]
        else:
            temp_list_exon = []
            temp_list_intron = []
            if strand == "+":
                for sublist in exon_list:
                    temp_list_exon.append(abs(int(sublist[1]) - int(sublist[0])) + 1)
                for i, sublist in enumerate(exon_list):
                    if i == 0:
                        temp_start = sublist[1]
                    else:
                        temp_list_intron.append(abs(int(temp_start) - int(sublist[0])))
                        temp_start = sublist[1]
            else:
                for sublist in exon_list:
                    temp_list_exon.insert(0, abs(int(sublist[1]) - int(sublist[0])) + 1)
                for i, sublist in enumerate(exon_list):
                    if i == 0:
                        temp_start = sublist[1]
                    else:
                        temp_list_intron.insert(0, abs(int(temp_start) - int(sublist[0])))
                        temp_start = sublist[1]
            xp_dict[key] = [temp_list_exon, temp_list_intron]
    return xp_dict


def fst_metadata(path):
    """Load population FST metadata into a dict."""
    temp_dict = {}
    with open(path, "r") as f:
        for line in f:
            temp_split = line.strip("\n").replace('"', "").split("\t")
            key = temp_split[2] + "-" + temp_split[1]
            temp_dict[key] = (
                temp_split[1], temp_split[5], temp_split[2],
                temp_split[6], temp_split[4],
            )
    return temp_dict


# OG filtering / isoform handling
# =============================================================================

def longest(xp_list, db_ft):
    """Return the accession with the highest CDS length from db_ft."""
    high = 0
    best_xp = ""
    for xp in xp_list:
        if db_ft[xp][3] > high:
            high = db_ft[xp][3]
            best_xp = xp
    return best_xp


def OG_isoform_filter(og_dict, ft_dbs):
    """Retain only the longest isoform per species within each orthogroup."""
    temp_dict = {}
    num_sp = len(ft_dbs)
    for key, sp_lists in og_dict.items():
        resolved = []
        for i, sp_list in enumerate(sp_lists[:num_sp]):
            if len(sp_list) == 1:
                resolved.append(sp_list[0])
            else:
                resolved.append(longest(sp_list, ft_dbs[i]))
        temp_dict[key] = tuple(resolved)
    return temp_dict


def overlap(gene_start, gene_stop, fst_start, fst_stop):
    """Return True if a gene and an FST window overlap."""
    return (
        max(
            max((int(fst_stop) - int(gene_start)), 0)
            - max((int(fst_stop) - int(gene_stop)), 0)
            - max((int(fst_start) - int(gene_start)), 0),
            0,
        )
        > 0
    )


def OG_fst_filter(og_dict, fst, cutoff, focal_ft_db):
    """Keep orthogroups overlapping high-FST windows for the focal species."""
    temp_dict = {}
    for key, og in og_dict.items():
        for region, fst_vals in fst.items():
            try:
                if focal_ft_db[og[1]][5] == fst_vals[0]:
                    if (
                        overlap(
                            focal_ft_db[og[1]][7],
                            focal_ft_db[og[1]][8],
                            fst_vals[1],
                            fst_vals[3],
                        )
                        and float(fst_vals[-1]) >= cutoff
                    ):
                        temp_dict[key] = og + (fst_vals[-1],)
            except KeyError:
                continue
    return temp_dict


# Intron pairwise comparison
# =============================================================================

def pairwise_intron(og_accessions, exon_dbs, min_diff, max_len, mode, focal_indices=None):
    """Check whether any intron pair exceeds the minimum length difference."""
    try:
        intron_lists = [exon_dbs[i][og_accessions[i]][1] for i in range(len(exon_dbs))]
        num_introns = len(intron_lists[0])
        if not all(len(il) == num_introns for il in intron_lists):
            return False

        if mode == "focal":
            idx_a, idx_b = focal_indices
            for a, b in zip(intron_lists[idx_a], intron_lists[idx_b]):
                if abs(int(a) - int(b)) > min_diff and int(a) < max_len and int(b) < max_len:
                    return True

        elif mode == "all":
            for idx_a, idx_b in combinations(range(len(exon_dbs)), 2):
                for a, b in zip(intron_lists[idx_a], intron_lists[idx_b]):
                    if (
                        abs(int(a) - int(b)) > min_diff
                        and int(a) < max_len
                        and int(b) < max_len
                    ):
                        return True

        return False
    except (KeyError, IndexError):
        return False


def OG_intron(og_dict, exon_dbs, min_diff, max_len, mode, focal_indices=None):
    """Filter orthogroups by pairwise intron length differences."""
    temp_dict = {}
    num_sp = len(exon_dbs)
    for key, og in og_dict.items():
        if any(og[i] == "NA" for i in range(num_sp)):
            continue
        if pairwise_intron(og, exon_dbs, min_diff, max_len, mode, focal_indices):
            intron_entry = tuple(
                [og[i], exon_dbs[i][og[i]][1]] for i in range(num_sp)
            )
            fst_val = og[num_sp] if len(og) > num_sp else "NA"
            temp_dict[key] = intron_entry + (fst_val,)
    return temp_dict


# Output
# =============================================================================

def intron_output(intron_dict, output_dir, genome_labels, ft_dbs):
    """Write one text file per retained orthogroup."""
    os.makedirs(output_dir, exist_ok=True)
    ref_db = ft_dbs[0]

    for key, og in intron_dict.items():
        out_path = os.path.join(output_dir, key + ".txt")
        with open(out_path, "w") as f:
            ref_accession = og[0][0]
            fst_val = og[-1]
            f.write(
                "# " + key + "\t"
                + ref_db[ref_accession][0] + "\t"
                + ref_db[ref_accession][6] + "\t"
                + str(fst_val) + "\n"
            )
            num_sp = len(genome_labels)
            length = len(og[0][1])
            for i in range(num_sp):
                label = genome_labels[i]
                accession = og[i][0]
                introns = og[i][1]
                f.write(label + "\t" + accession + "\t")
                for j, val in enumerate(introns):
                    if j < length - 1:
                        f.write(str(val) + "\t")
                    else:
                        f.write(str(val) + "\n")



# Main
# =============================================================================

def main():
    args = parse_args()

    # Load genome table.
    genomes = parse_genomes_tsv(args.genomes)
    genome_labels = [g["label"] for g in genomes]
    print(f"Genomes ({len(genomes)}): {', '.join(genome_labels)}")

    # Validate focal-pair labels against the genome table.
    if args.focal_pair:
        for label in args.focal_pair:
            if label not in genome_labels:
                sys.exit(
                    f"Error: focal-pair label '{label}' not found in genomes TSV. "
                    f"Available labels: {', '.join(genome_labels)}"
                )

    # Build feature-table and exon databases for every genome.
    print("Building feature-table databases...")
    ft_dbs = [build_db_ft(g["ft_path"], args.annotation_type) for g in genomes]

    print("Building exon/intron databases...")
    exon_dbs = [build_db_exon(g["gff_path"], args.annotation_type) for g in genomes]

    # Load Broccoli orthogroup table.
    print("Loading orthogroup names...")
    og_names = name(args.broccoli_table)

    # Isoform filtering -- keep longest isoform per species.
    print("Filtering isoforms...")
    og_names = OG_isoform_filter(og_names, ft_dbs)

    # Optional FST filtering.
    if args.use_fst:
        print("Applying FST filter...")
        fst_db = fst_metadata(args.fst_table)
        # Uses the second genome (index 1) as the FST reference species,
        # matching the original script's behaviour.
        og_names = OG_fst_filter(og_names, fst_db, args.fst_cutoff, ft_dbs[1])

    # Pairwise intron comparison.
    mode = "focal" if args.focal_pair else "all"
    focal_indices = None
    if mode == "focal":
        focal_indices = (
            genome_labels.index(args.focal_pair[0]),
            genome_labels.index(args.focal_pair[1]),
        )
        print(f"Focal pair: {args.focal_pair[0]} vs {args.focal_pair[1]}")

    print(f"Running pairwise intron comparison (mode='{mode}')...")
    intron = OG_intron(
        og_names,
        exon_dbs,
        args.min_intron_diff,
        args.max_intron_len,
        mode,
        focal_indices,
    )

    # Write output.
    print(f"Writing output to: {args.output_dir}")
    intron_output(intron, args.output_dir, genome_labels, ft_dbs)
    print("Done.")


if __name__ == "__main__":
    main()