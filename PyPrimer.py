import sys
import os
from re import split, findall
from itertools import combinations


# =============================================================================
# Control file parser
# =============================================================================

def parse_control(ctrl_path):
    """Parse the PyPrimer control file and return a configuration dictionary.

    Control file format (key = value, # for comments):

        num_genomes     = 5
        genome_0        = pine | path/to/feature_table.txt | path/to/genome.gff
        genome_1        = leco | path/to/feature_table.txt | path/to/genome.gff
        ...
        broccoli_table  = path/to/table_OGs_protein_names.txt
        use_fst         = True
        fst_table       = path/to/fst_file.txt
        fst_cutoff      = 0.75
        focal_pair      = pine,leco          # only used when all_genomes = False
        all_genomes     = True
        min_intron_diff = 100
        max_intron_len  = 1100
        output_dir      = path/to/output/
        annotation_type = CDS
    """
    config = {}
    with open(ctrl_path, "r") as ctrl:
        for line in ctrl:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, _, value = line.partition("=")
            config[key.strip()] = value.strip()

    # --- required fields -------------------------------------------------------
    required = [
        "num_genomes", "broccoli_table", "use_fst",
        "all_genomes", "min_intron_diff", "max_intron_len", "output_dir",
    ]
    for field in required:
        if field not in config:
            raise ValueError(f"Control file is missing required field: '{field}'")

    # --- genome entries ---------------------------------------------------------
    num_genomes = int(config["num_genomes"])
    genomes = []
    for i in range(num_genomes):
        gkey = f"genome_{i}"
        if gkey not in config:
            raise ValueError(f"Control file is missing genome entry: '{gkey}'")
        parts = [p.strip() for p in config[gkey].split("|")]
        if len(parts) != 3:
            raise ValueError(
                f"'{gkey}' must have exactly 3 pipe-separated fields: "
                "label | feature_table | gff"
            )
        genomes.append({"label": parts[0], "ft_path": parts[1], "gff_path": parts[2]})
    config["genomes"] = genomes

    # --- type coercions --------------------------------------------------------
    config["num_genomes"]     = num_genomes
    config["use_fst"]         = config["use_fst"].strip().lower() == "true"
    config["all_genomes"]     = config["all_genomes"].strip().lower() == "true"
    config["min_intron_diff"] = int(config["min_intron_diff"])
    config["max_intron_len"]  = int(config["max_intron_len"])
    config["annotation_type"] = config.get("annotation_type", "CDS").strip()

    if config["use_fst"]:
        for field in ("fst_table", "fst_cutoff"):
            if field not in config:
                raise ValueError(
                    f"'use_fst = True' requires field '{field}' in the control file."
                )
        config["fst_cutoff"] = float(config["fst_cutoff"])

    if not config["all_genomes"]:
        if "focal_pair" not in config:
            raise ValueError(
                "'all_genomes = False' requires a 'focal_pair' field "
                "(e.g. focal_pair = pine,leco)."
            )
        config["focal_pair"] = [s.strip() for s in config["focal_pair"].split(",")]
        if len(config["focal_pair"]) != 2:
            raise ValueError("'focal_pair' must name exactly 2 genome labels.")

    return config


# =============================================================================
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


# =============================================================================
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
    """Retain only the longest isoform per species within each orthogroup.

    Parameters
    ----------
    og_dict : dict
        Orthogroup dictionary produced by name().
    ft_dbs : list of dict
        Feature-table databases in genome order matching og_dict columns.
    """
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
    """Keep orthogroups overlapping high-FST windows for the focal species.

    Parameters
    ----------
    og_dict : dict
        Filtered orthogroup dictionary.
    fst : dict
        FST metadata produced by fst_metadata().
    cutoff : float
        Minimum FST value to retain an orthogroup.
    focal_ft_db : dict
        Feature-table database for the focal (FST-reference) species.
    """
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


# =============================================================================
# Intron pairwise comparison
# =============================================================================

def pairwise_intron(og_accessions, exon_dbs, min_diff, max_len, mode, focal_indices=None):
    """Check whether any intron pair exceeds the minimum length difference.

    Parameters
    ----------
    og_accessions : tuple
        Per-species accession IDs for one orthogroup.
    exon_dbs : list of dict
        Exon databases in genome order.
    min_diff : int
        Minimum intron length difference required to flag an OG.
    max_len : int
        Maximum intron length; both introns must be shorter than this.
    mode : str
        "focal" — compare only the focal species pair;
        "all"   — compare all pairwise species combinations.
    focal_indices : tuple of int, optional
        Indices of the two focal species when mode == "focal".
    """
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
    """Filter orthogroups by pairwise intron length differences.

    Parameters
    ----------
    og_dict : dict
        Isoform-filtered (and optionally FST-filtered) orthogroup dict.
    exon_dbs : list of dict
        Exon databases in genome order.
    min_diff : int
        Minimum intron length difference.
    max_len : int
        Maximum allowed intron length.
    mode : str
        "focal" or "all" — passed directly to pairwise_intron().
    focal_indices : tuple of int, optional
        Required when mode == "focal".
    """
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


# =============================================================================
# Output
# =============================================================================

def intron_output(intron_dict, output_dir, genome_labels, ft_dbs):
    """Write one text file per retained orthogroup.

    Parameters
    ----------
    intron_dict : dict
        Result of OG_intron().
    output_dir : str
        Directory for output files (created if absent).
    genome_labels : list of str
        Species labels in genome order.
    ft_dbs : list of dict
        Feature-table databases in genome order.
    """
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


# =============================================================================
# Main
# =============================================================================

def main():
    if len(sys.argv) != 2:
        print("Usage: python PyPrimer.py <control_file>")
        sys.exit(1)

    ctrl_path = sys.argv[1]
    print(f"Reading control file: {ctrl_path}")
    config = parse_control(ctrl_path)

    genome_labels = [g["label"] for g in config["genomes"]]
    annot = config["annotation_type"]
    print(f"Genomes ({config['num_genomes']}): {', '.join(genome_labels)}")

    # Build feature-table and exon databases for every genome.
    print("Building feature-table databases...")
    ft_dbs = [build_db_ft(g["ft_path"], annot) for g in config["genomes"]]

    print("Building exon/intron databases...")
    exon_dbs = [build_db_exon(g["gff_path"], annot) for g in config["genomes"]]

    # Load Broccoli orthogroup table.
    print("Loading orthogroup names...")
    og_names = name(config["broccoli_table"])

    # Isoform filtering — keep longest per species.
    print("Filtering isoforms...")
    og_names = OG_isoform_filter(og_names, ft_dbs)

    # Optional FST filtering.
    if config["use_fst"]:
        print("Applying FST filter...")
        fst_db = fst_metadata(config["fst_table"])
        # FST filtering uses the second genome (index 1) as the focal reference,
        # matching the original script's behaviour.
        og_names = OG_fst_filter(og_names, fst_db, config["fst_cutoff"], ft_dbs[1])

    # Pairwise intron comparison.
    mode = "all" if config["all_genomes"] else "focal"
    focal_indices = None
    if mode == "focal":
        fp = config["focal_pair"]
        focal_indices = (genome_labels.index(fp[0]), genome_labels.index(fp[1]))
        print(f"Focal pair: {fp[0]} vs {fp[1]}")

    print(f"Running pairwise intron comparison (mode='{mode}')...")
    intron = OG_intron(
        og_names,
        exon_dbs,
        config["min_intron_diff"],
        config["max_intron_len"],
        mode,
        focal_indices,
    )

    # Write output.
    print(f"Writing output to: {config['output_dir']}")
    intron_output(intron, config["output_dir"], genome_labels, ft_dbs)
    print("Done.")


if __name__ == "__main__":
    main()