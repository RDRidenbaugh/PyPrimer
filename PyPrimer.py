from re import split, findall

def name(protein_names):
    with open(protein_names, "r") as pn:
        temp_dict = {}
        for line in pn.readlines():
            line = line.strip("\n")
            if line.startswith("#"):
                continue
            else:
                temp_split = line.split("\t")
                temp_value = []
                if temp_split[1] == '':
                    temp_value = [["NA"]]
                elif len(temp_split[1]) == 14:
                    temp_value = [[temp_split[1]]]
                else:
                    temp_value = [temp_split[1].split(" ")]
                if temp_split[2] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[2]) == 14:
                    temp_value.append([temp_split[2]])
                else:
                    temp_value.append(temp_split[2].split(" "))
                if temp_split[3] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[3]) == 14:
                    temp_value.append([temp_split[3]])
                else:
                    temp_value.append(temp_split[3].split(" "))
                if temp_split[4] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[4]) == 14:
                    temp_value.append([temp_split[4]])
                else:
                    temp_value.append(temp_split[4].split(" "))
                if temp_split[5] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[5]) == 14:
                    temp_value.append([temp_split[5]])
                else:
                    temp_value.append(temp_split[5].split(" "))
                temp_dict[temp_split[0]] = temp_value
    pn.close()
    return temp_dict

# This function processes the orthogroups identfied by Broccoli. 
# Modify the argument below to the path of your file location. 
OG_names = name("../broccoli_gene_orthology/run1/dir_step3/table_OGs_protein_names.txt")

def build_db_ft(path, feature):
    temp_dict = {}
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith(feature):
                if feature == "mRNA":
                    temp_split = line.strip().split("\t")
                    temp_dict[temp_split[10]] = (temp_split[14], temp_split[10], temp_split[12], int(temp_split[18]), int(temp_split[17]), temp_split[5], temp_split[13])
                if feature == "CDS":
                    temp_split = line.strip().split("\t")
                    temp_dict[temp_split[10]] = (temp_split[14], temp_split[12], temp_split[10], int(temp_split[18]), int(temp_split[17]), temp_split[5], temp_split[13], temp_split[7], temp_split[8], temp_split[9])
    f.close()
    return temp_dict

# These functions assemble feature databases using NCBI feature tables. Each species requires its own line and unique variable name.
# The object names below are used within several functions. If you wish to modify these objects (add/remove species), modifications will be required in downstream functions. 
# This function can build databases using CDS or mRNA annotations. Modify the second argument for your desired output.
# Modify the first argument below to the path of your file location.
pine_db_ft = build_db_ft("../feature_table/GCF_021155775.1_iyNeoPine1.1_feature_table.txt", "CDS")
leco_db_ft = build_db_ft("../feature_table/GCF_021901455.1_iyNeoLeco1.1_feature_table.txt", "CDS")
virg_db_ft = build_db_ft("../feature_table/GCF_021901495.1_iyNeoVirg1.1_feature_table.txt", "CDS")
fabr_db_ft = build_db_ft("../feature_table/GCF_021155785.1_iyNeoFabr1.1_feature_table.txt", "CDS")
simi_db_ft = build_db_ft("../feature_table/GCF_021155765.1_iyDipSimi1.1_feature_table.txt", "CDS")

def build_db_exon(path, feature):
    XP_dict = {}
    with open(path, "r") as f:
        for line in f.readlines():
            if findall(feature, line) and not findall("partial=true", line):
                temp_split = line.strip("\n").split("\t")
                meta_split = split(r'[;=]',temp_split[-1])
                if meta_split[-1] not in XP_dict:
                    XP_dict[meta_split[-1]] = [temp_split[6], [[temp_split[3], temp_split[4]]]]
                else:
                    if XP_dict[meta_split[-1]][0] == "+":
                        XP_dict[meta_split[-1]][1].append([temp_split[3], temp_split[4]])
                    else:
                        XP_dict[meta_split[-1]][1].insert(0, [temp_split[3], temp_split[4]])
    f.close()
    for key in XP_dict.keys():
        exon_list = XP_dict[key][-1]
        exon_len = len(exon_list)
        if exon_len == 1:
            temp_list_exon = [abs(int(exon_list[0][1])-int(exon_list[0][0]))+1]
            XP_dict[key] = temp_list_exon
        else:
            temp_list_exon = []
            temp_list_intron = []
            if XP_dict[key][0] == "+":
                for sublist in exon_list:
                    temp_len = abs(int(sublist[1])-int(sublist[0]))+1
                    temp_list_exon.append(temp_len)
                for i, sublist in enumerate(exon_list):
                    if i == 0:
                        temp_start = sublist[1]
                    elif i < exon_len-1 and i != 0:
                        temp_len = abs(int(temp_start)-int(sublist[0]))
                        temp_list_intron.append(temp_len)
                        temp_start = sublist[1]
                    else:
                        temp_len = abs(int(temp_start)-int(sublist[0]))
                        temp_list_intron.append(temp_len)
            else:
                for sublist in exon_list:
                    temp_len = abs(int(sublist[1])-int(sublist[0]))+1
                    temp_list_exon.insert(0, temp_len)
                for i, sublist in enumerate(exon_list):
                    if i == 0:
                        temp_start = sublist[1]
                    elif i < exon_len-1 and i != 0:
                        temp_len = abs(int(temp_start)-int(sublist[0]))
                        temp_list_intron.insert(0, temp_len)
                        temp_start = sublist[1]
                    else:
                        temp_len = abs(int(temp_start)-int(sublist[0]))
                        temp_list_intron.insert(0, temp_len)
            XP_dict[key] = [temp_list_exon, temp_list_intron]
    return XP_dict

# These functions assemble exon databases using NCBI GFF files. Each species requires its own line and unique variable name.
# The object names below are used within several functions. If you wish to modify these objects (add/remove species), modifications will be required in downstream functions. 
# This function can build databases using CDS or mRNA annotations. Modify the second argument for your desired output.
# Modify the first argument below to the path of your file location.
pine_db_exon = build_db_exon("../genome_gff/GCF_021155775.1_iyNeoPine1.1_genomic.gff", "CDS")
leco_db_exon = build_db_exon("../genome_gff/GCF_021901455.1_iyNeoLeco1.1_genomic.gff", "CDS")
virg_db_exon = build_db_exon("../genome_gff/GCF_021901495.1_iyNeoVirg1.1_genomic.gff", "CDS")
fabr_db_exon = build_db_exon("../genome_gff/GCF_021155785.1_iyNeoFabr1.1_genomic.gff", "CDS")
simi_db_exon = build_db_exon("../genome_gff/GCF_021155765.1_iyDipSimi1.1_genomic.gff", "CDS")

def fst_metadata(path):
    temp_dict = {}
    with open(path, "r") as f:
        for line in f.readlines():
            temp_split = line.strip("\n").replace('"', "").split("\t")
            key = temp_split[2]+"-"+temp_split[1]
            temp_dict[key] = (temp_split[1], temp_split[5], temp_split[2], temp_split[6], temp_split[4])
    return temp_dict

# This function loads the optional population data for downstream filtering of orthogroups. 
# Modify the argument below to the path of your file location.
fst_db = fst_metadata("../HighSites_FILTERED_fst_GD_GC_pi_TD_RR_dxy_50kbp.txt")

def longest(list, species):
    if species == "pine":
        high = 0
        XP = ""
        for i in range(len(list)):
            if pine_db_ft[list[i]][3] > high:
                high = pine_db_ft[list[i]][3]
                XP = list[i]
    elif species == "leco":
        high = 0
        XP = ""
        for i in range(len(list)):
            if leco_db_ft[list[i]][3] > high:
                high = leco_db_ft[list[i]][3]
                XP = list[i]
    elif species == "virg":
        high = 0
        XP = ""
        for i in range(len(list)):
            if virg_db_ft[list[i]][3] > high:
                high = virg_db_ft[list[i]][3]
                XP = list[i]
    elif species == "fabr":
        high = 0
        XP = ""
        for i in range(len(list)):
            if fabr_db_ft[list[i]][3] > high:
                high = fabr_db_ft[list[i]][3]
                XP = list[i]
    elif species == "simi":
        high = 0
        XP = ""
        for i in range(len(list)):
            if simi_db_ft[list[i]][3] > high:
                high = simi_db_ft[list[i]][3]
                XP = list[i]
    return XP

def OG_isoform_filter(dict):
    temp_dict = {}
    for key, list in dict.items():
        if len(list[0]) == 1 and len(list[1]) == 1 and len(list[2]) == 1 and len(list[3]) == 1 and len(list[4]) == 1:
            temp_dict[key] = (list[0][0], list[1][0], list[2][0], list[3][0], list[4][0])
        else:
            if len(list[0]) == 1:
                if list[0][0] == "NA":
                    temp_dict[key] = ("NA",)
                else:
                    temp_dict[key] = (list[0][0],)
            elif len(list[0]) > 1:
                temp_dict[key] = (longest(list[0], "pine"),)  
            if len(list[1]) == 1:
                if list[1][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[1][0],)
            elif len(list[1]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[1], "leco"),)
            if len(list[2]) == 1:
                if list[2][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[2][0],)
            elif len(list[2]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[2], "virg"),)
            if len(list[3]) == 1:
                if list[3][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[3][0],)
            elif len(list[3]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[3], "fabr"),)
            if len(list[4]) == 1:
                if list[4][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[4][0],)
            elif len(list[4]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[4], "simi"),)
    return temp_dict

# This function filters isoforms by species within orthogroups. Only the longest isoform is retained. 
OG_names = OG_isoform_filter(OG_names)

def overlap(gene_start, gene_stop, fst_start, fst_stop):
    if max(max((int(fst_stop)-int(gene_start)), 0) - max((int(fst_stop)-int(gene_stop)), 0) - max((int(fst_start)-int(gene_start)), 0), 0) > 0:
        return True
    else:
        return False

def OG_fst_filter(OG, fst, cutoff):
    temp_dict = {}
    for key in OG.keys():
        for region in fst.keys():
            try:
                if leco_db_ft[OG[key][1]][5] == fst[region][0]:
                    if overlap(leco_db_ft[OG[key][1]][7], leco_db_ft[OG[key][1]][8], fst[region][1], fst[region][3]) and float(fst[region][-1]) >= cutoff:
                        temp_value = OG[key]
                        temp_dict[key] = temp_value + (fst[region][-1],)
            except KeyError:
                continue
    return temp_dict

# This function performs optional filtering of OGs using poulation data.
# In the example dataset, genomewide FST between Lexington, KY Neodiprion lecontei and N. pinetum was used with a 0.75 cutoff.
# Modify the third argument to change the FST cutoff value. 
OG_names = OG_fst_filter(OG_names, fst_db, 0.75)

def pairwise_intron(key, dif, mode, n):
    try:
        if len(pine_db_exon[key[0]][1]) == len(leco_db_exon[key[1]][1]) == len(virg_db_exon[key[2]][1]) == len(fabr_db_exon[key[3]][1]) == len(simi_db_exon[key[4]][1]):
            if mode == "pl":
                for p, l in zip(pine_db_exon[key[0]][1], leco_db_exon[key[1]][1]):
                    if abs(int(p)-int(l)) > dif:
                        if p < n and l < n:
                            return True
            elif mode == "all":
                for p, l in zip(pine_db_exon[key[0]][1], leco_db_exon[key[1]][1]):
                    if abs(int(p)-int(l)) > dif:
                        if p < n and l < n:
                            return True
                for p, v in zip(pine_db_exon[key[0]][1], virg_db_exon[key[2]][1]):
                    if abs(int(p)-int(v)) > dif:
                        if p < n and v < n:
                            return True
                for p, f in zip(pine_db_exon[key[0]][1], fabr_db_exon[key[3]][1]):
                    if abs(int(p)-int(f)) > dif:
                        if p < n and f < n:
                            return True
                for p, s in zip(pine_db_exon[key[0]][1], simi_db_exon[key[4]][1]):
                    if abs(int(p)-int(s)) > dif:
                        if p < n and s < n:
                            return True
                for l, v in zip(leco_db_exon[key[1]][1], virg_db_exon[key[2]][1]):
                    if abs(int(l)-int(v)) > dif:
                        if l < n and v < n:
                            return True
                for l, f in zip(leco_db_exon[key[1]][1], fabr_db_exon[key[3]][1]):
                    if abs(int(l)-int(f)) > dif:
                        if l < n and f < n:
                            return True
                for l, s in zip(leco_db_exon[key[1]][1], simi_db_exon[key[4]][1]):
                    if abs(int(l)-int(s)) > dif:
                        if l < n and s < n:
                            return True
                for v, f in zip(virg_db_exon[key[2]][1], fabr_db_exon[key[3]][1]):
                    if abs(int(v)-int(f)) > dif:
                        if v < n and f < n:
                            return True
                for v, s in zip(virg_db_exon[key[2]][1], simi_db_exon[key[4]][1]):
                    if abs(int(v)-int(s)) > dif:
                        if v < n and s < n:
                            return True
                for f, s in zip(fabr_db_exon[key[3]][1], simi_db_exon[key[4]][1]):
                    if abs(int(f)-int(s)) > dif:
                        if f < n and s < n:
                            return True
        else:
            return False
    except: 
        return False

def OG_intron(dict, dif, mode, n):
    temp_dict = {}
    for key in dict.keys():
        if dict[key][0] == "NA" or dict[key][1] == "NA" or dict[key][2] == "NA" or dict[key][3] == "NA" or dict[key][4] == "NA": 
            continue
        elif pairwise_intron(dict[key], dif, mode, n):
            temp_dict[key] = ([dict[key][0], pine_db_exon[dict[key][0]][1]], [dict[key][1],leco_db_exon[dict[key][1]][1]], [dict[key][2], virg_db_exon[dict[key][2]][1]], [dict[key][3], fabr_db_exon[dict[key][3]][1]], [dict[key][4], simi_db_exon[dict[key][4]][1]], dict[key][-1])
        else:
            continue
    return temp_dict

# This function performs pairwise comparisons of intron length for all retained OGs. 
# Three arguments can be modified in this function to fit your desired results. The first and third argument control your minimum and maximum for intron length differences.
# The second argument controls the paireise search type, either between two focal species ("pl" for N. pinetum and N. lecontei) or all species ("all")
intron = OG_intron(OG_names, 100, "pl", 1100)
print(len(intron))

def intron_output(dict, path):
    for key in dict:
        with open(path+key+".txt", "w") as f:
            f.write("# "+key+"\t"+pine_db_ft[dict[key][0][0]][0]+"\t"+pine_db_ft[dict[key][0][0]][6]+"\t"+dict[key][-1]+"\n")
            f.write("pine"+"\t"+dict[key][0][0]+"\t")
            length = len(dict[key][0][1])
            for i, p in enumerate(dict[key][0][1]):
                if i < length-1:
                    f.write(str(p)+"\t")
                else:
                    f.write(str(p)+"\n")
            f.write("leco"+"\t"+dict[key][1][0]+"\t")
            for i, l in enumerate(dict[key][1][1]):
                if i < length-1:
                    f.write(str(l)+"\t")
                else:
                    f.write(str(l)+"\n")
            f.write("virg"+"\t"+dict[key][2][0]+"\t")
            for i, v in enumerate(dict[key][2][1]):
                if i < length-1:
                    f.write(str(v)+"\t")
                else:
                    f.write(str(v)+"\n")
            f.write("fabr"+"\t"+dict[key][3][0]+"\t")
            for i, fa in enumerate(dict[key][3][1]):
                if i < length-1:
                    f.write(str(fa)+"\t")
                else:
                    f.write(str(fa)+"\n")
            f.write("simi"+"\t"+dict[key][4][0]+"\t")
            for i, s in enumerate(dict[key][4][1]):
                if i < length-1:
                    f.write(str(s)+"\t")
                else:
                    f.write(str(s)+"\n")
        f.close()

# This function outputs the results of the pairwise search. The length of each intron, for each species, is output in a seperate text file for each retained OG.
# Modify the argument below to the path of your file location.
intron_output(intron, "../intron_out/")