from Helpers.data_loader import load_dict_from_cvs, load_csv, get_feature_dict

gene_probe_dict = get_feature_dict('Data/Probes_L1000_metadata.csv', delimiter=',', key_index=1)
probe_gene_dict = get_feature_dict('Data/Probes_L1000_metadata.csv', delimiter=',', key_index=0)
lincs_dict = get_feature_dict('Data/GSE92742_Broad_LINCS_gene_info.txt', delimiter='\t', key_index=1)
landmark_dict = get_feature_dict('Data/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt', delimiter='\t', key_index=1)

def check_if_l1000_probes_are_in_lm_genes():
    for key in gene_probe_dict:
        if key not in landmark_dict:
            print(key)

# check_if_l1000_probes_are_in_lm_genes()


def check_if_tag_genes_are_in_lm_probes(): # only 24 of them are there
    probes = load_csv('Data/up-probes.grp.txt')
    for probe in probes:
        if probe[0] not in probe_gene_dict:
            print(probe)

# check_if_tags_are_in_lm_probes()

# def check_if_tags_genes_are_in_lm_probes():
# dict = load_dict_from_cvs('Data/Probes_full_metadata.csv')

def  check_if_tag_genes_are_in_L1000_total_genes():
    #  1023 tags and their genes can't be located in lincs data
    # 11694 can be found
    # 22269 rows in the tag set but only 12717 can be mapped to a gene
    # these genes will be listed out in this method
    gene_tag = get_feature_dict('Data/Probes_full_metadata.csv', delimiter=',', key_index=1)
    for gene in gene_tag:
        if gene not in lincs_dict:
            print(gene)

check_if_tag_genes_are_in_L1000_total_genes()

# these genes in the L1000 probes dict can't be found in the landmark genes
# ERO1L
# CHP
# KIAA0494
# KIAA1279
# B3GNT1
# CTSL1
# BRP44
# MTERFD1
# CCDC90A
# KIAA0528
# CD97
# GPR56
# WDR67
# PTPLAD1
# PHF15
# GPER