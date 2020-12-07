from Helpers.data_loader import get_feature_dict, load_gene_expression_data, printProgressBar, load_csv, load_dict_from_cvs
import json
from datetime import datetime
import numpy as np


def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start


class ScoreItem:
    def __init__(self,item_score,item_string):
        self.item_score= item_score
        self.item_string = item_string

    def __hash__(self):
        return hash((self.item_score, self.item_string))

    def __eq__(self, other):
        return (self.item_score, self.item_string) == (other.item_score, other.item_string)

    def __lt__(self, other):
        if self.item_score < other.item_score:
            return True
        return self.item_string < other.item_string


def get_ranked_instances(up_ids, down_ids):
    top50 = {}
    top50min = ScoreItem(0, "")

    tagged_gene_ids = up_ids + down_ids
    experiments_dose_dict = get_feature_dict('Data/GSE70138_Broad_LINCS_sig_info.txt', '\t', 0)

    print("Loading gene expressions from gctx")
    level_5_gctoo = load_gene_expression_data("Data/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx",
                                              tagged_gene_ids)
    length = len(level_5_gctoo.col_metadata_df.index)
    # length = 100
    total_genes = len(tagged_gene_ids)
    start_time = datetime.now()

    up_len = len(up_ids)

    counter = 0
    for col_name_obj in level_5_gctoo.col_metadata_df.itertuples():
        col_name = col_name_obj[0]
        counter = counter + 1
        printProgressBar(counter, length, prefix='Load experiments progress, length: ' + str(length))
        column = level_5_gctoo.data_df[col_name]
        column_arr = column.to_numpy()

        # score it
        up_count = np.sum(column_arr[:up_len] > 0)
        down_count = np.sum(column_arr[up_len:] < 0)
        score = ScoreItem((up_count + down_count) / total_genes, col_name)

        if score < top50min and len(top50) >= 50:
            continue

        # get drug features
        col_name_key = col_name  # [2:-1]
        if col_name_key not in experiments_dose_dict:
            continue
        experiment_data = experiments_dose_dict[col_name_key]
        drug_name = experiment_data[1]

        # parse the time
        start = col_name.rfind("_")
        end = find_nth(col_name, ":", 1)
        exposure_time = col_name[start + 1:end]

        # parse the dosage unit and value
        dose_unit = experiment_data[4][-2:]
        dose_amt = float(experiment_data[4][:-2])

        # parse the cell name
        start = find_nth(col_name, "_", 1)
        end = find_nth(col_name, "_", 2)
        cell_name = col_name[start + 1:end]

        # get the new min, and remove the old min
        if len(top50) >= 50:
            top50.pop(top50min)

        print_str = drug_name + " " + cell_name + " " + str(dose_amt) + dose_unit + " " + exposure_time
        top50[score] = print_str
        top50min = min(list(top50.keys()))
        print(datetime.now(), "{:.3f}".format(score.item_score), print_str)

    elapsed = datetime.now() - start_time

    print("RESULTS, running time: " + str(elapsed))

    sorted_keys = sorted(top50.keys(), key=lambda key: key.item_score, reverse=True)
    for key in sorted_keys:
        print("{:.3f}".format(key.item_score), top50[key])


def get_gene_id_dict():
    lm_genes = json.load(open('Data/landmark_genes.json'))
    dict = {}
    for lm_gene in lm_genes:
        # dict[lm_gene['gene_symbol']] = lm_gene['entrez_id']
        dict[lm_gene['entrez_id']] = lm_gene['gene_symbol']
    return dict


def tags2entrez_list(tags, probes):
    id_list = []
    for tag in tags:
        if tag[0] in probes:
            entrez_id = probes[tag[0]]
            if entrez_id != "":
                if entrez_id not in id_list:
                    id_list.append(entrez_id)
    return id_list


prob2entrez = json.load(open('Data/probe2entrez.json'))
up_gene_ids = tags2entrez_list(load_csv('Data/up-probes.grp.txt'), prob2entrez)
down_gene_ids = tags2entrez_list(load_csv('Data/down-probes.grp.txt'), prob2entrez)
get_ranked_instances(up_gene_ids, down_gene_ids)
