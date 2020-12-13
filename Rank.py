from Helpers.data_loader import get_feature_dict, load_gene_expression_data, printProgressBar, load_csv, load_dict_from_cvs
import json
from datetime import datetime
from Helpers.cmap import find_nth, tags2entrez_list
from Scoring import ScoreItem, get_gene_set_score, get_connectivity_score, get_WTCS


def get_ranked_instances(up_ids, down_ids):
    top50 = {}
    top50min = ScoreItem(0, "")
    experiments_dose_dict = get_feature_dict('Data/GSE70138_Broad_LINCS_sig_info.txt', '\t', 0)

    print("Loading gene expressions from gctx")
    # tagged_gene_ids = up_ids + down_ids
    tagged_gene_ids = None
    level_5_gctoo = load_gene_expression_data("Data/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx",
                                              tagged_gene_ids)
    length = len(level_5_gctoo.col_metadata_df.index)
    length = 1000
    start_time = datetime.now()

    counter = 0
    for col_name_obj in level_5_gctoo.col_metadata_df.itertuples():
        col_name = col_name_obj[0]
        counter = counter + 1
        if counter > length:
            break
        printProgressBar(counter, length, prefix='Load experiments progress, length: ' + str(length))
        column = level_5_gctoo.data_df[col_name]

        # score it
        # score_val = get_gene_set_score(column, up_ids, down_ids)
        # score_val = get_connectivity_score(column, up_ids, down_ids)
        score_val = get_WTCS(column, up_ids, down_ids)
        score = ScoreItem(score_val, col_name)

        if score < top50min and len(top50) >= 50:
            continue

        # get drug name
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


prob2entrez = json.load(open('Data/probe2entrez.json'))
up_gene_ids = tags2entrez_list(load_csv('Data/up-probes.grp.txt'), prob2entrez)
down_gene_ids = tags2entrez_list(load_csv('Data/down-probes.grp.txt'), prob2entrez)
get_ranked_instances(up_gene_ids, down_gene_ids)
