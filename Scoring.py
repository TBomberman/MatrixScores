import numpy as np


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


def get_gene_set_score(column, up_ids, down_ids):
    up_len = len(up_ids)
    total_len = up_len + len(down_ids)
    column_arr = column.to_numpy()
    up_count = np.sum(column_arr[:up_len] > 0)
    down_count = np.sum(column_arr[up_len:] < 0)
    return (up_count + down_count) / total_len


def get_ks(ids, rank_dict, n):
    max_a = 0
    max_b = 0
    t = len(ids)
    for i in range(0, t):
        gene_entrez = ids[i]
        j = i + 1
        Vj = rank_dict[gene_entrez]
        val = j / t - Vj / n
        if val > max_a:
            max_a = val
        val = Vj / n - (j - 1) / t
        if val > max_b:
            max_b = val

    return max_a if max_a > max_b else -max_b


def get_connectivity_score(column, up_ids, down_ids):
    n = len(column)
    gene_z = column.to_dict()
    sorted_genes = sorted(gene_z, key=gene_z.get, reverse=True)
    rank_dict = {key: rank for rank, key in enumerate(sorted_genes, 1)}
    ksup = get_ks(up_ids, rank_dict, n)
    ksdown = get_ks(down_ids, rank_dict, n)
    return 0 if ksup * ksdown >= 0 else ksup - ksdown


def get_es(ids, N, gene_z, sorted_gene_list):
    max_dev = 0
    NR = 0
    Phit = 0  # should always be 0-1
    NH = 0

    for i in range(0, N):
        gene_entrez = sorted_gene_list[i]
        if gene_entrez in ids:  # hit
            NR = NR + abs(gene_z[gene_entrez])

    for i in range(0, N):
        gene_entrez = sorted_gene_list[i]

        if gene_entrez in ids:  # hit
            NH = NH + 1
            Phit = 0  # reset to recalc Phit

            for j in range(0, i+1):
                gene_entrez_inner = sorted_gene_list[j]
                if gene_entrez_inner in ids:  # hit
                    zi_inner = gene_z[gene_entrez_inner]
                    Phit = Phit + abs(zi_inner)/NR

        Pmiss = (i+1)/(N-NH)

        dev = Phit-Pmiss
        if abs(dev) > abs(max_dev):
            max_dev = dev

    return max_dev


def get_WTCS(column, up_ids, down_ids):
    n = len(column)
    gene_z = column.to_dict()
    sorted_genes_list = sorted(gene_z, key=gene_z.get, reverse=True)

    esup = get_es(up_ids, n, gene_z, sorted_genes_list)
    esdown = get_es(down_ids, n, gene_z, sorted_genes_list)
    print(esup, esdown)
    return 0 if esup * esdown >= 0 else (esup - esdown)/2



