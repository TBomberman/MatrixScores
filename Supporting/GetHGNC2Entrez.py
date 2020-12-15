import urllib3
http = urllib3.PoolManager()
from Helpers.data_loader import load_csv
import json
import requests

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

url = "https://api.clue.io/api/genes"
# "https://api.clue.io/api/genes?filter={%22where%22:{%22gene_symbol%22:%22TP53%22}}&user_key=7bc5e4093b38cae1deb3c3e0d160b5f2"

headers = {
    'user_key': '7bc5e4093b38cae1deb3c3e0d160b5f2'
}

up_hgnc = load_csv("../Data/MEK_inhibitor_up.txt")
down_hgnc = load_csv("../Data/MEK_inhibitor_down.txt")


def print_list(hgnc_genes):
    for hgnc in hgnc_genes:
        hgnc_str = hgnc[0]
        parameters = "filter={\"where\":{\"gene_symbol\":\"" + hgnc_str + "\"}}"
        response = requests.get(url + "?" + parameters, headers=headers, verify=False)
        obj = json.loads(response.text)
        if obj:
            print(obj[0]["entrez_id"])
            # print(hgnc_str, obj[0]["entrez_id"])
        # else:
            # print(hgnc_str)

print("UP ==========================")
print_list(up_hgnc)
print("DOWN ==========================")
print_list(down_hgnc)