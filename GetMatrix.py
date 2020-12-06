from Helpers.data_loader import load_dict_from_cvs, load_csv

dict = load_dict_from_cvs('Data/Probes_full_metadata.csv')


def get_file_genes(filename):
    list = []
    probes = load_csv(filename)
    for probe in probes:
        if probe[0] in dict:
            if dict[probe[0]] != "":
                if dict[probe[0]] not in list:
                    list.append(dict[probe[0]])

    for gene in list:
        print(gene)
print("UP==============================")
get_file_genes('Data/up-probes.grp.txt')

print("DOWN==============================")
get_file_genes('Data/down-probes.grp.txt')
