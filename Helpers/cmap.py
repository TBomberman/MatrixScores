def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

def tags2entrez_list(tags, probes):
    id_list = []
    for tag in tags:
        if tag[0] in probes:
            entrez_id = probes[tag[0]]
            if entrez_id != "":
                if entrez_id not in id_list:
                    id_list.append(entrez_id)
    return id_list