from Helpers.data_loader import load_dict_from_cvs, load_csv
# import urllib3
# http = urllib3.PoolManager()
# import json
#
url = "https://api.clue.io/api/probeset_to_entrez_id/convert"
# data = "[\"1007_s_at\",\"121_at\",\"1255_g_at\",\"1438_at\",\"1487_at\"]"
data = "[121_at]"
# data = '["1007_s_at","121_at","1255_g_at","1438_at","1487_at"]'
headers = {
    'user_key': '29e948d88629d2eea6ea26bd6518a4bb'
}
# encoded_body = json.dumps(data)
# response = http.request('POST', url, fields=encoded_body, headers=headers)
# for thing in response:
#     print(thing)
# # print(response.read())

#
# import requests
# response = requests.post(url, headers=headers, data=data)
# print(response.json())

dict = load_dict_from_cvs('Data/Probes_full_metadata.csv')
urlstr = 'curl -X POST --header "Content-Type: application/json" -r "Accept: application/json" --header "user_key: ' \
         '29e948d88629d2eea6ea26bd6518a4bb" -d "[\\"'
for probe in dict:
    if probe != 'pr_id':
        urlstr = urlstr + probe + '\\",\\"'

urlstr = urlstr[:-3] + "]\" \"https://api.clue.io/api/probeset_to_entrez_id/convert\" > probe2entrez"
print(urlstr)





