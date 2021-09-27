import pynmrstar
import json
import requests

def meta_stat():
    url = "http://api.bmrb.io/v2/list_entries?database=metabolomics"
    r = requests.get(url).json()
    stata={}

    for bid in r:
        url2='http://api.bmrb.io/v2/entry/{}?tag=Entry_src.Organization_full_name'.format(bid)
        r2 = requests.get(url2).json()
        k=r2[bid]['Entry_src.Organization_full_name'][0]
        print (k,bid)
        if 'bmst' in bid:
            k = 'Theory'
        if k not in stata.keys():
            stata[k]=0
        stata[k]+=1
    for k in stata:
        print (k,stata[k])

if __name__ == "__main__":
    meta_stat()