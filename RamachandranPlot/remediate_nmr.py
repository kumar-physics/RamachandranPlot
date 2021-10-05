import csv
import os
import sys



def remediate_nmr(fname):
    # Filters can be introduced later here.
    csv_data = {}
    k = []
    write_file=False
    with open(fname) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        fo = open('/tmp/nmr.csv','w')
        for row in spamreader:
            if 'NMR' in row[4]:
                write_file=True
            else:
                break
            k1 = '{}-{}'.format(row[0],row[1])
            if k1 not in k:
                k.append(k1)
                fo.write('{}\n'.format(','.join(row)))
                #print (k1,','.join(row))
            else:
                break
        fo.close()
    if write_file:
        os.system('mv /tmp/nmr.csv {}'.format(fname))



if __name__=="__main__":
    in_path = sys.argv[1]
    # in_path = '/Users/kumaran/ramachandran'
    flist = [_ for _ in os.listdir(in_path) if _.endswith('.csv')]
    for fname in flist:
        print(fname)
        remediate_nmr('{}/{}'.format(in_path,fname))
    # read_csv('/Users/kumaran/ramachandran/2kpu.csv')