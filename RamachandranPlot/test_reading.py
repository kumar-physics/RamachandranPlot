import csv
import sys
import numpy

aa = ['SER', 'HIS', 'GLU', 'GLY', 'LYS',
      'ALA', 'LEU', 'GLN', 'PRO', 'MET',
      'ASP', 'PHE', 'VAL', 'THR', 'ILE',
      'ASN', 'ARG', 'TYR', 'CYS', 'TRP']
def test_reading(fname,n):
    if n==1:
        x=read_txt(fname)
    else:
        x=read_csv(fname)
    print (len(x))


def read_txt(fname):
    csv_data={}
    f=open(fname,'r').read().split("\n")[:-1]
    for l in f:
        row = l.split(",")
        if row[1] in aa:
            if row[1] not in csv_data.keys():
                csv_data[row[1]] = {'phi': numpy.array([]),
                                    'psi': numpy.array([])}
            csv_data[row[1]]['phi'] = numpy.append(csv_data[row[1]]['phi'], float(row[2]))
            csv_data[row[1]]['psi'] = numpy.append(csv_data[row[1]]['psi'], float(row[3]))
    return csv_data



def read_csv(fname):
    # Filters can be introduced later here.
    csv_data = {}
    with open(fname, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            if row[1] in aa:
                if row[1] not in csv_data.keys():
                    csv_data[row[1]] = {'phi': numpy.array([]),
                                        'psi': numpy.array([])}
                csv_data[row[1]]['phi'] = numpy.append(csv_data[row[1]]['phi'], float(row[2]))
                csv_data[row[1]]['psi'] = numpy.append(csv_data[row[1]]['psi'], float(row[3]))
    return csv_data
if __name__ == "__main__":
    fname = sys.argv[1]
    n = int(sys.argv[2])
    test_reading(fname,n)
