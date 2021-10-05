import plotly.express as px
import sys
import csv

def plot_scatter(csv_file):
    info = []
    phi =[]
    psi = []
    with open(csv_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            info.append('{}-{}-{}'.format(row[0],row[1],row[-1]))
            phi.append(float(row[2]))
            psi.append(float(row[3]))
    fig = px.scatter(x=phi,y=psi,hover_name=info)
    fig.show()

if __name__ == "__main__":
    plot_scatter('/Users/kumaran/af_dihearal/af/af_test.csv')
