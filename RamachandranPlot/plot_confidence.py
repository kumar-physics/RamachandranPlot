import plotly.express as px
import csv
import numpy



proteomes = {
    'UP000006548':'Arabidopsis',
    'UP000001940': 'Nematode worm',
    'UP000000559': 'C. albicans',
    'UP000000437': 'Zebrafish',
    'UP000002195': 'Dictyostelium',
    'UP000000803': 'Fruit fly',
    'UP000000625': 'E. coli',
    'UP000008827': 'Soybean',
    'UP000005640': 'Human',
    'UP000008153': 'L. infantum',
    'UP000000805': 'M. jannaschii',
    'UP000000589': 'Mouse',
    'UP000001584': 'M. tuberculosis',
    'UP000059680': 'Asian rice',
    'UP000001450': 'P. falciparum',
    'UP000002494': 'Rat',
    'UP000002311': 'Budding yeast',
    'UP000002485': 'Fission yeast',
    'UP000008816': 'S. aureus',
    'UP000002296': 'T. cruzi',
    'UP000007305': 'Maize'
}

def confidence_hist(inpath):
    x=[]
    y=[]
    c=[]
    grid_size = 5
    s = list(numpy.arange(0, 101, grid_size))
    ds = grid_size
    for proteome in proteomes:
        csv_file = '{}/{}.txt'.format(inpath,proteome)
        with open(csv_file) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            p=[]
            for row in spamreader:
                p.append(float(row[4]))
        x1, y1 = cal_hist(p,s,ds)
        for i in range(len(x1)):
            x.append(x1[i])
            y.append(y1[i])
            c.append(proteomes[proteome])
    fig=px.bar( x=x,y=y, color=c, barmode="overlay",opacity=0.7,labels={'x':'Confidence', 'y': 'Density'})
    fig.show()
    fig.write_html('/Users/kumaran/af_dihearal/af_all/conf_hist.heml')


def cal_hist(c, s, ds):
    x = list([0] * len(s))
    n=len(c)
    for i in c:
        x[s.index(round(float(i)/ds)*ds)]+=1
    x1=[]
    for i in x:
        x1.append(float(i)/float(n))
    print (sum(x1))
    return s,x1

if __name__ == "__main__":
    confidence_hist('/Users/kumaran/af_dihearal/af_all')