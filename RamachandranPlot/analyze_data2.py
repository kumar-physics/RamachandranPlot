import numpy
import plotly.express as px
import csv
from scipy import stats
import sys

def analyze_data2(data_path,res):
    # aa = ['SER', 'HIS', 'GLU', 'GLY', 'LYS',
    #       'ALA', 'LEU', 'GLN', 'PRO', 'MET',
    #       'ASP', 'PHE', 'VAL', 'THR', 'ILE',
    #       'ASN', 'ARG', 'TYR', 'CYS', 'TRP']
    #aa = ['XXX']
    aa=[res]
    x=[]
    y=[]
    for res in aa:
        print ('Srarted {}'.format(res))
        fname = '{}/af/{}.csv'.format(data_path, res)
        with open(fname, newline='') as csvfile:
            print('Reading {}'.format(fname))
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                try:
                    con=float(row[4])
                    if con > 80.0:
                        phi = float(row[2])
                        psi = float(row[3])
                        x.append(phi)
                        y.append(psi)
                except ValueError:
                    pass
        xx, yy = numpy.mgrid[-180:180:100j,180:-180:100j]
        positions = numpy.vstack([xx.ravel(),yy.ravel()])
        values = numpy.vstack([x,y])
        kernel = stats.gaussian_kde(values)
        f1 = numpy.reshape(kernel(positions).T,xx.shape)
        lab1 = []
        lab2=[]
        for i in range(len(xx)):
            lab1.append(xx[i][i])
            lab2.append(xx[len(xx)-i-1][len(xx)-i-1])
        fig = px.imshow(f1.T,x=lab1,y=lab2,labels={'x':u"\u03D5",'y':u"\u03A8"})
        fig.update_yaxes(autorange=True)
        #fig.update_yaxes(autorange="reversed")
        fig.write_html('{}/plots2/{}_af.html'.format(data_path,res))
        fig = px.imshow(numpy.log(f1.T), x=lab1, y=lab2, labels={'x': u"\u03D5", 'y': u"\u03A8"})
        fig.update_yaxes(autorange=True)
        # fig.update_yaxes(autorange="reversed")
        fig.write_html('{}/plots2/{}_af_log.html'.format(data_path, res))
        fname = '{}/pdb/{}.csv'.format(data_path, res)
        x=[]
        y=[]
        with open(fname, newline='') as csvfile:
            print('Reading {}'.format(fname))
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                try:
                    resol=float(row[5])
                    if resol < 1.5:
                        phi = float(row[2])
                        psi = float(row[3])
                        x.append(phi)
                        y.append(psi)
                except ValueError:
                    pass
        xx, yy = numpy.mgrid[-180:180:100j, 180:-180:100j]
        positions = numpy.vstack([xx.ravel(), yy.ravel()])
        values = numpy.vstack([x, y])
        kernel = stats.gaussian_kde(values)
        f2 = numpy.reshape(kernel(positions).T, xx.shape)
        lab1 = []
        lab2 = []
        for i in range(len(xx)):
            lab1.append(xx[i][i])
            lab2.append(xx[len(xx) - i - 1][len(xx) - i - 1])
        fig = px.imshow(f2.T, x=lab1, y=lab2, labels={'x': u"\u03D5", 'y': u"\u03A8"})
        fig.update_yaxes(autorange=True)
        # fig.update_yaxes(autorange="reversed")
        fig.write_html('{}/plots2/{}_pdb.html'.format(data_path, res))
        fig = px.imshow(numpy.log(f2.T), x=lab1, y=lab2, labels={'x': u"\u03D5", 'y': u"\u03A8"})
        fig.update_yaxes(autorange=True)
        # fig.update_yaxes(autorange="reversed")
        fig.write_html('{}/plots2/{}_pdb_log.html'.format(data_path, res))
        f=f1-f2
        fig = px.imshow(f.T, x=lab1, y=lab2, labels={'x': u"\u03D5", 'y': u"\u03A8"})
        fig.update_yaxes(autorange=True)
        # fig.update_yaxes(autorange="reversed")
        fig.write_html('{}/plots2/{}_diff.html'.format(data_path, res))
        fig = px.imshow(numpy.log(f.T), x=lab1, y=lab2, labels={'x': u"\u03D5", 'y': u"\u03A8"})
        fig.update_yaxes(autorange=True)
        # fig.update_yaxes(autorange="reversed")
        fig.write_html('{}/plots2/{}_diff_log.html'.format(data_path, res))
        f=numpy.log(f1)-numpy.log(f2)
        fig = px.imshow(f.T, x=lab1, y=lab2, labels={'x': u"\u03D5", 'y': u"\u03A8"})
        fig.update_yaxes(autorange=True)
        # fig.update_yaxes(autorange="reversed")
        fig.write_html('{}/plots2/{}_diff_log2.html'.format(data_path, res))


def analyze_data(data_path):
    aa = ['SER','HIS','GLU','GLY','LYS',
          'ALA','LEU','GLN','PRO','MET',
          'ASP','PHE','VAL','THR','ILE',
          'ASN','ARG','TYR','CYS','TRP']
    grid = []
    for i in range(361):
        grid.append([])
        for j in range(361):
            grid[-1].append(0)
    for res in aa:
        fname = '{}/{}.csv'.format(data_path,res)
        with open(fname, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                phi=float(row[2])
                psi=float(row[3])
                grid = update_count(grid,phi,psi)
        x=[]
        y=[]
        z=[]
        grid = normalize(grid)
        for i in range(len(grid)):
            for j in range(len(grid[i])):
                x.append(i-180)
                y.append(j-180)
                z.append(grid[i][j])
        fig = px.scatter_3d(x=x,y=y,z=z)
        fig.write_html('{}/{}.html'.format(data_path,res))
        print ('Done {}'.format(res))




    x = update_count(x,-0.3)

def normalize(z):
    c=0
    for i in z:
        for j in i:
            c+=j
    for i in range(len(z)):
        for j in range(len(z[i])):
            z[i][j]=float(z[i][j])/float(c)
    return z
def update_count(z,phi,psi):
    x_index = round(phi)+180
    y_index = round(psi)+180
    #print (x_index,y_index)
    z[x_index][y_index]+=1
    return z


if __name__=="__main__":
    res=sys.argv[1]
    analyze_data2('/home/nmrbox/kbaskaran/af_dihedral',res)