import numpy.random
import plotly.express as px
import csv
from scipy import stats
import plotly.graph_objects as go
import os

def plot_data():
    flist = [_ for _ in os.listdir('/Users/kumaran/pdb') if _.endswith('.csv')]
    x = []
    y = []
    res = []
    m = []
    rr = []
    id = []
    for fname2 in flist:
        fname = '/Users/kumaran/pdb/{}'.format(fname2)
        print (fname)
        with open(fname, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')

            for r in spamreader:
                if r[1] in ['SER','HIS','GLU','GLY','LYS',
                            'ALA','LEU','GLN','PRO','MET',
                            'ASP','PHE','VAL','THR','ILE',
                            'ASN','ARG','TYR','CYS','TRP']:
                    x.append(float(r[2]))
                    y.append(float(r[3]))
                    res.append(r[1])
                    m.append(r[4])
                    id.append(r[6])
                    try:
                        rv=float(r[5])
                        if rv<=1.0:
                            rr.append('UltraHigh')
                        elif 1.0<rv<=2.0:
                            rr.append('High')
                        elif 2.0<rv<=4.0:
                            rr.append('Medium')
                        else:
                            rr.append('Low')
                    except ValueError:
                        rr.append('NotFound')
    mm = list(set(m))
    for method in mm:
        x1=[]
        y1=[]
        r1=[]
        res1=[]
        for i in range(len(m)):
            if m[i]==mm:
                x1.append(x[i])
                y1.append(y[i])
                r1.append(rr[i])
                res1.append(res[i])
        if len(set(r1))==1:
            fig = px.scatter(x=x1, y=y1, color=res1, facet_col=res1, facet_col_wrap=5)
        else:
            fig = px.scatter(x=x1,y=y1,color=r1,facet_col=res1,facet_col_wrap=5)
        fig.update_traces(marker={'size': 1})
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_layout(
            xaxis_title=r'$\phi$',
            yaxis_title=r'$\psi}$'
        )
        # fig.write_image('/Users/kumaran/homo_sapiens.jpg',width=1600,height=1400)
        # fig.write_image('/Users/kumaran/homo_sapiens.pdf',width=1600,height=1400)
        fig.write_image('/Users/kumaran/{}.jpg'.format(m), width=1600, height=1400)
        fig.write_image('/Users/kumaran/{}.pdf'.format(m), width=1600, height=1400)
    #fig.write_html('/Users/kumaran/pdb.html')

def get_kernel(r_cutoff=1.0):
    aa=['SER','HIS','GLU','GLY','LYS',
        'ALA','LEU','GLN','PRO','MET',
        'ASP','PHE','VAL','THR','ILE',
        'ASN','ARG','TYR','CYS','TRP']
    aa_data={}
    aa_pdf ={}
    x=[]
    y=[]
    res=[]
    c=0
    flg=False
    flist = [_ for _ in os.listdir('/Users/kumaran/pdb') if _.endswith('.csv')]
    for fname2 in flist:
        fname = '/Users/kumaran/pdb/{}'.format(fname2)
        print (fname)
        with open(fname, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:

                try:
                    flg=False
                    if row[4] != 'X-RAY DIFFRACTION' or float(row[5])>r_cutoff:
                        break
                    if row[1] in aa:
                        flg=True
                        if float(row[5])<r_cutoff:
                            if row[1] not in aa_data.keys():
                                aa_data[row[1]]=[[],[]]
                            x.append(float(row[2]))
                            y.append(float(row[3]))
                            res.append(row[1])
                            aa_data[row[1]][0].append(float(row[2]))
                            aa_data[row[1]][1].append(float(row[3]))
                except ValueError:
                    pass
            if flg:
                c+=1
    print ("Number of entries = {}".format(c))
    for k in aa_data.keys():
        aa_pdf[k]=stats.gaussian_kde(numpy.vstack(aa_data[k]))
    x1=[]
    y1=[]
    z1=[]
    res1=[]
    grid = numpy.arange(-180,180,5)
    for i in grid:
        for j in grid:
            for k in aa_pdf.keys():
                p=numpy.array([numpy.array([i]),numpy.array([j])])
                x1.append(i)
                y1.append(j)
                z1.append(aa_pdf[k].evaluate(p)[0])
                res1.append(k)
    fig=px.scatter_3d(x=x1,y=y1,z=z1,color=res1)
    fig.update_traces(marker={'size': 2})
    fig.write_image('/Users/kumaran/ramachandran/pdf_15A.jpg')
    fig.write_image('/Users/kumaran/ramachandran/pdf_15A.pdf')
    fig.write_html('/Users/kumaran/ramachandran/pdf_15A.html')
    fig.update_xaxes(range=[-180, 180])
    fig.update_yaxes(range=[-180, 180])
    fig2 = px.scatter(x=x, y=y, color=res,facet_col=res,facet_col_wrap=5)
    print (len(x))
    fig2.update_traces(marker={'size': 2})
    fig2.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig2.update_xaxes(range=[-180, 180])
    fig2.update_yaxes(range=[-180, 180])
    fig2.write_image('/Users/kumaran/ramachandran/pdb_15A.jpg')
    fig2.write_image('/Users/kumaran/ramachandran/pdb_15A.pdf')
    fig2.write_html('/Users/kumaran/ramachandran/pdb_15A.html')



def test_pdf():
    def measure(n):
        m1= numpy.random.normal(size=n)
        m2=numpy.random.normal(scale=0.5, size=n)
        return m1+m2, m1-m2
    m1, m2 = measure(2000)
    xmin = m1.min()
    xmax = m1.max()
    ymin = m2.min()
    ymax = m2.max()
    X,Y = numpy.mgrid[xmin:xmax:100j,ymin:ymax:100j]
    positions = numpy.vstack([X.ravel(), Y.ravel()])
    values = numpy.vstack([m1,m2])
    kernel = stats.gaussian_kde(values)
    fig = px.scatter(x=m1,y=m2)
    gap=numpy.arange(-4.0,4.0,0.2)
    x1=[]
    x2=[]
    x3=[]
    for i in gap:
        for j in gap:
            x1.append(i)
            x2.append(j)
            p=numpy.array([numpy.array([i]),numpy.array([j])])
            x3.append(kernel.evaluate(p)[0])
    x=numpy.array(x1)
    y=numpy.array(x2)
    z=numpy.array(x3)
    # fig = px.scatter_3d(x=x1,y=x2,z=x3)
    # fig.show()

    fig = go.Figure(data=[go.Surface(x=x, y=y, z=z)])
    fig.show()
if __name__ == "__main__":
    #plot_data('/Users/kumaran/homo_sapiens.csv')
    #plot_data('/Users/kumaran/pdb.csv')
    #test_pdf()
    #get_kernel(1.5)
    plot_data()