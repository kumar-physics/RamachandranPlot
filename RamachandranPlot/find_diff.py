import csv
import numpy
import plotly.express as px

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

def plot_data():
    out_dir='/Users/kumaran/af_dihearal/test'
    grid_size = 1
    s = list(numpy.arange(-180.0, 181.0, grid_size))
    ds = grid_size / 2.0
    lab_x = s
    lab_y = s
    pm={}
    for p in proteomes:
        out_file= '{}/{}'.format(out_dir,proteomes[p])
        fname='/Users/kumaran/af_dihearal/output_af/{}/ALL/ALL_af.csv'.format(p)
        m=read_csv(fname)
        pm[p]=m
        #plot_density(numpy.log(m),lab_x,lab_y,out_file,proteomes[p])

    dm=[]
    lab=[]
    all_out='{}/diff'.format(out_dir)
    for p in proteomes:
        lab.append(proteomes[p])
        dm.append([])
        for q in proteomes:
            out_file = '{}/{}_{}_2'.format(out_dir, proteomes[p],proteomes[q])
            df = pm[p]-pm[q]
            #ddf= 1+df
            plot_density(numpy.log(df),lab_x,lab_y,out_file,'{}-{}'.format(proteomes[p],proteomes[q]))
            dm[-1].append(sum(sum(abs(df))))
    plot_density2(dm,lab,lab,all_out,'Differnce')

def plot_density(df, lab_x, lab_y, out_file,tit):
    fig = px.imshow(df, x=lab_x, y=lab_y,labels={'x': u"\u03D5", 'y': u"\u03A8"})
    fig.update_yaxes(autorange=True)
    fig.update_layout(title_text=tit, title_x=0.5)
    fig.write_html('{}.html'.format(out_file))
    fig.write_image('{}.jpeg'.format(out_file), width=800, height=800)
    fig.write_image('{}.pdf'.format(out_file), width=800, height=800)

def plot_density2(df, lab_x, lab_y, out_file,tit):
    fig = px.imshow(df,x=lab_x,y=lab_y)
    fig.update_yaxes(autorange=True)
    fig.update_layout(title_text=tit, title_x=0.5)
    fig.write_html('{}.html'.format(out_file))
    fig.write_image('{}.jpeg'.format(out_file), width=800, height=800)
    fig.write_image('{}.pdf'.format(out_file), width=800, height=800)


def read_csv(csv_file):
    m=[]
    xi=765765675
    with open(csv_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            x=row[0]
            y=row[1]
            z=float(row[2])
            if x!=xi:
                m.append([])
                xi=x
            m[-1].append(z)
    return numpy.array(m)


if __name__ == "__main__":
    plot_data()