import plotly.express as px
import csv


def plot_data(fname):
    with open(fname, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        x=[]
        y=[]
        res=[]
        for r in spamreader:
            x.append(float(r[2]))
            y.append(float(r[3]))
            res.append(r[1])
    fig = px.scatter(x=x,y=y,color=res,facet_col=res,facet_col_wrap=5)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.write_image('/Users/kumaran/homo_sapiens.jpg',width=1600,height=1400)
    fig.write_image('/Users/kumaran/homo_sapiens.pdf',width=1600,height=1400)

    #fig.show()

if __name__ == "__main__":
    plot_data('/Users/kumaran/homo_sapiens.csv')