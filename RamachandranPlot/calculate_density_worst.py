import numpy
import plotly.express as px
import csv
from scipy import stats
import sys
import os
import logging

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

logging.getLogger().setLevel(logging.INFO)

aa = ['SER', 'HIS', 'GLU', 'GLY', 'LYS',
      'ALA', 'LEU', 'GLN', 'PRO', 'MET',
      'ASP', 'PHE', 'VAL', 'THR', 'ILE',
      'ASN', 'ARG', 'TYR', 'CYS', 'TRP']
amino_acid_dict = {_: True for _ in aa}


def read_csv(fname,ds):
    # Filters can be introduced later here.
    csv_data = {}
    with open(fname) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            if row[1] in amino_acid_dict:
                if row[1] not in csv_data:
                    csv_data[row[1]] = {'phi': [], 'psi': []}
                if ds=='pdb':
                    try:
                        if float(row[5])>3.0:
                            csv_data[row[1]]['phi'].append(float(row[2]))
                            csv_data[row[1]]['psi'].append(float(row[3]))
                    except ValueError:
                        pass
                elif ds=='af':
                    if float(row[4])<40:
                        csv_data[row[1]]['phi'].append(float(row[2]))
                        csv_data[row[1]]['psi'].append(float(row[3]))
        for amino_acid in csv_data:
            csv_data[amino_acid]['phi'] = numpy.array(csv_data[amino_acid]['phi'])
            csv_data[amino_acid]['psi'] = numpy.array(csv_data[amino_acid]['psi'])
    return csv_data

def cal_hist(phi,psi,s,ds):
    n=len(s)
    hist = numpy.array([[0]*n]*n)
    for i in range(len(phi)):
        xid=get_index(phi[i],s,ds)
        yid=get_index(psi[i],s,ds)
        hist[xid][yid]+=1
    hist = hist/sum(sum(hist))
    return hist


def get_index(angle,glist,ds):
    for i in glist:
        if i-ds<angle<=i+ds:
            c = glist.index(i)
            break
    return c

def calculate_density2(af_file, pdb_file, out_dir):
    logging.info("Reading af data..")
    af_data = read_csv(af_file,'af')
    logging.info("Reading pdb data..")
    pdb_data = read_csv(pdb_file,'pdb')
    logging.info('Done Reading')
    grid_size = 1
    s = list(numpy.arange(-180.0, 181.0, grid_size))
    ds = grid_size / 2.0
    lab_x = s
    lab_y = s
    for res in aa:
        out_path = '{}/{}'.format(out_dir, res)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_base_file_name = '{}/{}'.format(out_path, res)
        af_density = cal_hist(af_data[res]['phi'], af_data[res]['psi'], s, ds).T
        logging.info('Plotting {}'.format(res))
        plot_density(af_density, lab_x, lab_y, '{}_af'.format(out_base_file_name))
        write_density(af_density, lab_x, lab_y, '{}_af.csv'.format(out_base_file_name))

        log_af_density = numpy.log(af_density)
        plot_density(log_af_density, lab_x, lab_y, '{}_log_af'.format(out_base_file_name))
        write_density(log_af_density, lab_x, lab_y, '{}_log_af.csv'.format(out_base_file_name))

        pdb_density = cal_hist(pdb_data[res]['phi'], pdb_data[res]['psi'], s, ds).T
        plot_density(pdb_density, lab_x, lab_y, '{}_pdb'.format(out_base_file_name))
        write_density(pdb_density, lab_x, lab_y, '{}_pdb.csv'.format(out_base_file_name))

        log_pdb_density = numpy.log(pdb_density)
        plot_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb'.format(out_base_file_name))
        write_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb.csv'.format(out_base_file_name))

        af_pdb = af_density - pdb_density
        plot_density(af_pdb, lab_x, lab_y, '{}_af_pdb'.format(out_base_file_name))
        write_density(af_pdb, lab_x, lab_y, '{}_af_pdb.csv'.format(out_base_file_name))

        log_af_pdb = numpy.log(af_pdb)
        plot_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb'.format(out_base_file_name))
        write_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb.csv'.format(out_base_file_name))

        log_af_log_pdb = log_af_density - log_pdb_density
        plot_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb'.format(out_base_file_name))
        write_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb.csv'.format(out_base_file_name))

        pdb_af = pdb_density - af_density
        plot_density(pdb_af, lab_x, lab_y, '{}_pdb_af'.format(out_base_file_name))
        write_density(pdb_af, lab_x, lab_y, '{}_pdb_af.csv'.format(out_base_file_name))

        log_pdb_af = numpy.log(pdb_af)
        plot_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af'.format(out_base_file_name))
        write_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af.csv'.format(out_base_file_name))

        log_pdb_log_af = log_pdb_density - log_af_density
        plot_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af'.format(out_base_file_name))
        write_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af.csv'.format(out_base_file_name))

        af_data_points = len(af_data[res]['psi'])
        pdb_data_points = len(pdb_data[res]['psi'])

        fo = open('{}_info.txt'.format(out_base_file_name), 'w')
        fo.write('Number alpha fold data points : {}\n'.format(af_data_points))
        fo.write('Number pdb data points : {}\n'.format(pdb_data_points))
        fo.close()

    # density comparison for all the residues except GLY.
    af_all = {'phi': numpy.array([]), 'psi': numpy.array([]), 'q': numpy.array([])}
    pdb_all = {'phi': numpy.array([]), 'psi': numpy.array([]), 'q': numpy.array([])}
    for res in aa:
        if res != 'GLY':  # Its better to exclude glycine
            af_all['phi'] = numpy.append(af_all['phi'], af_data[res]['phi'])
            af_all['psi'] = numpy.append(af_all['psi'], af_data[res]['psi'])
            pdb_all['phi'] = numpy.append(pdb_all['phi'], pdb_data[res]['phi'])
            pdb_all['psi'] = numpy.append(pdb_all['psi'], pdb_data[res]['psi'])

    out_path = '{}/ALL'.format(out_dir)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    out_base_file_name = '{}/ALL'.format(out_path)

    af_density=cal_hist(af_all['phi'], af_all['psi'], s, ds).T
    plot_density(af_density, lab_x, lab_y, '{}_af'.format(out_base_file_name))
    write_density(af_density, lab_x, lab_y, '{}_af.csv'.format(out_base_file_name))

    log_af_density = numpy.log(af_density)
    plot_density(log_af_density, lab_x, lab_y, '{}_log_af'.format(out_base_file_name))
    write_density(log_af_density, lab_x, lab_y, '{}_log_af.csv'.format(out_base_file_name))


    pdb_density=cal_hist(pdb_all['phi'], pdb_all['psi'], s, ds).T
    plot_density(pdb_density, lab_x, lab_y, '{}_pdb'.format(out_base_file_name))
    write_density(pdb_density, lab_x, lab_y, '{}_pdb.csv'.format(out_base_file_name))

    log_pdb_density = numpy.log(pdb_density)
    plot_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb'.format(out_base_file_name))
    write_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb.csv'.format(out_base_file_name))

    af_pdb = af_density - pdb_density
    plot_density(af_pdb, lab_x, lab_y, '{}_af_pdb'.format(out_base_file_name))
    write_density(af_pdb, lab_x, lab_y, '{}_af_pdb.csv'.format(out_base_file_name))

    log_af_pdb = numpy.log(af_pdb)
    plot_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb'.format(out_base_file_name))
    write_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb.csv'.format(out_base_file_name))

    log_af_log_pdb = log_af_density - log_pdb_density
    plot_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb'.format(out_base_file_name))
    write_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb.csv'.format(out_base_file_name))

    pdb_af = pdb_density - af_density
    plot_density(pdb_af, lab_x, lab_y, '{}_pdb_af'.format(out_base_file_name))
    write_density(pdb_af, lab_x, lab_y, '{}_pdb_af.csv'.format(out_base_file_name))

    log_pdb_af = numpy.log(pdb_af)
    plot_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af'.format(out_base_file_name))
    write_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af.csv'.format(out_base_file_name))

    log_pdb_log_af = log_pdb_density - log_af_density
    plot_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af'.format(out_base_file_name))
    write_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af.csv'.format(out_base_file_name))

    af_data_points = len(af_all['psi'])
    pdb_data_points = len(pdb_all['psi'])

    fo = open('{}_info.txt'.format(out_base_file_name), 'w')
    fo.write('Number alpha fold data points : {}\n'.format(af_data_points))
    fo.write('Number pdb data points : {}\n'.format(pdb_data_points))
    fo.close()




def calculate_density(af_file, pdb_file, out_dir):
    logging.info("Reading af data..")
    af_data = read_csv(af_file)
    logging.info("Reading pdb data..")
    pdb_data = read_csv(pdb_file)
    logging.info('Done Reading')
    xx, yy = numpy.mgrid[-180:180:361j, 180:-180:361j]  # grid size
    positions = numpy.vstack([xx.ravel(), yy.ravel()])
    lab_x = []
    lab_y = []
    for i in range(len(xx)):
        lab_x.append(xx[i][i])
        lab_y.append(xx[len(xx) - i - 1][len(xx) - i - 1])
    # density comparison for each residue
    for res in aa:
        out_path = '{}/{}'.format(out_dir, res)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_base_file_name = '{}/{}'.format(out_path, res)
        values = numpy.vstack([af_data[res]['phi'], af_data[res]['psi']])
        logging.info('Calculating KDE for {}'.format(res))
        kernel = stats.gaussian_kde(values)
        logging.info('Done Calculating KDE for {}'.format(res))
        af_density = numpy.reshape(kernel(positions).T, xx.shape).T
        logging.info('Plotting {}'.format(res))
        plot_density(af_density, lab_x, lab_y, '{}_af'.format(out_base_file_name))
        write_density(af_density, lab_x, lab_y, '{}_af.csv'.format(out_base_file_name))

        log_af_density = numpy.log(af_density)
        plot_density(log_af_density, lab_x, lab_y, '{}_log_af'.format(out_base_file_name))
        write_density(log_af_density, lab_x, lab_y, '{}_log_af.csv'.format(out_base_file_name))

        values = numpy.vstack([pdb_data[res]['phi'], pdb_data[res]['psi']])
        kernel = stats.gaussian_kde(values)
        pdb_density = numpy.reshape(kernel(positions).T, xx.shape).T
        plot_density(pdb_density, lab_x, lab_y, '{}_pdb'.format(out_base_file_name))
        write_density(pdb_density, lab_x, lab_y, '{}_pdb.csv'.format(out_base_file_name))

        log_pdb_density = numpy.log(pdb_density)
        plot_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb'.format(out_base_file_name))
        write_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb.csv'.format(out_base_file_name))

        af_pdb = af_density - pdb_density
        plot_density(af_pdb, lab_x, lab_y, '{}_af_pdb'.format(out_base_file_name))
        write_density(af_pdb, lab_x, lab_y, '{}_af_pdb.csv'.format(out_base_file_name))

        log_af_pdb = numpy.log(af_pdb)
        plot_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb'.format(out_base_file_name))
        write_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb.csv'.format(out_base_file_name))

        log_af_log_pdb = log_af_density - log_pdb_density
        plot_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb'.format(out_base_file_name))
        write_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb.csv'.format(out_base_file_name))

        pdb_af = pdb_density - af_density
        plot_density(pdb_af, lab_x, lab_y, '{}_pdb_af'.format(out_base_file_name))
        write_density(pdb_af, lab_x, lab_y, '{}_pdb_af.csv'.format(out_base_file_name))

        log_pdb_af = numpy.log(pdb_af)
        plot_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af'.format(out_base_file_name))
        write_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af.csv'.format(out_base_file_name))

        log_pdb_log_af = log_pdb_density - log_af_density
        plot_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af'.format(out_base_file_name))
        write_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af.csv'.format(out_base_file_name))

        af_data_points = len(af_data[res]['psi'])
        pdb_data_points = len(pdb_data[res]['psi'])

        fo = open('{}_info.txt'.format(out_base_file_name), 'w')
        fo.write('Number alpha fold data points : {}\n'.format(af_data_points))
        fo.write('Number pdb data points : {}\n'.format(pdb_data_points))
        fo.close()
    # density comparison for all the residues except GLY.
    af_all = {'phi': numpy.array([]), 'psi': numpy.array([]), 'q': numpy.array([])}
    pdb_all = {'phi': numpy.array([]), 'psi': numpy.array([]), 'q': numpy.array([])}
    for res in aa:
        if res != 'GLY':  # Its better to exclude glycine
            af_all['phi'] = numpy.append(af_all['phi'], af_data[res]['phi'])
            af_all['psi'] = numpy.append(af_all['psi'], af_data[res]['psi'])
            pdb_all['phi'] = numpy.append(pdb_all['phi'], pdb_data[res]['phi'])
            pdb_all['psi'] = numpy.append(pdb_all['psi'], pdb_data[res]['psi'])

    out_path = '{}/ALL'.format(out_dir)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    out_base_file_name = '{}/ALL'.format(out_path)
    values = numpy.vstack([af_all['phi'], af_all['psi']])
    kernel = stats.gaussian_kde(values)

    af_density = numpy.reshape(kernel(positions).T, xx.shape).T
    plot_density(af_density, lab_x, lab_y, '{}_af'.format(out_base_file_name))
    write_density(af_density, lab_x, lab_y, '{}_af.csv'.format(out_base_file_name))

    log_af_density = numpy.log(af_density)
    plot_density(log_af_density, lab_x, lab_y, '{}_log_af'.format(out_base_file_name))
    write_density(log_af_density, lab_x, lab_y, '{}_log_af.csv'.format(out_base_file_name))

    values = numpy.vstack([pdb_all['phi'], pdb_all['psi']])
    kernel = stats.gaussian_kde(values)
    pdb_density = numpy.reshape(kernel(positions).T, xx.shape).T
    plot_density(pdb_density, lab_x, lab_y, '{}_pdb'.format(out_base_file_name))
    write_density(pdb_density, lab_x, lab_y, '{}_pdb.csv'.format(out_base_file_name))

    log_pdb_density = numpy.log(pdb_density)
    plot_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb'.format(out_base_file_name))
    write_density(log_pdb_density, lab_x, lab_y, '{}_log_pdb.csv'.format(out_base_file_name))

    af_pdb = af_density - pdb_density
    plot_density(af_pdb, lab_x, lab_y, '{}_af_pdb'.format(out_base_file_name))
    write_density(af_pdb, lab_x, lab_y, '{}_af_pdb.csv'.format(out_base_file_name))

    log_af_pdb = numpy.log(af_pdb)
    plot_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb'.format(out_base_file_name))
    write_density(log_af_pdb, lab_x, lab_y, '{}_log_af_pdb.csv'.format(out_base_file_name))

    log_af_log_pdb = log_af_density - log_pdb_density
    plot_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb'.format(out_base_file_name))
    write_density(log_af_log_pdb, lab_x, lab_y, '{}_log_af_log_pdb.csv'.format(out_base_file_name))

    pdb_af = pdb_density - af_density
    plot_density(pdb_af, lab_x, lab_y, '{}_pdb_af'.format(out_base_file_name))
    write_density(pdb_af, lab_x, lab_y, '{}_pdb_af.csv'.format(out_base_file_name))

    log_pdb_af = numpy.log(pdb_af)
    plot_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af'.format(out_base_file_name))
    write_density(log_pdb_af, lab_x, lab_y, '{}_log_pdb_af.csv'.format(out_base_file_name))

    log_pdb_log_af = log_pdb_density - log_af_density
    plot_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af'.format(out_base_file_name))
    write_density(log_pdb_log_af, lab_x, lab_y, '{}_log_pdb_log_af.csv'.format(out_base_file_name))

    af_data_points = len(af_all['psi'])
    pdb_data_points = len(pdb_all['psi'])

    fo = open('{}_info.txt'.format(out_base_file_name), 'w')
    fo.write('Number alpha fold data points : {}\n'.format(af_data_points))
    fo.write('Number pdb data points : {}\n'.format(pdb_data_points))
    fo.close()


def plot_density(df, lab_x, lab_y, out_file):
    fig = px.imshow(df, x=lab_x, y=lab_y, labels={'x': u"\u03D5", 'y': u"\u03A8"})
    fig.update_yaxes(autorange=True)
    fig.write_html('{}.html'.format(out_file))
    fig.write_image('{}.jpeg'.format(out_file), width=800, height=800)
    fig.write_image('{}.pdf'.format(out_file), width=800, height=800)


def write_density(df, lab_x, lab_y, out_file):
    fo = open(out_file, 'w')
    for i in range(len(lab_x)):
        for j in range(len(lab_y)):
            fo.write('{},{},{}\n'.format(round(lab_x[i], 2), round(lab_y[j], 2), format(df[i][j], ".4E")))
    fo.close()


#
# def read_csv(fname):
#     # Filters can be introduced later here.
#     csv_data = {}
#     with open(fname, newline='') as csvfile:
#         spamreader = csv.reader(csvfile, delimiter=',')
#         for row in spamreader:
#             if row[1] in aa:
#                 if row[1] not in csv_data.keys():
#                     csv_data[row[1]] = {'phi': numpy.array([]),
#                                         'psi': numpy.array([])}
#                 csv_data[row[1]]['phi'] = numpy.append(csv_data[row[1]]['phi'], float(row[2]))
#                 csv_data[row[1]]['psi'] = numpy.append(csv_data[row[1]]['psi'], float(row[3]))
#     return csv_data


if __name__ == "__main__":
    af_file = sys.argv[1]
    pdb_file = sys.argv[2]
    out_dir = sys.argv[3]
    calculate_density2(af_file, pdb_file, out_dir)
    # calculate_density2('/Users/kumaran/af_dihearal/af/af_test.csv', '/Users/kumaran/af_dihearal/pdb/pdb_test.csv',
    #                   '/Users/kumaran/af_dihearal/output3')
