import logging
import numpy
import parse_cif_file
import os
import sys
from operator import itemgetter

def get_dihedral_angle1(p0,p1,p2,p3):
    """http://stackoverflow.com/q/20305272/1128289"""
    p = [p0, p1, p2, p3]
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = numpy.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= numpy.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = numpy.dot(v[0], v[1])
    m = numpy.cross(v[0], b1)
    y = numpy.dot(m, v[1])
    return numpy.degrees(np.arctan2( y, x ))


def get_dihedral_angle2(p0,p1,p2,p3):
    """formula from Wikipedia article on "Dihedral angle"; formula was removed
    from the most recent version of article (no idea why, the article is a
    mess at the moment) but the formula can be found in at this permalink to
    an old version of the article:
    https://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=689165217#Angle_between_three_vectors
    uses 1 sqrt, 3 cross products"""
    # p0 = p[0]
    # p1 = p[1]
    # p2 = p[2]
    # p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b0xb1 = numpy.cross(b0, b1)
    b1xb2 = numpy.cross(b2, b1)

    b0xb1_x_b1xb2 = numpy.cross(b0xb1, b1xb2)

    y = numpy.dot(b0xb1_x_b1xb2, b1)*(1.0/numpy.linalg.norm(b1))
    x = numpy.dot(b0xb1, b1xb2)

    return numpy.degrees(numpy.arctan2(y, x))

def get_dihedral_angle(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= numpy.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - numpy.dot(b0, b1) * b1
    w = b2 - numpy.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = numpy.dot(v, w)
    y = numpy.dot(numpy.cross(b1, v), w)
    return numpy.degrees(numpy.arctan2(y, x))


def calculate_dihedral_angles(cif_file_name, in_dir, out_dir):
    cif_file = '{}/{}'.format(in_dir, cif_file_name)
    cif,bf,ent_id = parse_cif_file.get_coordinates(cif_file)
    #cif= parse_cif_file.get_coordinates(cif_file)
    outfilename = '{}/{}.csv'.format(out_dir, cif_file_name.split(".cif")[0])
    fo = open(outfilename, 'w')
    for model in cif.keys():
        seq = sorted(list(set([(i[0], i[1], i[2]) for i in cif[model].keys()])),key=itemgetter(1, 0))
        for r in range(1, len(seq) - 1):
            phi_atoms = ((seq[r - 1][0], seq[r - 1][1], seq[r - 1][2], 'C'),
                         (seq[r][0], seq[r][1], seq[r][2], 'N'),
                         (seq[r][0], seq[r][1], seq[r][2], 'CA'),
                         (seq[r][0], seq[r][1], seq[r][2], 'C'))
            psi_atoms = ((seq[r][0], seq[r][1], seq[r][2], 'N'),
                         (seq[r][0], seq[r][1], seq[r][2], 'CA'),
                         (seq[r][0], seq[r][1], seq[r][2], 'C'),
                         (seq[r + 1][0], seq[r + 1][1], seq[r + 1][2], 'N'))
            try:
                phi = get_dihedral_angle2(cif[model][phi_atoms[0]],
                                         cif[model][phi_atoms[1]],
                                         cif[model][phi_atoms[2]],
                                         cif[model][phi_atoms[3]])
                psi = get_dihedral_angle2(cif[model][psi_atoms[0]],
                                         cif[model][psi_atoms[1]],
                                         cif[model][psi_atoms[2]],
                                         cif[model][psi_atoms[3]])
                b=bf[model][phi_atoms[1]]
                if seq[r+1][2] == 'PRO':
                    rtype='XPR'
                elif r==1 or r==(len(seq)-2):
                    rtype='TER'
                elif seq[r][2] == 'GLY':
                    rtype='GLY'
                else:
                    rtype='REG'
                fo.write('{},{},{},{},{},{},{}\n'.format(seq[r][0], seq[r][2], round(phi,4), round(psi,4),b,rtype,ent_id))
            except KeyError:
                logging.warning('Coordinate data not found for {}/{}'.format(phi_atoms, psi_atoms))


if __name__ == "__main__":
    # calculate_dihedral_angles('4txr.cif','/Users/kumaran/Downloads','/Users/kumaran')
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    flist = [_ for _ in os.listdir(in_path) if _.endswith('.cif')]
    for fname in flist:
        print (fname)
        logging.info('Working on {}'.format(fname))
        calculate_dihedral_angles(fname, in_path, out_path)
