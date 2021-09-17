import logging
import numpy
import parse_cif_file
import os
import sys

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
    cif = parse_cif_file.get_coordinates(cif_file)
    outfilename = '{}/{}.csv'.format(out_dir, cif_file_name.split(".cif")[0])
    fo = open(outfilename, 'w')
    for model in cif.keys():
        seq = sorted(list(set([(i[0], i[1], i[2]) for i in cif[model].keys()])))
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
                phi = get_dihedral_angle(cif[model][phi_atoms[0]],
                                         cif[model][phi_atoms[1]],
                                         cif[model][phi_atoms[2]],
                                         cif[model][phi_atoms[3]])
                psi = get_dihedral_angle(cif[model][psi_atoms[0]],
                                         cif[model][psi_atoms[1]],
                                         cif[model][psi_atoms[2]],
                                         cif[model][psi_atoms[3]])
                fo.write('{},{},{},{}\n'.format(seq[r][0], seq[r][2], phi, psi))
            except KeyError:
                logging.warning('Coordinate data not found for {}/{}'.format(phi_atoms, psi_atoms))


if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    flist = [_ for _ in os.listdir(in_path) if _.endswith('.cif')]
    for fname in flist:
        print (fname)
        logging.info('Working on {}'.format(fname))
        calculate_dihedral_angles(fname, in_path, out_path)
