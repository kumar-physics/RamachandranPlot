
from mmcif.io.PdbxReader import PdbxReader
import numpy
import gzip

def get_coordinates(cif_file, use_auth_tag=False):
    """
    Extract coordinate information from cif file as a dictionary
    {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
    :param cif_file: Input coordinate file
    :return: dictionary
    """
    cif_data = []
    if cif_file.endswith('.gz'):
        ifh = gzip.open(cif_file,'rt')
    else:
        ifh = open(cif_file, 'r')
    pRd = PdbxReader(ifh)
    pRd.read(cif_data)
    ifh.close()
    c0 = cif_data[0]
    ent = c0.getObj('entry')
    ent_id = ent.getValue('id')
    atom_site = c0.getObj('atom_site')
    max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
    col_names = atom_site.getAttributeList()
    model_id = col_names.index('pdbx_PDB_model_num')
    x_id = col_names.index('Cartn_x')
    y_id = col_names.index('Cartn_y')
    z_id = col_names.index('Cartn_z')
    atom_id = col_names.index('label_atom_id')
    comp_id = col_names.index('label_comp_id')
    asym_id = col_names.index('label_asym_id')
    #entity_id = col_names.index('label_entity_id')
    seq_id = col_names.index('label_seq_id')
    #icode_id = col_names.index('pdbx_PDB_ins_code')
    #alt_id = col_names.index('label_alt_id')
    aut_seq_id = col_names.index('auth_seq_id')
    aut_asym_id = col_names.index('auth_asym_id')
    aut_atom_id = col_names.index('auth_atom_id')
    aut_comp_id = col_names.index('auth_comp_id')
    bf_id=col_names.index('B_iso_or_equiv')
    pdb_models = {}
    bfactor={}
    for model in range(1, max_models + 1):
        pdb = {}
        pdb1={}
        for dat in atom_site.getRowList():
            if int(dat[model_id]) == model:
                try:
                    if use_auth_tag:
                        pdb[(int(dat[aut_seq_id]), dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                        pdb1[(int(dat[aut_seq_id]), dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                            float(dat[bf_id])
                    else:
                        #print (dat[x_id],dat[y_id],dat[z_id],dat[seq_id],dat[asym_id], dat[comp_id], dat[atom_id])
                        pdb[(int(dat[seq_id]), dat[asym_id], dat[comp_id], dat[atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                        pdb1[(int(dat[seq_id]), dat[asym_id], dat[comp_id], dat[atom_id])] = \
                            float(dat[bf_id])
                except ValueError:
                    pass
        pdb_models[model] = pdb
        bfactor[model]=pdb1
    return pdb_models,bfactor,ent_id

def get_coordinates_pdb(cif_file, use_auth_tag=False):
    """
    Extract coordinate information from cif file as a dictionary
    {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
    :param cif_file: Input coordinate file
    :return: dictionary
    """

    cif_data = []
    if cif_file.endswith('.gz'):
        ifh = gzip.open(cif_file,'rt')
    else:
        ifh = open(cif_file, 'r')
    pRd = PdbxReader(ifh)
    pRd.read(cif_data)
    ifh.close()
    c0 = cif_data[0]
    atom_site = c0.getObj('atom_site')
    exptl=c0.getObj('exptl')
    m=exptl.getValue('method')
    reflns = c0.getObj('reflns')
    emr = c0.getObj('em_3d_reconstruction')
    ent=c0.getObj('entry')
    pdb_id = ent.getValue('id')
    refine=c0.getObj('refine')
    print (refine)
    if refine is not None:
        resolution = refine.getValue('ls_d_res_high')
        print (resolution)
    elif reflns is not None:
        resolution = reflns.getValue('d_resolution_high')
    elif emr is not None:
        resolution = emr.getValue('resolution')
    else:
        resolution = '.'
    print (resolution)
    print (m,resolution,pdb_id)
    max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
    col_names = atom_site.getAttributeList()
    model_id = col_names.index('pdbx_PDB_model_num')
    x_id = col_names.index('Cartn_x')
    y_id = col_names.index('Cartn_y')
    z_id = col_names.index('Cartn_z')
    atom_id = col_names.index('label_atom_id')
    comp_id = col_names.index('label_comp_id')
    asym_id = col_names.index('label_asym_id')
    #entity_id = col_names.index('label_entity_id')
    seq_id = col_names.index('label_seq_id')
    #icode_id = col_names.index('pdbx_PDB_ins_code')
    #alt_id = col_names.index('label_alt_id')
    aut_seq_id = col_names.index('auth_seq_id')
    aut_asym_id = col_names.index('auth_asym_id')
    aut_atom_id = col_names.index('auth_atom_id')
    aut_comp_id = col_names.index('auth_comp_id')
    pdb_models = {}
    for model in range(1, max_models + 1):
        pdb = {}
        for dat in atom_site.getRowList():
            if int(dat[model_id]) == model:
                try:
                    if use_auth_tag:
                        pdb[(int(dat[aut_seq_id]), dat[aut_asym_id], dat[aut_comp_id], dat[aut_atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                    else:
                        #print (dat[x_id],dat[y_id],dat[z_id],dat[seq_id],dat[asym_id], dat[comp_id], dat[atom_id])
                        pdb[(int(dat[seq_id]), dat[asym_id], dat[comp_id], dat[atom_id])] = \
                            numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
                except ValueError:
                    pass
        pdb_models[model] = pdb
    return pdb_models,m,resolution,pdb_id

if __name__ == "__main__":
    get_coordinates_pdb('/Users/kumaran/Downloads/4w9i.cif')