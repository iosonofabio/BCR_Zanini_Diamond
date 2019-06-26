# vim: fdm=indent
'''
author:     Fabio Zanini
date:       13/06/19
content:    Load mRNA gene expression from 10X Genomics HDF5 files.
'''
import os
import sys
import numpy as np
import pandas as pd
import loompy
import collections
import scipy.sparse as sp_sparse
import tables

import matplotlib.pyplot as plt
import seaborn as sns


def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        barcodes = np.array([x.decode() for x in barcodes])

        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()

        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        # Convert to dense
        matrix = matrix.todense().astype(np.float32)

        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_ref['id'] = np.array([x.decode() for x in feature_ids])
        feature_ref['name'] = np.array([x.decode() for x in feature_names])

        sample_ref = {'barcode': barcodes}

        return (feature_ref, sample_ref, matrix)


if __name__ == '__main__':

    fdn = '../../data/sequencing/'
    samplenames = ['anti-CD137_7dpi', 'isotype_control', 'M_GV-Na_ve-Na_ve']

    for sn in samplenames:
        print(sn)

        print('Load from H5')
        fn_filtered_matrix = fdn+sn+'-mRNA/filtered_feature_bc_matrix.h5'
        ra, ca, matrix = get_matrix_from_h5(fn_filtered_matrix)

        print('Save to loom')
        fn_filtered_loom = fdn+sn+'-mRNA/filtered_feature_bc_matrix.loom'
        loompy.create(fn_filtered_loom, matrix, ra, ca)
