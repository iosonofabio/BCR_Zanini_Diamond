# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/19
content:    Check the output of the 10X pipeline
'''
import os
import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


def load_annotations(fdn):
    samplenames = [sfdn[:-4] for sfdn in os.listdir(fdn) if sfdn.endswith('VDJ')]

    df = []
    for sn in samplenames:
        fn = fdn+sn+'-VDJ/filtered_contig_annotations.csv'
        dfi = pd.read_csv(fn, sep=',')
        dfi['samplename'] = sn
        dfi['id'] = dfi['samplename'] + '_' + dfi['contig_id']
        df.append(dfi)
    df = pd.concat(df, axis=0)
    df.set_index('id', inplace=True)
    return df


if __name__ == '__main__':

    fdn = '../../data/sequencing/'
    df = load_annotations(fdn)
    samplenames = list(df['samplename'].unique())

    print('Count clonotypes')
    clonotypes_u = df['raw_clonotype_id'].unique()
    print('# clonal families: {:}'.format(len(clonotypes_u) - 1))

    print('Plot histogram of heavy and light chain clonotype sizes')
    kinds = ['IGH', 'IGK', 'IGL']
    fig, axs = plt.subplots(1, 3, figsize=(10, 4), sharex=True, sharey=True)
    colors = sns.color_palette('Set1', n_colors=3)
    for isn, sn in enumerate(samplenames):
        ax = axs[isn]
        if 'M_GV' in sn:
            title = 'uninfected_control'
        else:
            title = sn
        ax.set_title(title)
        dfs = df.loc[df['samplename'] == sn]
        for ik, kind in enumerate(kinds):
            dfi = dfs.loc[dfs['c_gene'].str.startswith(kind)]
            sizes = dfi.groupby('raw_clonotype_id').count().iloc[:, 0]
            sizes.drop('None', inplace=True)

            # rank plot
            y = np.sort(sizes.values)[::-1]
            rank = np.arange(len(y)) + 1

            ax.plot(rank, y, lw=2, color=colors[ik], label=kind)

        ax.set_ylabel('Clone sizes')
        ax.set_xlabel('Clone rank')
        ax.legend(title='Chain type:', loc='upper right')
        ax.set_xlim(left=0.9)
        ax.grid(True)
        ax.set_xscale('log')
        ax.set_yscale('log')
    fig.tight_layout()


    print('Find the names of the big clones')
    sn = 'anti-CD137_7dpi'
    dfs = df.loc[df['samplename'] == sn]
    big_clonotypes = {}
    sizes_kinds = []
    for kind in kinds:
        dfi = dfs.loc[dfs['c_gene'].str.startswith(kind)]
        sizes = dfi.groupby('raw_clonotype_id').count().iloc[:, 0]
        sizes.drop('None', inplace=True)
        sizes.name = kind
        sizes_kinds.append(sizes)
        big_clonotypes[kind] = sizes.nlargest(2).index
    sizes_kinds = pd.concat(sizes_kinds, axis=1, sort=False).fillna(0).astype(int)

    # NOTE: the names do not match, so they could be different?

    #print('Find the cells within the big clones')
    #from collections import defaultdict
    #clonotypes_cells = defaultdict(list)
    #for kind in kinds:
    #    for clonotype in big_clonotypes[kind]:
    #        ind = (dfs['raw_clonotype_id'] == clonotype) & (dfs['c_gene'].str.startswith(kind))
    #        bcs = dfs.loc[ind, 'barcode'].unique()
    #        for bc in bcs:
    #            clonotypes_cells[bc].append((kind, clonotype))
    #clonodouble = {} 
    #for kp in 

    plt.ion()
    plt.show()
