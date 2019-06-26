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

    df['isotype'] = df['c_gene'].str.split('*').str[0]
    df = df.loc[df['isotype'].str.startswith('IG')]
    df = df.loc[df['isotype'] != 'None']

    return df


if __name__ == '__main__':

    fdn = '../../data/sequencing/'
    df = load_annotations(fdn)
    samplenames = list(df['samplename'].unique())
    isotypes = [
            'IGHM', 'IGHD', 'IGHG1', 'IGHG2A', 'IGHG2B',
            'IGHG2C', 'IGHG3', 'IGHA', 'IGHE',
            'IGKC', 'IGLC1', 'IGLC2', 'IGLC3', 'IGLC4',
            ]

    def get_vj_data(df, vj):
        vjm = vj.lower()
        col = vjm+'_gene'
        dfi = df[['samplename', col, 'isotype', 'barcode']].groupby(
                ['samplename', col, 'isotype'],
                ).count()['barcode'].unstack(col).fillna(0).astype(int)
        del dfi['None']
        dfi = dfi.stack().unstack('samplename').fillna(0).astype(int)
        return dfi

    n_vs = get_vj_data(df, vj='V')
    n_js = get_vj_data(df, vj='J')

    print('Plot heatmap of V genes')
    kinds = ['IGH', 'IGK', 'IGL']
    fig, axs = plt.subplots(
            3, 3, figsize=(26, 10),
            gridspec_kw={'width_ratios': [1, 1, 0.05]},
            )
    for isn, sn in enumerate(samplenames):
        data = n_vs[sn].unstack('isotype')
        for ik, kind in enumerate(kinds):
            ax = axs[isn, ik]
            datai = data.loc[data.index.str.startswith(kind)]
            datai = datai.loc[:, datai.columns.str.startswith(kind)]
            datai.index = [x[2:] for x in datai.index]
            sns.heatmap(datai.T, ax=ax, xticklabels=True, yticklabels=True)
            if ik == 0:
                if sn == 'M_GV-Na_ve-Na_ve':
                    t = 'uninfected'
                else:
                    t = sn
                ax.text(-0.085, 0.5, t, transform=ax.transAxes, clip_on=False, ha='right')
            ax.set_ylabel('')
            for tk in ax.get_xticklabels():
                tk.set_fontsize(8)

    fig.suptitle('V gene usage')
    fig.tight_layout(rect=(0.01, 0, 0.98, 0.95), w_pad=0, h_pad=0.1)

    print('Plot heatmap of J genes')
    kinds = ['IGH', 'IGK', 'IGL']
    fig, axs = plt.subplots(
            3, 3, figsize=(11, 10),
            gridspec_kw={'width_ratios': [1, 0.7, 0.55]},
            )
    for isn, sn in enumerate(samplenames):
        data = n_js[sn].unstack('isotype')
        for ik, kind in enumerate(kinds):
            ax = axs[isn, ik]
            datai = data.loc[data.index.str.startswith(kind)]
            datai = datai.loc[:, datai.columns.str.startswith(kind)]
            datai.index = [x[2:] for x in datai.index]
            sns.heatmap(datai.T, ax=ax, xticklabels=True, yticklabels=True)
            if ik == 0:
                if sn == 'M_GV-Na_ve-Na_ve':
                    t = 'uninfected'
                else:
                    t = sn
                ax.text(-0.35, 0.5, t, transform=ax.transAxes, clip_on=False, ha='right')
            ax.set_ylabel('')

    fig.suptitle('J gene usage')
    fig.tight_layout(rect=(0.03, 0, 0.98, 0.95), w_pad=0, h_pad=0.1)


    plt.ion()
    plt.show()
