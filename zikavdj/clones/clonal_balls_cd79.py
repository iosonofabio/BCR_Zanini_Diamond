# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/06/19
content:    Explore mRNA datasets
'''
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd
import argparse

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
import singlet


def load_vdj_annotations(fdn):
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

    datasetnames = ['cd137', 'isotype_control', 'uninfected']
    datasetd = {
            'cd137': 'anti-CD137_7dpi',
            'isotype_control': 'isotype_control',
            'uninfected': 'M_GV-Na_ve-Na_ve',
            }

    pa = argparse.ArgumentParser()
    pa.add_argument('--sample', required=True, choices=datasetnames)
    args = pa.parse_args()
    dn = args.sample

    print('Load {:} data from loom file'.format(dn))
    sn = datasetd[dn]
    ds = singlet.Dataset(dataset=dn+'_cpm')
    # FIXME
    ds.counts = singlet.CountsTable(ds.counts)
    ds.counts._normalized = 'counts_per_million'

    print('Exclude genes with more than one id')
    from collections import Counter
    cou = Counter(ds.featuresheet['name'].values)
    good_genes = [key for key, val in cou.items() if val == 1]
    ind = ds.featuresheet['name'].isin(good_genes)
    good_ids = ds.featurenames[ind]
    ds.query_features_by_name(good_ids, inplace=True)

    print('Convert to gene names')
    ds.reindex('features', 'name', inplace=True)

    print('Get list of cycling B cells, plasmablasts, and naive/memory')
    dsB = ds.query_samples_by_metadata('cellType == "B cell"')
    cycling_B_cells = dsB.samplenames[dsB.samplesheet['cellSubtype'] == 'Proliferating']
    plasmablasts = dsB.samplenames[dsB.samplesheet['cellSubtype'] == 'Plasmablasts']
    naive_memory = dsB.samplenames[dsB.samplesheet['cellSubtype'] == 'Naive_memory']

    print('Load VDJ info')
    fdn = '../../data/sequencing/'
    df = load_vdj_annotations(fdn)
    df = df.loc[df['samplename'] == sn]
    vdj_barcodes = df['barcode'].unique()

    print('Take intersection (most are there)')
    cycling_good_B_cells = np.intersect1d(cycling_B_cells, vdj_barcodes)
    dsgb = dsB.query_samples_by_name(cycling_good_B_cells)
    dfgb = df.loc[df['barcode'].isin(cycling_good_B_cells)]

    if len(plasmablasts):
        good_plasmablasts = np.intersect1d(plasmablasts, vdj_barcodes)
        dspb = dsB.query_samples_by_name(good_plasmablasts)
        dfpb = df.loc[df['barcode'].isin(good_plasmablasts)]

    good_naive = np.intersect1d(naive_memory, vdj_barcodes)
    dsnm = dsB.query_samples_by_name(good_naive)
    dfnm = df.loc[df['barcode'].isin(good_naive)]

    print('Make clone graphs')
    cell_types = ['cycling', 'plasmablasts']
    dfd = {'cycling': dfgb}
    if len(plasmablasts):
        dfd['plasmablasts'] = dfpb
    dfd['naive_memory'] = dfnm
    graphs = {}
    kinds = ['IGH', 'IGK', 'IGL']
    for cellt in dfd.keys():
        dfs = dfd[cellt]
        for kind in kinds:
            dfi = dfs.loc[dfs['c_gene'].str.startswith(kind)]
            sizes = dfi.groupby('raw_clonotype_id').count().iloc[:, 0]
            clones_2p = sizes.index[(sizes >= 2) & (sizes.index != 'None')]

            if len(clones_2p) == 0:
                print(sn, kind, 'no clones')
                continue

            dfi = dfi.loc[dfi['raw_clonotype_id'].isin(clones_2p)]
            noden = 0
            edges = []
            nodes = []
            node_ids = []
            nct = 0
            for ct, data in dfi.groupby('raw_clonotype_id'):
                # FIXME: only for testing
                if nct > 200:
                    break
                #if len(data) > 30:
                #    continue

                nodes.append(noden)
                node_ids.append('clonotype_center_{:}_{:}_{:}'.format(sn, kind, ct))

                for ir, (sid, row) in enumerate(data.iterrows()):
                    edges.append((noden, noden + ir + 1))
                    nodes.append(noden + ir + 1)
                    node_ids.append(sid)
                noden += len(data) + 1
                nct += 1
            graphs[(sn, cellt, kind)] = {
                    'edges': edges,
                    'nodes': nodes,
                    'node_ids': node_ids,
                    }

    print('Calculate graph layout of clonotypes with > 1 sequence')
    import igraph as ig
    layout_alg = 'fr'
    layout_kwargs = {}
    kinds = ['IGH', 'IGK', 'IGL']
    layouts = {}
    for cellt in dfd.keys():
        for kind in kinds:
            # NOTE: by definition there are no singletons in this graph
            if (sn, cellt, kind) not in graphs:
                continue
            g = ig.Graph(graphs[(sn, cellt, kind)]['edges'])
            if 'grid' not in layout_alg:
                lo = g.layout(layout_alg, **layout_kwargs)
            else:
                lo = g.layout_fruchterman_reingold(grid=True)
            layouts[(sn, cellt, kind)] = np.array(lo.coords)

    print('Plot pretty graphs of clonotypes with > 1 sequence')
    kinds = ['IGH', 'IGK', 'IGL']
    isotypes = [
            'IGHM', 'IGHD', 'IGHG1', 'IGHG2A', 'IGHG2B',
            'IGHG2C', 'IGHG3', 'IGHA', 'IGHE',
            'IGKC', 'IGLC1', 'IGLC2', 'IGLC3', 'IGLC4',
            ]
    colors = sns.color_palette('Set1', n_colors=9)
    colors_d = {
        'IGHM': colors[1],
        'IGHD': colors[3],
        'IGHG1': colors[8],
        'IGHG2A': colors[2],
        'IGHG2B': colors[7],
        'IGHG2C': colors[6],
        'IGHG3': colors[0],
        'IGHA': colors[5],
        'IGHE': colors[4],
        'IGKC': colors[2],
        'IGLC1': colors[0],
        'IGLC2': colors[4],
        'IGLC3': colors[5],
        'IGLC4': colors[1],
        }
    fig, axs = plt.subplots(
            3, 3 + bool(len(plasmablasts)),
            figsize=(7.7 + 3.3 * bool(len(plasmablasts)), 8.7),
            gridspec_kw={'width_ratios': ([1] * (2 + bool(len(plasmablasts)))) + [0.3]},
            )
    for icellt, cellt in enumerate(dfd.keys()):
        for ik, kind in enumerate(kinds):
            ax = axs[ik, icellt]

            # scatter vertices
            if (sn, cellt, kind) not in graphs:
                ax.set_axis_off()
                continue
            lo = layouts[(sn, cellt, kind)]
            graph = graphs[(sn, cellt, kind)]
            sdata = []
            for iid, nid in enumerate(graph['node_ids']):
                sd = {'number': iid}
                if nid.startswith('clonotype_center'):
                    sd['s'] = 3
                    sd['c'] = 'k'
                else:
                    sd['s'] = 15
                    sd['c'] = colors_d[df.at[nid, 'isotype']]

                sdata.append(sd)
            sdata = pd.DataFrame(sdata)
            s = sdata['s']
            c = sdata['c']

            ax.scatter(lo[:, 0], lo[:, 1], s=s, c=c, alpha=0.8)
            # plot edges
            for n1, n2 in graph['edges']:
                ax.plot(
                    lo[[n1, n2], 0], lo[[n1, n2], 1],
                    lw=1, color='k', alpha=0.6,
                    )
            if ik == 0:
                if sn == 'M_GV-Na_ve-Na_ve':
                    t = 'uninfected'
                else:
                    t = sn
                ax.set_title(t+', '+cellt)
            ax.set_ylabel(kind, rotation=0, labelpad=20)
            ax.set_yticks([])
            ax.get_xaxis().set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            if icellt == len(dfd.keys()) - 1:
                from matplotlib.patches import Rectangle
                for ik, kind in enumerate(kinds):
                    ax = axs[ik, -1]
                    ax.set_axis_off()
                    handles = []
                    isotypes_kind = []
                    for isotype, color in colors_d.items():
                        if isotype.startswith(kind):
                            isotypes_kind.append(isotype)
                            r = Rectangle((0.1, 0.1), 0.2, 0.6, ec='none', fc=color)
                            handles.append(r)
                    ax.legend(handles, isotypes_kind, loc='upper left', title='Isotype:')
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, len(isotypes_kind))

    fig.tight_layout()
    fig.savefig('../../figures/pretty_graph_{:}_specific_subtypes_Bcellonly.png'.format(dn))

    plt.ion()
    plt.show()
