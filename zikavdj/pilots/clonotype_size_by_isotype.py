# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/06/19
content:    Check clonotype size by isotype
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

    print('Annotate isotypes')
    df['isotype'] = df['c_gene'].str.split('*').str[0]
    isotypes = [
            'IGHM', 'IGHD', 'IGHG1', 'IGHG2A', 'IGHG2B',
            'IGHG2C', 'IGHG3', 'IGHA', 'IGHE',
            'IGKC', 'IGLC1', 'IGLC2', 'IGLC3', 'IGLC4',
            ]
    colors = sns.color_palette('Set1', n_colors=len(isotypes))
    # Swap some colors to avoid repeats
    tmp = colors[isotypes.index('IGLC3')]
    colors[isotypes.index('IGLC3')] = colors[isotypes.index('IGHA')]
    colors[isotypes.index('IGHA')] = tmp

    print('Filter out T cells and other contaminants')
    df = df.loc[df['isotype'].str.startswith('IG')]
    df = df.loc[df['isotype'] != 'None']

    print('Singletons versus families')
    kinds = ['IGH', 'IGK', 'IGL']
    for isn, sn in enumerate(samplenames):
        dfs = df.loc[df['samplename'] == sn]
        for ik, kind in enumerate(kinds):
            dfi = dfs.loc[dfs['c_gene'].str.startswith(kind)]
            sizes = dfi.groupby('raw_clonotype_id').count().iloc[:, 0]
            if 'None' in sizes.index:
                sizes.drop('None', inplace=True)
            n_singletons = sizes.loc[sizes == 1].sum()
            n_2p = sizes.loc[sizes > 1].sum()
            fr_singletons = 1.0 * n_singletons / (n_singletons + n_2p)
            fr_2p = 1.0 * n_2p / (n_singletons + n_2p)
            print(sn, kind, 'singletons / in_families: {:} / {:} ({:.0%} / {:.0%})'.format(
                n_singletons, n_2p, fr_singletons, fr_2p))


    print('Plot clone sizes by sample and isotype')
    isotypes = df['isotype'].unique()
    kinds = isotypes
    fig, axs = plt.subplots(1, 3, figsize=(10, 4), sharex=True, sharey=True)
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
            if 'None' in sizes.index:
                sizes.drop('None', inplace=True)

            # rank plot
            y = np.sort(sizes.values)[::-1]
            rank = np.arange(len(y)) + 1

            if kind.startswith('IGH'):
                ls = '-'
            else:
                ls = '--'

            ax.plot(rank, y, lw=2, color=colors[ik], label=kind, ls=ls)

        ax.set_ylabel('Clone sizes')
        ax.set_xlabel('Clone rank')
        ax.legend(title='Chain type:', loc='upper right', ncol=2, fontsize=8)
        ax.set_xlim(left=0.9)
        ax.grid(True)
        ax.set_xscale('log')
        ax.set_yscale('log')
    fig.tight_layout()

    print('Make clone graphs')
    graphs = {}
    kinds = ['IGH', 'IGK', 'IGL']
    for isn, sn in enumerate(samplenames):
        dfs = df.loc[df['samplename'] == sn]
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
                if nct > 500:
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
            graphs[(sn, kind)] = {
                    'edges': edges,
                    'nodes': nodes,
                    'node_ids': node_ids,
                    }

    print('Calculate graph layout of clonotypes with > 1 sequence')
    import igraph as ig
    #layout_alg = 'drl'
    #layout_alg = 'grid_fr'  # FIXME: not implemented??
    layout_alg = 'fr'
    layout_kwargs = {}
    kinds = ['IGH', 'IGK', 'IGL']
    layouts = {}
    for isn, sn in enumerate(samplenames):
        for kind in kinds:
            # NOTE: by definition there are no singletons in this graph
            g = ig.Graph(graphs[(sn, kind)]['edges'])
            if 'grid' not in layout_alg:
                lo = g.layout(layout_alg, **layout_kwargs)
            else:
                lo = g.layout_fruchterman_reingold(grid=True)
            layouts[(sn, kind)] = np.array(lo.coords)

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
            3, 4,
            figsize=(12, 10),
            gridspec_kw={'width_ratios': [1, 1, 1, 0.3]},
            )
    for isn, sn in enumerate(samplenames):
        for ik, kind in enumerate(kinds):
            ax = axs[ik, isn]
            #ax.set_axis_off()

            # scatter vertices
            lo = layouts[(sn, kind)]

            graph = graphs[(sn, kind)]
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
                ax.set_title(t)
            if isn == 0:
                ax.set_ylabel(kind, rotation=0, labelpad=20)
                ax.set_yticks([])
            else:
                ax.get_yaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            if isn == len(samplenames) - 1:
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

    plt.ion()
    plt.show()
