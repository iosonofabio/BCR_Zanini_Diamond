# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/06/19
content:    Try out a bunch of figure panels
'''
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd
import argparse
import igraph as ig

import matplotlib.pyplot as plt
import seaborn as sns
import loompy

sys.path.append('/home/fabio/university/postdoc/singlet')
os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
import singlet
from singlet import config as singlet_config


fig_fdn = '../../figures/for_paper/'


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


def get_clonality_graph(df, dn, samplesheet, cats=['plasma'], max_cells=500):
    print('Load VDJ info')
    sn = datasetd[dn]
    df = df.loc[df['samplename'] == sn]
    vdj_barcodes = df['barcode'].unique()

    print('Take intersection (most are there)')
    dfd = {}
    if 'proliferating' in cats:
        cycling_B_cells = samplesheet.query('cellSubtype == "Proliferating"').index
        cycling_good_B_cells = np.intersect1d(cycling_B_cells, vdj_barcodes)
        dfgb = df.loc[df['barcode'].isin(cycling_good_B_cells)]
        dfd['proliferating'] = dfgb

    if 'plasma' in cats:
        plasmablasts = samplesheet.query('cellSubtype == "Plasmablasts"').index
        good_plasmablasts = np.intersect1d(plasmablasts, vdj_barcodes)
        dfpb = df.loc[df['barcode'].isin(good_plasmablasts)]
        dfd['plasma'] = dfpb

    if 'other' in cats:
        naive_memory = samplesheet.query('cellSubtype == "Naive_memory"').index
        good_naive = np.intersect1d(naive_memory, vdj_barcodes)
        dfnm = df.loc[df['barcode'].isin(good_naive)]
        dfd['other'] = dfnm

    print('Make clone graphs')
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
                if nct > max_cells:
                    break
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

    return {
        'graphs': graphs,
        'layouts': layouts,
        }


def get_clonality_graphs(df, samplesheets, **kwargs):
    graphs = {}
    layouts = {}
    for dn, samplesheet in samplesheets.items():
        res = get_clonality_graph(df, dn, samplesheet, **kwargs)
        graphs.update(res['graphs'])
        layouts.update(res['layouts'])
    return {
        'graphs': graphs,
        'layouts': layouts,
        }


if __name__ == '__main__':

    datasetnames = ['cd137', 'isotype_control', 'uninfected']
    datasetd = {
            'cd137': 'anti-CD137_7dpi',
            'isotype_control': 'isotype_control',
            'uninfected': 'M_GV-Na_ve-Na_ve',
            }

    pa = argparse.ArgumentParser()

    print('Load samplesheets')
    samplesheets = {}
    for dsn in datasetd:
        loomfn = singlet_config.config['io']['datasets'][dsn+'_cpm']['path']
        with loompy.connect(loomfn) as dsl:
            cid = dsl.ca['barcode']
            ct = dsl.ca['cellType']
            cst = dsl.ca['cellSubtype']
        df = pd.DataFrame([], index=cid)
        df['cellType'] = ct
        df['cellSubtype'] = cst
        samplesheets[dsn] = df


    print('Calculate fractions')
    fractions = pd.DataFrame([], columns=['plasma', 'proliferating', 'other', 'total_number'])
    for dsn in datasetd:
        ntot = (samplesheets[dsn]['cellType'] == 'B cell').sum()
        n_plasma = (samplesheets[dsn]['cellSubtype'] == 'Plasmablasts').sum()
        n_prolif = (samplesheets[dsn]['cellSubtype'] == 'Proliferating').sum()
        n_other = ntot - n_plasma - n_prolif
        fractions.loc[dsn] = [
                1. * n_plasma / ntot,
                1. * n_prolif / ntot,
                1. * n_other / ntot,
                ntot,
                ]

    print('Plot fractions')
    fig, ax = plt.subplots(figsize=(2, 3))
    order = ['uninfected', 'isotype_control', 'cd137']
    colors = {'plasma': 'tomato', 'proliferating': 'steelblue'}
    for i, cat in enumerate(['plasma', 'proliferating']):
        xl = 0.3 * (-1 + i) + np.arange(3)
        w = 0.33
        y = 100. * fractions[cat].loc[order].values
        ax.bar(
            xl, y, w, align='edge',
            color=colors[cat],
            label=cat.capitalize(),
            zorder=10,
            )
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(
            ['Uninfected\nctrl', 'Isotype\nctrl', 'Anti-CD137'],
            rotation=90, ha='center',
            )
    ax.set_ylabel('% of B cells')
    ax.grid(True, axis='y')
    ax.legend(loc='lower right', bbox_to_anchor=(1.04, 1.01), bbox_transform=ax.transAxes)
    plt.tight_layout()

    fig.savefig(fig_fdn+'fractions.png', dpi=300)
    fig.savefig(fig_fdn+'fractions.pdf')
    fig.savefig(fig_fdn+'fractions.svg')

    if False:
        print('Load {:} data from loom file'.format(dn))
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

    fdn = '../../data/sequencing/'
    df = load_vdj_annotations(fdn)
    res = get_clonality_graphs(df, samplesheets)
    graphs = res['graphs']
    layouts = res['layouts']

    print('Plot pretty graphs of clonotypes with > 1 sequence')
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
            1, 4,
            figsize=(10, 3),
            gridspec_kw={'width_ratios': [1, 1, 1, 0.3]},
            )
    cellt = 'plasma'
    kind = 'IGH'
    for isn, dn in enumerate(order):
        sn = datasetd[dn]
        ax = axs[isn]

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
        if sn == 'M_GV-Na_ve-Na_ve':
            t = 'uninfected'
        else:
            t = sn
        ax.set_title(t)
        ax.set_yticks([])
        ax.get_xaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        if isn == 2:
            from matplotlib.patches import Rectangle
            ax = axs[-1]
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

    fig.savefig(fig_fdn+'balls_plasma.png', dpi=300)
    fig.savefig(fig_fdn+'balls_plasma.pdf')
    fig.savefig(fig_fdn+'balls_plasma.svg')

    print('Study entropy within the clusters')
    kind = 'IGH'
    cellt = 'plasma'
    fracs = {}
    for dn in order:
        sn = datasetd[dn]
        graph = graphs[(sn, cellt, kind)]

        # Figure out the connected components
        cc = {}
        for e1, e2 in graph['edges']:
            if e1 in cc:
                cc[e2] = cc[e1]
            else:
                cc[e1] = cc[e2] = len(cc)
        # Invert dictionary
        cci = [[] for x in range(len(cc))]
        for v, c in cc.items():
            cci[c].append(v)

        # Find fraction of families that are not homogeneous
        n_nh = 0
        for cf in cci:
            its = set()
            for v in cf:
                nid = graph['node_ids'][v]
                if 'center' in nid:
                    continue
                its.add(df.at[nid, 'isotype'])
            if len(its) > 1:
                n_nh += 1
        frac = 1.0 * n_nh / len(cci)
        fracs[dn] = frac
    fracs = pd.Series(fracs).loc[order]

    print('Plot entropy of isotypes')
    fig, ax = plt.subplots(figsize=(2.2, 3))
    ax.bar(-0.35 + np.arange(3), 100. * fracs.values, 0.7, align='edge', zorder=10)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(
            ['Uninfected\nctrl', 'Isotype\nctrl', 'Anti-CD137'],
            rotation=90, ha='center',
            )
    ax.set_ylabel('heterotypic CFs\n[% of all CFs]')
    ax.grid(True, axis='y')
    fig.tight_layout(rect=(0.03, 0, 1, 1))

    fig.savefig(fig_fdn+'plasmablast_fraction_heterotypic_CFs.png', dpi=300)
    fig.savefig(fig_fdn+'plasmablast_fraction_heterotypic_CFs.pdf')
    fig.savefig(fig_fdn+'plasmablast_fraction_heterotypic_CFs.svg')

    print('Get VJ usage')
    # Restrict to real B cells
    idx = []
    for dn, samplesheet in samplesheets.items():
        sn = datasetd[dn]
        bcells = samplesheet.query('cellType == "B cell"').index
        idxi = df.query(
            '(samplename == @sn) & (barcode in @bcells)',
            local_dict=locals(),
            ).index
        idx.extend(list(idxi))
    dfb = df.loc[idx]

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

    n_vs = get_vj_data(dfb, vj='V')
    n_js = get_vj_data(dfb, vj='J')

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
            ax.set_ylim(0, datai.shape[1])
            for tk in ax.get_yticklabels():
                tk.set_rotation(0)
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

    fig.savefig(fig_fdn+'V_usage_with_IgG_subtypes.png', dpi=300)
    fig.savefig(fig_fdn+'V_usage_with_IgG_subtypes.pdf')
    fig.savefig(fig_fdn+'V_usage_with_IgG_subtypes.svg')

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
            ax.set_ylim(0, datai.shape[1])
            for tk in ax.get_yticklabels():
                tk.set_rotation(0)
            if ik == 0:
                if sn == 'M_GV-Na_ve-Na_ve':
                    t = 'uninfected'
                else:
                    t = sn
                ax.text(-0.35, 0.5, t, transform=ax.transAxes, clip_on=False, ha='right')
            ax.set_ylabel('')

    fig.suptitle('J gene usage')
    fig.tight_layout(rect=(0.03, 0, 0.98, 0.95), w_pad=0, h_pad=0.1)

    fig.savefig(fig_fdn+'J_usage_with_IgG_subtypes.png', dpi=300)
    fig.savefig(fig_fdn+'J_usage_with_IgG_subtypes.pdf')
    fig.savefig(fig_fdn+'J_usage_with_IgG_subtypes.svg')

    print('Coarse grain IgG isotypes')
    n_vs_cg = n_vs.copy().unstack()
    n_vs_cg.loc['IGHG'] = 0
    for isot in n_vs_cg.index:
        if isot.startswith('IGHG') and (isot != 'IGHG'):
            n_vs_cg.loc['IGHG'] += n_vs_cg.loc[isot]
    itypes = ['IGHA', 'IGHD', 'IGHE', 'IGHG', 'IGHM', 'IGKC', 'IGLC1', 'IGLC2', 'IGLC3', 'IGLC4']
    n_vs_cg = n_vs_cg.loc[itypes].stack()

    n_js_cg = n_js.copy().unstack()
    n_js_cg.loc['IGHG'] = 0
    for isot in n_js_cg.index:
        if isot.startswith('IGHG') and (isot != 'IGHG'):
            n_js_cg.loc['IGHG'] += n_js_cg.loc[isot]
    itypes = ['IGHA', 'IGHD', 'IGHE', 'IGHG', 'IGHM', 'IGKC', 'IGLC1', 'IGLC2', 'IGLC3', 'IGLC4']
    n_js_cg = n_js_cg.loc[itypes].stack()

    print('Plot heatmap of V genes, coarse grained')
    kinds = ['IGH', 'IGK', 'IGL']
    fig, axs = plt.subplots(
            3, 3, figsize=(26, 10),
            gridspec_kw={'width_ratios': [1, 1, 0.05]},
            )
    for isn, sn in enumerate(samplenames):
        data = n_vs_cg[sn].unstack('isotype')
        for ik, kind in enumerate(kinds):
            ax = axs[isn, ik]
            datai = data.loc[data.index.str.startswith(kind)]
            datai = datai.loc[:, datai.columns.str.startswith(kind)]
            datai.index = [x[2:] for x in datai.index]
            sns.heatmap(datai.T, ax=ax, xticklabels=True, yticklabels=True)
            ax.set_ylim(0, datai.shape[1])
            for tk in ax.get_yticklabels():
                tk.set_rotation(0)
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

    fig.savefig(fig_fdn+'V_usage_without_IgG_subtypes.png', dpi=300)
    fig.savefig(fig_fdn+'V_usage_without_IgG_subtypes.pdf')
    fig.savefig(fig_fdn+'V_usage_without_IgG_subtypes.svg')

    print('Plot heatmap of J genes, coarse grained')
    kinds = ['IGH', 'IGK', 'IGL']
    fig, axs = plt.subplots(
            3, 3, figsize=(11, 10),
            gridspec_kw={'width_ratios': [1, 0.7, 0.55]},
            )
    for isn, sn in enumerate(samplenames):
        data = n_js_cg[sn].unstack('isotype')
        for ik, kind in enumerate(kinds):
            ax = axs[isn, ik]
            datai = data.loc[data.index.str.startswith(kind)]
            datai = datai.loc[:, datai.columns.str.startswith(kind)]
            datai.index = [x[2:] for x in datai.index]
            sns.heatmap(datai.T, ax=ax, xticklabels=True, yticklabels=True)
            ax.set_ylim(0, datai.shape[1])
            for tk in ax.get_yticklabels():
                tk.set_rotation(0)
            if ik == 0:
                if sn == 'M_GV-Na_ve-Na_ve':
                    t = 'uninfected'
                else:
                    t = sn
                ax.text(-0.35, 0.5, t, transform=ax.transAxes, clip_on=False, ha='right')
            ax.set_ylabel('')

    fig.suptitle('J gene usage')
    fig.tight_layout(rect=(0.03, 0, 0.98, 0.95), w_pad=0, h_pad=0.1)

    fig.savefig(fig_fdn+'J_usage_without_IgG_subtypes.png', dpi=300)
    fig.savefig(fig_fdn+'J_usage_without_IgG_subtypes.pdf')
    fig.savefig(fig_fdn+'J_usage_without_IgG_subtypes.svg')

    plt.ion()
    plt.show()
