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
    pa.add_argument('--normalized', action='store_true')
    args = pa.parse_args()
    dn = args.sample

    print('Load {:} data from loom file'.format(dn))
    sn = datasetd[dn]
    if args.normalized:
        ds = singlet.Dataset(dataset=dn+'_cpm')
        # FIXME
        ds.counts = singlet.CountsTable(ds.counts)
        ds.counts._normalized = 'counts_per_million'
    else:
        ds = singlet.Dataset(dataset=dn)

        ds.samplesheet['coverage'] = ds.counts.sum(axis=0)

        print('Normalize cpm')
        normg = ds.samplesheet['coverage'] / 1000000
        ds.counts = singlet.CountsTable(ds.counts / normg)
        ds.counts._normalized = 'counts_per_million'

        print('Save to normalized loom file')
        ds.to_dataset_file('../../data/sequencing/{:}-mRNA/filtered_feature_bc_matrix_normalized.loom'.format(sn))

    print('Exclude genes with more than one id')
    from collections import Counter
    cou = Counter(ds.featuresheet['name'].values)
    good_genes = [key for key, val in cou.items() if val == 1]
    ind = ds.featuresheet['name'].isin(good_genes)
    good_ids = ds.featurenames[ind]
    ds.query_features_by_name(good_ids, inplace=True)

    print('Convert to gene names')
    ds.reindex('features', 'name', inplace=True)

    print('Feature selection')
    features = ds.feature_selection.overdispersed(500, inplace=False)
    dsf = ds.query_features_by_name(features)

    print('PCA')
    dsc = dsf.dimensionality.pca(n_dims=25, robust=False, return_dataset='samples')

    print('Knn graph')
    edges0 = dsc.graph.knn(axis='samples', return_sparse=False)[0]
    edges = set()
    for edgesi in edges0:
        for edge in edgesi:
            if edge in edges:
                continue
            if (edge[1], edge[0]) in edges:
                continue
            edges.add(edge)
    edges = list(edges)

    fn_anno = '../../data/sequencing/{:}-mRNA/samplesheet_with_Leiden_community_and_coverage.tsv'.format(sn)
    if os.path.isfile(fn_anno):
        print('Load clusters from file')
        ds.samplesheet['community'] = pd.read_csv(
            fn_anno,
            sep='\t', index_col=0)['community']
    else:
        print('Unsupervised clustering')
        import igraph as ig
        sys.path.insert(0, os.path.abspath('../../packages/'))
        import leidenalg
        G = ig.Graph(edges=edges)
        partition = partition = leidenalg.CPMVertexPartition(
                G,
                resolution_parameter=0.01,
                )
        opt = leidenalg.Optimiser()
        opt.optimise_partition(partition)
        communities = partition.membership
        print('n. communities: {:}'.format(len(np.unique(communities))))
        ds.samplesheet['community'] = communities

    print('Unsupervised clustering, rough')
    import igraph as ig
    sys.path.insert(0, os.path.abspath('../../packages/'))
    import leidenalg
    G = ig.Graph(edges=edges)
    partition = partition = leidenalg.CPMVertexPartition(
            G,
            resolution_parameter=0.002,
            )
    opt = leidenalg.Optimiser()
    opt.optimise_partition(partition)
    communities = partition.membership
    n_communities = len(np.unique(communities))
    print('n. communities: {:}'.format(n_communities))
    ds.samplesheet['community_rough'] = communities

    fn_anno = '../../data/sequencing/{:}-mRNA/samplesheet_with_Leiden_community_and_coverage.tsv'.format(sn)
    if os.path.isfile(fn_anno):
        print('Load tSNE and UMAP from file')
        tmp = pd.read_csv(
            fn_anno,
            sep='\t', index_col=0)
        vst = tmp[['tsne1', 'tsne2']].rename(columns={'tsne1': 'dimension 1', 'tsne2': 'dimension 2'})
        if 'umap1' in tmp.columns:
            vsu = tmp[['umap1', 'umap2']].rename(columns={'umap1': 'dimension 1', 'umap2': 'dimension 2'})
    else:
        print('tSNE and UMAP')
        vst = dsc.dimensionality.tsne(perplexity=30)
        vsu = dsc.dimensionality.umap()

    if False:
        print('Get cell type markers')
        comp = None
        for com in np.arange(n_communities):
            print('cell type: {:}'.format(com))
            ds.samplesheet['is_celltype_{:}'.format(com)] = ds.samplesheet['community_rough'] == com
            dsp = ds.split('is_celltype_{:}'.format(com))
            if comp is None:
                comp = dsp[True].compare(dsp[False])[['P-value']]
            else:
                comp['P-value'] = dsp[True].compare(dsp[False])['P-value']
            comp['log2_fc'] = np.log2(dsp[True].counts.mean(axis=1) + 0.1) - np.log2(dsp[False].counts.mean(axis=1))
            nla = comp.loc[comp['log2_fc'] > 0].nsmallest(10, 'P-value').index
            print(nla)

    if False:
        print('Plot expression of markers')
        vs = vst
        marker_genes = [
            'Elane',
            'Fcer1g',
            'Sell',
            'Fcrla',
            'Fcmr',
            'Cd28',
            'Cd48',
            'Abcb6',
            'mt-Nd1',
            'Tmcc2',
            'Bpgm',
            'Aqp1',
            'Cst3',
            'Cxcr2',
            'Chil1',
            'Npm1',
            'Car2',
            'Tyrobp',
            'Ptprs',
            'Eng',
            'Cd3e',
            'Xbp1',
            'Ppbp',
            'Cd9',
            'C1qa',
            'C1qb',
            'Cd79a',
            'Cd79b',
            'Fcmr',
            'Cox6a2',
            'Siglech',
            'Rnase6',
            'Cd37',
            'Tagln',
            'Cxcl12',
            'Ifitm3',
            'Epcam',
            'Lyz2',
            'Cd63',
            'Gata2',
            ('Pecam1', 'Cd31'),
            'Col6a2',
            ('Ptprc', 'Cd45'),
            'Cd19',
            ('Ms4a1', 'Cd20'),
            'Ms4a2',
            ('Trac', 'Trac (T cells)'),
            ('Gzma', 'Gzma (NK cells)'),
            'Cd14',
            ('Fcgr3', 'Cd16'),
            'Ptprs',
            'Cd68',
            'Tcl1',
            'Mki67',
            'Prdm1',
            'Jchain',
            'Cd38',
            'Ighm',
            'Ighd',
            'Ighg1',
            #'Ighg2a',
            'Ighg2b',
            'Ighg2c',
            'Ighg3',
            'Igha',
            #'Ighe',
            ]
        markers = marker_genes + [
            #('community', 'cluster'),
            ('community_rough', 'cell type'),
            ]
        fig, axs = plt.subplots(6, 11, figsize=(24, 11), sharex=True, sharey=True)
        axs = axs.ravel()
        for ipl, (gene, ax) in enumerate(zip(markers, axs)):
            print('Plotting gene {:} of {:}'.format(ipl+1, len(markers)))
            if isinstance(gene, str):
                gene, title = gene, gene
            else:
                gene, title = gene
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=10,
                    alpha=0.15 + 0.2 * (gene not in ['annotated', 'community', 'community_rough', 'Mousename', 'Treatment', 'Timepoint']),
                    color_by=gene,
                    color_log=(gene not in ['annotated', 'community', 'community_rough', 'Mousename', 'Treatment', 'Time [days] + 2', 'doublet', 'Timepoint']),
                    )
            ax.grid(False)
            ax.set_title(title)

            if gene in ['community', 'community_rough']:
                for com in ds.samplesheet[gene].unique():
                    vsc = vs.loc[ds.samplesheet[gene] == com]
                    xc, yc = vsc.values.mean(axis=0)
                    ax.scatter([xc], [yc], s=10, facecolor='none', edgecolor='red', marker='^')
                    ax.text(xc, yc, str(com), fontsize=8, ha='center', va='bottom')

        fig.tight_layout()

    print('Plot few markers onto embedding')
    vs = vst
    marker_genes = [
        'Tmcc2',
        'Cd48',
        'Aqp1',
        ('Ms4a1', 'Cd20'),
        'Cd3e',
        'Prdm1',
        'Mki67',
        'Cd14',
        ]
    markers = marker_genes + [
        ('community', 'cluster'),
        ('community_rough', 'cell type'),
        ]
    fig, axs = plt.subplots(2, 5, figsize=(10, 4), sharex=True, sharey=True)
    axs = axs.ravel()
    for ipl, (gene, ax) in enumerate(zip(markers, axs)):
        print('Plotting gene {:} of {:}'.format(ipl+1, len(markers)))
        if isinstance(gene, str):
            gene, title = gene, gene
        else:
            gene, title = gene
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=10,
                alpha=0.15 + 0.2 * (gene not in ['annotated', 'community', 'community_rough', 'Mousename', 'Treatment', 'Timepoint']),
                color_by=gene,
                color_log=(gene not in ['annotated', 'community', 'community_rough', 'Mousename', 'Treatment', 'Time [days] + 2', 'doublet', 'Timepoint']),
                )
        ax.grid(False)
        ax.set_title(title)

        if gene in ['community', 'community_rough']:
            for com in ds.samplesheet[gene].unique():
                vsc = vs.loc[ds.samplesheet[gene] == com]
                xc, yc = vsc.values.mean(axis=0)
                ax.scatter([xc], [yc], s=10, facecolor='none', edgecolor='red', marker='^')
                ax.text(xc, yc, str(com), fontsize=8, ha='center', va='bottom')

    fig.tight_layout()
    fig.savefig('../../figures/tsne_{:}_markers.png'.format(dn))

    cluster_numbers = {
            'cd137': {'cycling': 30, 'plasmablasts': 11},
            'isotype_control': {'cycling': 22, 'plasmablasts': 21},
            'uninfected': {'cycling': 25, 'plasmablasts': None},
            }

    cni = cluster_numbers[dn]
    cncy = cni['cycling']
    cnpb = cni['plasmablasts']

    print('Cluster {:} is the one with Mki67+ Cd20+ cells'.format(cncy))
    dsgb = ds.query_samples_by_metadata('community == @cncy', local_dict=locals())
    cycling_B_cells = dsgb.samplenames

    if cnpb is not None:
        print('Cluster {:} is the plasmablasts (roughly)'.format(cnpb))
        dspb = ds.query_samples_by_metadata('community == @cnpb', local_dict=locals())
        plasmablasts = dspb.samplenames
    else:
        print('No cluster with plasmablasts found')

    print('Load VDJ info')
    fdn = '../../data/sequencing/'
    df = load_vdj_annotations(fdn)
    df = df.loc[df['samplename'] == sn]
    vdj_barcodes = df['barcode'].unique()

    print('Take intersection (most are there)')
    cycling_good_B_cells = np.intersect1d(cycling_B_cells, vdj_barcodes)
    dsgb.query_samples_by_name(cycling_good_B_cells, inplace=True)
    dfgb = df.loc[df['barcode'].isin(cycling_good_B_cells)]

    if cnpb is not None:
        good_plasmablasts = np.intersect1d(plasmablasts, vdj_barcodes)
        dspb.query_samples_by_name(good_plasmablasts, inplace=True)
        dfpb = df.loc[df['barcode'].isin(good_plasmablasts)]

    print('Make clone graphs')
    cell_types = ['cycling', 'plasmablasts']
    dfd = {'cycling': dfgb}
    if cnpb is not None:
        dfd['plasmablasts'] = dfpb
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
            3, 2 + int(cnpb is not None),
            figsize=(4.7 + 3.3 * int(cnpb is not None), 8.7),
            gridspec_kw={'width_ratios': ([1] * (1 + int(cnpb is not None))) + [0.3]},
            )
    for icellt, cellt in enumerate(dfd.keys()):
        for ik, kind in enumerate(kinds):
            ax = axs[ik, icellt]

            # scatter vertices
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
    fig.savefig('../../figures/pretty_graph_{:}_specific_subtypes.png'.format(dn))

    plt.ion()
    plt.show()
