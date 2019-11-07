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
    print('Load clusters from file')
    ds.samplesheet['community'] = pd.read_csv(
        fn_anno,
        sep='\t', index_col=0)['community']

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


    print('Functions to assign cell types')
    ds.samplesheet['cellType'] = ''
    def assign_ct(cns, ct):
        ds.samplesheet.loc[ds.samplesheet['community_rough'].isin(cns), 'cellType'] = ct
    ds.samplesheet['cellSubtype'] = ''
    def assign_cst(cns, cst):
        ds.samplesheet.loc[ds.samplesheet['community_rough'].isin(cns), 'cellSubtype'] = cst

    def assign_tsne_fig(cst):
        x0, x1 = plt.gcf().get_axes()[0].get_xlim()
        y0, y1 = plt.gcf().get_axes()[0].get_ylim()
        print(x0, x1, y0, y1)
        ind = (vs['dimension 1'] >= x0) & (vs['dimension 1'] < x1)
        ind &= (vs['dimension 2'] >= y0) & (vs['dimension 2'] < y1)
        ind = vs.index[ind]
        ds.samplesheet.loc[ind, 'cellSubtype'] = cst

    print('Plot few markers onto embedding')
    vs = vst
    marker_genes = [
        'Tmcc2',
        'Cd48',
        'Aqp1',
        ('Ms4a1', 'Cd20'),
        'Cd19',
        'Cd79b',
        'Jchain',
        'Prdm1',
        'Mki67',
        'Aicda',
        'S1pr2',
        'Mef2b',
        'Cdc123',
        'Ighm',
        'Cd3e',
        'Cd14',
        ]
    markers = marker_genes + [
        ('community', 'cluster'),
        ('community_rough', 'cell type'),
        ]
    fig, axs = plt.subplots(2, 9, figsize=(15, 4), sharex=True, sharey=True)
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
    #fig.savefig('../../figures/tsne_{:}_markers.png'.format(dn))

    plt.ion()
    plt.show()

    if False:
        print('Save cell type annotations to loom file')
        import loompy
        for tmp in ['', '_cpm']:
            loom_fn = singlet.config.config['io']['datasets'][args.sample+tmp]['path']
            with loompy.connect(loom_fn) as dsl:
                dsl.ca['cellType'] = ds.samplesheet['cellType'].values
                dsl.ca['cellSubtype'] = ds.samplesheet['cellSubtype'].values
                dsl.ca['tSNE1'] = vs.iloc[:, 0].values
                dsl.ca['tSNE2'] = vs.iloc[:, 1].values


