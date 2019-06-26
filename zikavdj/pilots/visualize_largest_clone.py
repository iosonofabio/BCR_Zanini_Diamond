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

    df['isotype'] = df['c_gene'].str.split('*').str[0]
    df['chain_kind'] = [x[:3] for x in df['isotype'].values]
    df = df.loc[df['isotype'].str.startswith('IG')]
    df = df.loc[df['isotype'] != 'None']

    return df


if __name__ == '__main__':

    fdn = '../../data/sequencing/'
    df = load_annotations(fdn)
    samplenames = list(df['samplename'].unique())

    print('Annotate isotypes')
    sn = 'anti-CD137_7dpi'
    kind = 'IGH'
    clonotype = 'clonotype1'

    dfl = df.query(
        '(samplename == @sn) & (chain_kind == @kind) & (raw_clonotype_id == @clonotype)',
        )

    # Exclude a few None in the CDR3
    dfl = dfl.loc[dfl['cdr3'] != 'None']

    print('Get the sequences')
    from Bio import SeqIO
    fn = fdn+'{:}-VDJ/filtered_contig.fasta'.format(sn)
    seqsd = {}
    for seq in SeqIO.parse(fn, 'fasta'):
        seqid = seq.id
        if seqid in dfl['contig_id'].values:
            seqsd[seqid] = str(seq.seq)
    fn_seqs = '/tmp/largest_clone_seqs.fasta'
    with open(fn_seqs, 'wt') as f:
        for seqid in dfl['contig_id'].values:
            f.write('> {:}\n{:}\n'.format(seqid, seqsd[seqid]))

    print('Align the sequences')
    import subprocess as sp
    fn_ali = '/tmp/largest_clone_ali.fasta'
    sp.run('muscle -in {:} -out {:} -diags'.format(fn_seqs, fn_ali), shell=True)

    print('Load alignment')
    from Bio import AlignIO
    ali = AlignIO.read(fn_ali, 'fasta')

    print('Trim alignment')
    # FIXME: use a better criterion
    ali = ali[:, 50:600]

    print('Reduce to unique sequences in here')
    from collections import defaultdict
    seqs_red = defaultdict(list)
    for seq in ali:
        seqs_red[str(seq.seq)].append(seq.id)
    seqs_unique = list(seqs_red.keys())
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import ambiguous_dna
    ali_unique = []
    for isu, su in enumerate(seqs_unique):
        ali_unique.append(
            SeqRecord(
                Seq(su, alphabet=ambiguous_dna),
                id='seq_'+str(isu+1),
                name='seq_'+str(isu+1),
                description='',
                ))
    ali_unique = MultipleSeqAlignment(ali_unique)

    print('Save alignment in PHYLIP formats')
    fn_ali_trim = '/tmp/largest_clone_ali_trim.phy'
    AlignIO.write(ali_unique, fn_ali_trim, 'phylip')

    print('Infer tree')
    res = sp.run(
            'phyml -i {:} -d nt'.format(fn_ali_trim),
            shell=True, stdout=sp.PIPE,
            )

    print('Rename and bush up leaves')
    fn_tree = '/tmp/largest_clone_ali_trim.phy_phyml_tree.txt'
    from Bio import Phylo
    tree = Phylo.read('/tmp/largest_clone_ali_trim.phy_phyml_tree.txt', format='newick')
    leaves = tree.get_terminals()
    for leaf in leaves:
        lid = seqs_unique[int(leaf.name.split('_')[-1]) - 1]
        if len(seqs_red[lid]) == 1:
            leaf.name = seqs_red[lid][0]
        else:
            # Grow a hanging subtree
            for lid in seqs_red[lid]:
                leaf.clades.append(
                        leaf.__class__(
                            branch_length=1e-9,
                            name=lid,
                            )
                        )
    tree.root.name = 'root'

    print('Make graph out of tree')
    noden = 0
    node_ids = ['clonotype_center']
    edges = []
    # Depth first search starting from the root
    # The children make it to the dict at the time of their parent
    for node in tree.find_clades(order='preorder'):
        # Weird base of the tree
        if not hasattr(node, 'clades'):
            continue

        for ic, child in enumerate(node.clades):
            if child.name is not None:
                node_ids.append(child.name)
            else:
                node_ids.append('internal')
            edges.append((noden, noden + ic + 1))
        noden += len(node.clades)
    graph = {
        'nodes': list(range(len(node_ids))),
        'node_ids': node_ids,
        'edges': edges,
        }

    import igraph as ig
    g = ig.Graph(edges)

    print('Calculate graph layout')
    layout_alg = 'fr'
    layout_kwargs = {}
    lo = np.array(g.layout(layout_alg, **layout_kwargs).coords)

    print('Plot pretty graph of clonotype')
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
            1, 2, figsize=(6, 4),
            gridspec_kw={'width_ratios': [4, 1]})
    ax = axs[0]
    sdata = []
    for iid, nid in enumerate(graph['node_ids']):
        sd = {'number': iid}
        if sn+'_'+nid not in df.index:
            sd['s'] = 3
            sd['c'] = 'k'
        else:
            sd['s'] = 15
            sd['c'] = colors_d[df.at[sn+'_'+nid, 'isotype']]

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
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

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

    plt.ion()
    plt.show()
