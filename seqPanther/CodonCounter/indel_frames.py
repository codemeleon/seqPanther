#!/usr/bin/env python

import pandas as pd

from Bio.Seq import Seq


def indel_frames(indel_pos_type_size, params):
    gff_data = params["gff_data"]
    alt_codon_frac = params["alt_codon_frac"]
    rid = params["rid"]
    sample = params["sample"]
    indels_changes = indel_pos_type_size.copy()
    coors = set(indel_pos_type_size["coor"])
    indel_pos_type_size_fragmented = []
    indel_pos_type_size["r_shift"] = 0

    for coor in coors:
        # TODO: use df.to_dict('records'). more detail https://towardsdatascience.com/heres-the-most-efficient-way-to-iterate-through-your-pandas-dataframe-4dad88ac92ee
        t_gff_data = gff_data[(gff_data["start"] <= coor)
                              & (gff_data["end"] >= coor)].to_dict('records')

        if not len(t_gff_data):
            print(f'No CDS gff data found for {coor}')
            continue
        for gff_row in t_gff_data:
            t_indel_pos_type_size = indel_pos_type_size[
                indel_pos_type_size["coor"] == coor]

            shift = (coor - gff_row["start"] +
                     1) % 3  # original reference is zero based

            t_indel_pos_type_size["shift"] = shift
            # t_indel_pos_type_size.loc[t_indel_pos_type_size.indel < 0,
            # "shift"] = t_indel_pos_type_size.loc[
            # t_indel_pos_type_size.indel < 0,
            # "shift"].apply(lambda x: (x + 2) % 3)

            # t_indel_pos_type_size.loc[(t_indel_pos_type_size["coor"] == coor) &
            # (t_indel_pos_type_size["indel"] < 0),
            # "shift"] += 1
            t_indel_pos_type_size.loc[(
                t_indel_pos_type_size["coor"] == coor
            ), "r_shift"] = (3 - t_indel_pos_type_size.loc[
                (t_indel_pos_type_size["coor"] == coor), "shift"].values) % 3
            t_indel_pos_type_size.loc[
                t_indel_pos_type_size["coor"] == coor,
                ["shift", "r_shift"]] = t_indel_pos_type_size.loc[
                    t_indel_pos_type_size["coor"] == coor,
                    ["shift", "r_shift"]].applymap(lambda x: 3 - x)

            amino_pos = (coor - gff_row["start"]) // 3 + 1

            if gff_row["strand"] == "-":
                amino_len = (gff_row["end"] - gff_row["start"] + 1) // 3
                amino_pos = amino_len - amino_pos - 1

            t_indel_pos_type_size.loc[t_indel_pos_type_size["coor"] == coor,
                                      "amino_pos"] = amino_pos
            t_indel_pos_type_size.loc[
                t_indel_pos_type_size["coor"] == coor,
                "amino_pos"] = t_indel_pos_type_size[
                    t_indel_pos_type_size["coor"] == coor].apply(
                        lambda x: x["amino_pos"] + (1
                                                    if x['indel'] > 0 else 2),
                        axis=1)
            t_indel_pos_type_size.loc[
                t_indel_pos_type_size["coor"] == coor,
                "codon_pos"] = coor - t_indel_pos_type_size.loc[
                    t_indel_pos_type_size["coor"] == coor, "shift"]
            t_indel_pos_type_size.loc[
                t_indel_pos_type_size["coor"] == coor,
                'gene_range'] = (f"{gff_row['start']+1}:{gff_row['end']+1}"
                                 f"({gff_row['strand']})")
            indel_pos_type_size_fragmented.append(t_indel_pos_type_size)
    indel_pos_type_size = pd.concat(indel_pos_type_size_fragmented)
    del indel_pos_type_size_fragmented
    cols = indel_pos_type_size.columns.tolist()
    cols.remove('count')
    indel_pos_type_size = indel_pos_type_size.groupby(
        cols)['count'].sum().reset_index()

    indel_pos_type_size_ref_support = indel_pos_type_size.groupby([
        'coor', 'depth', 'indel', 'ref', 'amino_pos', 'codon_pos', 'gene_range'
    ])['count'].sum().reset_index().rename(columns={'count': 'ref_count'})
    indel_pos_type_size = indel_pos_type_size.merge(
        indel_pos_type_size_ref_support,
        on=[
            'coor', 'depth', 'indel', 'ref', 'amino_pos', 'codon_pos',
            'gene_range'
        ])
    del indel_pos_type_size_ref_support
    indel_pos_type_size['ref_count'] = indel_pos_type_size[
        'depth'] - indel_pos_type_size['ref_count']
    indel_pos_type_size = indel_pos_type_size[
        indel_pos_type_size['count'] > indel_pos_type_size['depth'] *
        alt_codon_frac]

    if indel_pos_type_size.empty:
        return pd.DataFrame(), pd.DataFrame()  # TODO: Add dummy columns

    # indels_changes["coor"] = indels_changes['codon_pos'] + shift
    # indels_changes.apply(
    # lambda x: 2 if len(x['ref']) > len(x['read']) else 1, axis=1)

    indels_changes["Nucleotide Percent"] = indels_changes.apply(
        lambda x: x['count'] / x['depth'], axis=1)

    indels_changes["Nucleotide Frequency"] = indels_changes.apply(
        lambda x: 'del' + x['ref']
        if (len(x['ref']) > len(x['read'])) else 'ins' + x['read'],
        axis=1)
    indel_pos_type_size.loc[indel_pos_type_size.indel > 0,
                            "amino_pos"] = indel_pos_type_size[
                                indel_pos_type_size.indel > 0].apply(
                                    lambda x: x['amino_pos'] -
                                    (len(x['ref']) // 3 - 1),
                                    axis=1)
    indel_pos_type_size.loc[
        (indel_pos_type_size.indel < 0)
        & indel_pos_type_size.gene_range.str.contains('-'),
        "amino_pos"] = indel_pos_type_size[
            (indel_pos_type_size.indel < 0)
            & indel_pos_type_size.gene_range.str.contains('-')].apply(
                lambda x: x['amino_pos'] - (len(x['ref']) // 3 - 2), axis=1)
    indel_pos_type_size.loc[
        (indel_pos_type_size.indel < 0)
        & ~indel_pos_type_size.gene_range.str.contains('-'),
        "amino_pos"] = indel_pos_type_size[
            (indel_pos_type_size.indel < 0)
            & ~indel_pos_type_size.gene_range.str.contains('-')].apply(
                lambda x: x['amino_pos'] - 1, axis=1)
    indel_pos_type_size["Nucleotide Change"] = indel_pos_type_size.apply(
        lambda x: f"{x['coor']+(1 if x['indel'] < 0 else 0)+1}:"
        f"{x['ref'][3:-3] }"
        f">{x['read'][3:-3]}",
        axis=1)
    indel_pos_type_size["ref"] = indel_pos_type_size.apply(
        lambda x: x['ref'][x['shift']:-x['r_shift']], axis=1)
    indel_pos_type_size["read"] = indel_pos_type_size.apply(
        lambda x: x['read'][x['shift']:-x['r_shift']], axis=1)

    indel_pos_type_size.loc[
        indel_pos_type_size.gene_range.str.contains('-'),
        ["ref", "read"]] = indel_pos_type_size.loc[
            indel_pos_type_size.gene_range.str.contains('-'),
            ["ref", "read"]].applymap(
                lambda x: str(Seq(x).reverse_complement()))
    indel_pos_type_size["Codon Change"] = indel_pos_type_size.apply(
        lambda x: f"{x['codon_pos']+1}:{x['ref']}>{x['read']}", axis=1)

    indel_pos_type_size["Amino Acid Change"] = indel_pos_type_size.apply(
        lambda x:
        f"{Seq(x['ref']).translate()}{x['amino_pos']}{Seq(x['read']).translate()}",
        axis=1)
    indel_pos_type_size["codon_count"] = indel_pos_type_size.apply(
        lambda x: f"{x['ref']}-{x['ref_count']};"
        f"{x['read']}-{x['count']}",
        axis=1)
    indel_pos_type_size["codon_percent"] = indel_pos_type_size.apply(
        lambda x: f"{x['ref']}-"
        f"{'%0.2f' % (x['ref_count']*100/x['depth'])};"
        f"{x['read']}-"
        f"{'%0.2f' % (x['count']*100/x['depth']) }",
        axis=1)
    indel_nuc = indel_pos_type_size[[
        'coor',
        'indel',
        'depth',
        'ref',
        'read',
        'count',
    ]]
    indel_nuc[['ref',
               'read']] = indel_nuc[['ref',
                                     'read']].applymap(lambda x: x[3:-3])
    indel_nuc[['ref',
               'read']] = indel_nuc[['ref',
                                     'read']].applymap(lambda x: ''.join(x))
    indel_nuc[['ref_len',
               'read_len']] = indel_nuc[['ref',
                                         'read']].applymap(lambda x: len(x))
    indel_nuc['Nucleotide Frequency'] = indel_nuc.apply(
        lambda x: 'del' + x['ref']
        if x['ref_len'] > x['read_len'] else 'ins' + x['read'],
        axis=1)
    indel_nuc['Nucleotide Percent'] = indel_nuc.apply(
        lambda x: '%s:%0.2f' %
        (x["Nucleotide Frequency"], x['count'] * 100 / x['depth']),
        axis=1)
    indel_nuc['Nucleotide Frequency'] = indel_nuc.apply(
        lambda x: '%s:%d' % (x["Nucleotide Frequency"], x['count']), axis=1)
    indel_nuc = indel_nuc.groupby([
        'coor', 'depth', "indel"
    ]).apply(lambda x:
             [list(x['Nucleotide Frequency']),
              list(x['Nucleotide Percent'])]).reset_index()
    indel_nuc["Nucleotide Frequency"] = indel_nuc.apply(
        lambda x: ','.join(x[0][0]) + f',read_count:{x["depth"]}', axis=1)

    indel_nuc["Nucleotide Percent"] = indel_nuc.apply(
        lambda x: ','.join(x[0][0]), axis=1)
    indel_nuc["coor"] = indel_nuc.apply(lambda x: (x['coor'] + 1) if
                                        (x['indel'] > 0) else (x['coor'] + 2),
                                        axis=1).values
    indel_nuc = indel_nuc.drop(columns=['depth', 0, 'indel'])

    indel_pos_type_size = indel_pos_type_size.drop([
        "coor", "ref_count", "amino_pos", "codon_pos", "read", "ref", "indel",
        "count", "depth"
    ],
                                                   axis=1)
    indel_pos_type_size.insert(0, 'Reference ID', rid)
    indel_pos_type_size.insert(0, 'Sample', sample)
    indel_pos_type_size = indel_pos_type_size.drop(
        columns=['shift', 'r_shift'])
    indel_nuc.insert(0, 'Reference ID', rid)
    indel_nuc.insert(0, 'Sample', sample)
    return indel_pos_type_size, indel_nuc  # TODO: retrun only one table
