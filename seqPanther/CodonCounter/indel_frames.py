#!/usr/bin/env python

import pandas as pd

from Bio import Seq


def indel_frames(indel_pos_type_size, params):
    gff_data = params["gff_data"]
    alt_codon_frac = params["alt_codon_frac"]
    rid = params["rid"]
    sample = params["sample"]
    coors = set(indel_pos_type_size["coor"])
    indel_pos_type_size["amino_pos"] = 0
    indel_pos_type_size["codon_pos"] = 0
    indel_pos_type_size["shift"] = 0
    indel_pos_type_size["r_shift"] = 0
    indel_pos_type_size['gene_range'] = '-'
    shift, r_shift = 0, 0
    indel_pos_type_size_fragmented = []

    for coor in coors:
        # TODO: use df.to_dict('records'). more detail https://towardsdatascience.com/heres-the-most-efficient-way-to-iterate-through-your-pandas-dataframe-4dad88ac92ee
        t_gff_data = gff_data[
            (gff_data["start"] <= coor)
            & (gff_data["end"] > coor)].to_dict(
                'record'
            )  # TODO: Need to make some correction when indel contain end site of CDS or stop codon
        t_indel_pos_type_size = indel_pos_type_size[indel_pos_type_size["coor"]
                                                    == coor]
        if not len(t_gff_data):
            print(f'No CDS gff data found for {coor}')
            continue
        for gff_row in t_gff_data:

            adjusted_coor = coor + 1
            shift = (adjusted_coor - gff_row["start"]) % 3  # + 1
            r_shift = (3 - shift) % 3
            t_indel_pos_type_size.loc[
                t_indel_pos_type_size["coor"] == coor,
                ["ref", "read"]] = t_indel_pos_type_size.loc[
                    t_indel_pos_type_size["coor"] == coor,
                    ["ref", "read"]].applymap(
                        lambda x: x[3 - shift:-(3 - r_shift)])
            t_indel_pos_type_size.loc[t_indel_pos_type_size["coor"] == coor,
                                      ["shift", "r_shift"]] = [shift, r_shift]

            amino_pos = (adjusted_coor -
                         gff_row["start"]) // 3  # - 1 * row["indel"]
            amino_pos = amino_pos + 1 if shift else amino_pos
            codon_pos = adjusted_coor - shift

            # TODO: Get the correct sequence

            if gff_row["strand"] == "-":
                t_indel_pos_type_size.loc[
                    t_indel_pos_type_size["coor"] == coor,
                    ["ref", "read"]] = t_indel_pos_type_size.loc[
                        t_indel_pos_type_size["coor"] == coor,
                        ["ref", "read"]].applymap(
                            lambda x: Seq.Seq(x).reverse_complement())
                amino_len = (gff_row["end"] - gff_row["start"]) // 3
                amino_pos = amino_len - amino_pos + 1  # TODO: Change in case of deletion
                amino_pos = amino_pos + 1 if shift else amino_pos
                # in case of insertions
                # TODO: Fix codon Position in reverse direction

                shift, r_shift = r_shift, shift
            t_indel_pos_type_size.loc[t_indel_pos_type_size["coor"] == coor,
                                      "amino_pos"] = amino_pos
            t_indel_pos_type_size.loc[t_indel_pos_type_size["coor"] == coor,
                                      "codon_pos"] = codon_pos
            t_indel_pos_type_size.loc[
                t_indel_pos_type_size["coor"] == coor,
                'gene_range'] = f"{gff_row['start']}:{gff_row['end']}({gff_row['strand']})"
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

    indels_changes = indel_pos_type_size.copy()
    indels_changes["coor"] = indels_changes['codon_pos'] + shift
    # indels_changes.apply(
    # lambda x: 2 if len(x['ref']) > len(x['read']) else 1, axis=1)

    indels_changes["Nucleotide Percent"] = indels_changes.apply(
        lambda x: x['count'] / x['depth'], axis=1)

    indels_changes["Nucleotide Frequency"] = indels_changes.apply(
        lambda x: 'del' + x['ref'][r_shift:-shift] if (len(x['ref']) > len(x[
            'read'])) else 'ins' + x['read'][r_shift:-shift],
        axis=1)

    indel_pos_type_size["Amino Acid Change"] = indel_pos_type_size.apply(
        lambda x:
        f"{Seq.Seq(x['ref']).translate()}{x['amino_pos']}{Seq.Seq(x['read']).translate()}",
        axis=1)
    indel_pos_type_size["Nucleotide Change"] = indel_pos_type_size.apply(
        lambda x:
        f"{x['codon_pos']+x['shift']+(1 if x['indel'] < 0 else 0)}:{x['ref'][x['shift']:-x['r_shift'] if x['r_shift']  else None]}>{x['read'][x['shift']:-x['r_shift'] if x['r_shift']  else None]}",
        axis=1)
    indel_pos_type_size["Codon Change"] = indel_pos_type_size.apply(
        lambda x: f"{x['codon_pos']}:{x['ref']}>{x['read']}", axis=1)

    indel_pos_type_size["codon_count"] = indel_pos_type_size.apply(
        lambda x: f"{x['ref']}-{x['ref_count']};{x['read']}-{x['count']}",
        axis=1)
    indel_pos_type_size["codon_percent"] = indel_pos_type_size.apply(
        lambda x:
        f"{x['ref']}-{'%0.2f' % (x['ref_count']*100/x['depth'])};{x['read']}-{'%0.2f' % (x['count']*100/x['depth']) }",
        axis=1)
    # TODO: Make Change to nucleide chan showing
    indel_nuc = indel_pos_type_size[['coor', 'depth', 'ref', 'read', 'count']]
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
        'coor', 'depth'
    ]).apply(lambda x:
             [list(x['Nucleotide Frequency']),
              list(x['Nucleotide Percent'])]).reset_index()
    indel_nuc["Nucleotide Frequency"] = indel_nuc.apply(
        lambda x: ','.join(x[0][0]) + f',read_count:{x["depth"]}', axis=1)

    indel_nuc["Nucleotide Percent"] = indel_nuc.apply(
        lambda x: ','.join(x[0][0]), axis=1)
    indel_nuc = indel_nuc.drop(columns=['depth', 0])

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
