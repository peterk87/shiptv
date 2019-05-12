# -*- coding: utf-8 -*-

"""Main module."""
import jinja2
import logging
import pandas as pd
import pkg_resources
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from typing import List, Dict, Optional

html_template = pkg_resources.resource_filename('shiptv', 'tmpl/phylocanvas.html')


def read_lines(filepath: str) -> List[str]:
    out = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            out.append(line)
    return out


def genbank_source_metadata(rec: SeqRecord) -> Dict[str, str]:
    """Get source feature metadata dictionary for a SeqRecord"""
    return {k: v[0] if v is not None and len(v) == 1 else v for k, v in rec.features[0].qualifiers.items()}


def genbank_metadata(genbank: str) -> pd.DataFrame:
    """Parse genome metadata from Genbank file into Pandas DataFrame.
    """
    id_to_rec = {r.id: r for r in SeqIO.parse(genbank, 'genbank')}
    df_metadata = pd.DataFrame({gid: genbank_source_metadata(rec) for gid, rec in id_to_rec.items()}).transpose()
    if 'isolate' in df_metadata and 'strain' in df_metadata:
        df_metadata['strain'] = df_metadata['isolate'].combine_first(df_metadata['strain'])
    return df_metadata


def fix_host_metadata(df_metadata: pd.DataFrame) -> None:
    cattle_syn = '''
    Bos taurus
    bovine
    cattle (Ankole cow, sentinel herd)
    Ankole cow
    Cattle
    '''.strip().split('\n')
    df_metadata.host[df_metadata.host.isin(cattle_syn)] = 'cattle'
    sheep_syn = '''
    ovine
    '''.strip().split('\n')
    df_metadata.host[df_metadata.host.isin(sheep_syn)] = 'sheep'
    pig_syn = '''
    sus scrofa domesticus
    swine
    porcine'''.strip().split('\n')
    df_metadata.host[df_metadata.host.isin(pig_syn)] = 'pig'


def fix_collection_date(df_metadata):
    years = [x.year if not pd.isnull(x) else x for x in pd.to_datetime(df_metadata.collection_date)]
    df_metadata.collection_date = pd.to_datetime(df_metadata.collection_date)
    df_metadata['collection_year'] = years
    df_metadata.collection_date = [str(x).split()[0] if not pd.isnull(x) else None for x in df_metadata.collection_date]


def fix_country_region(df_metadata):
    df_metadata['region'] = df_metadata.country.str.extract(r'.*:\s*(.*)\s*')
    df_metadata['country'] = df_metadata.country.str.extract(r'([^:]+)(:\s*.*\s*)?')[0]


def add_user_metadata(df_metadata, user_sample_metadata):
    df_user_metadata = pd.read_csv(user_sample_metadata, sep='\t', index_col=0)
    logging.info(f'Read user sample metadata table from '
                 f'"{user_sample_metadata}" with '
                 f'{df_user_metadata.shape[0]} rows and columns: '
                 f'{";".join(df_user_metadata.columns)}')
    for user_column in df_user_metadata.columns:
        if user_column not in df_metadata.columns:
            df_metadata[user_column] = None
    for idx, row in df_user_metadata.iterrows():
        if idx not in df_metadata.index:
            continue
        original_row = df_metadata.loc[idx, :]
        row_dict = original_row[~pd.isnull(original_row)].to_dict()
        row_dict.update(row.to_dict())
        df_metadata.loc[idx, :] = pd.Series(row_dict)


def parse_tree(newick):
    # Read phylogenetic tree newick file using ete3
    tree = Tree(newick=newick)
    # Calculate the midpoint node
    midpoint_node = tree.get_midpoint_outgroup()
    # Set midpoint as output
    tree.set_outgroup(midpoint_node)
    return tree


def write_html_tree(df_metadata, output_html, tree):
    with open(html_template) as fh, open(output_html, 'w') as fout:
        tmpl = jinja2.Template(fh.read())
        fout.write(tmpl.render(newick_string=tree.write(),
                               metadata_json_string=df_metadata.to_json(orient='index')))


def parse_leaf_list(leaflist: str) -> Optional[List[str]]:
    if leaflist:
        leaflist_filepath = leaflist
        leaflist = read_lines(leaflist)
        logging.info(f'Read {len(leaflist)} leaf names from "{leaflist_filepath}"')
    return leaflist


def prune_tree(df_metadata, leaflist, tree):
    if leaflist is not None and len(leaflist) > 0:
        n_nodes_before_prune = len(tree)
        tree.prune(leaflist)
        logging.info(f'Pruned tree to {len(tree)} leaves from {n_nodes_before_prune} leaves.')
        df_metadata = df_metadata.loc[leaflist, :]
    return df_metadata


def reorder_metadata_fields(df_metadata, metadata_fields_in_order):
    if metadata_fields_in_order:
        df_metadata = df_metadata[read_lines(metadata_fields_in_order)]
    return df_metadata


def get_metadata_fields(genbank_metadata_fields):
    if genbank_metadata_fields:
        metadata_columns = read_lines(genbank_metadata_fields)
    else:
        metadata_columns = ['strain',
                            'serotype',
                            'note',
                            'country',
                            'collection_date',
                            'host',
                            'isolation_source']
    return metadata_columns


def highlight_user_samples(df_metadata: pd.DataFrame, metadata_columns: List[str], tree_leaf_names: List[str]) -> pd.DataFrame:
    set_user_samples = set(tree_leaf_names) - set(df_metadata.index)
    logging.info(f'Found {len(set_user_samples)} user samples. {";".join(set_user_samples)}')
    df_metadata = df_metadata.reindex(index=tree_leaf_names)
    df_metadata = df_metadata.loc[:, metadata_columns]
    df_metadata['user_sample'] = None
    df_metadata.loc[set_user_samples, 'user_sample'] = 'Yes'
    logging.info(f'Highlighting {len(set_user_samples)} user samples with '
                 f'new field "user_sample" where user\'s samples have '
                 f'value "Yes"')
    return df_metadata


def try_fix_serotype_metadata(df_metadata: pd.DataFrame) -> None:
    if 'serotype' in df_metadata.columns:
        df_metadata.serotype = df_metadata.serotype.str.replace(r'\s', '')
        logging.info(f'Stripped internal whitespace from serotype metadata field')


def try_fix_country_metadata(df_metadata: pd.DataFrame) -> None:
    if 'country' in df_metadata.columns:
        fix_country_region(df_metadata)
        logging.info(f'Split country metadata field into country and region on ":"')


def try_fix_collection_date_metadata(df_metadata: pd.DataFrame) -> None:
    if 'collection_date' in df_metadata.columns:
        fix_collection_date(df_metadata)
        logging.info(f'Normalized collection_date metadata field to ISO date '
                     f'format and extracted year integer into '
                     f'collection_year field')


def try_fix_host_metadata(df_metadata: pd.DataFrame) -> None:
    if 'host' in df_metadata.columns:
        before_host_types = df_metadata.host.unique().size
        fix_host_metadata(df_metadata)
        after_host_types = df_metadata.host.unique().size
        logging.info(f'Harmonized some host info for cattle, sheep and pig '
                     f'to use more consistent host type categories. '
                     f'Before: {before_host_types} host types. '
                     f'After: {after_host_types} host types.')
