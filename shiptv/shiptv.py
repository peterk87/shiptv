# -*- coding: utf-8 -*-

"""Main module."""
import logging
from pathlib import Path
from typing import List, Dict, Optional

import jinja2
import pandas as pd
import requests
from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade
from Bio.SeqRecord import SeqRecord
from pkg_resources import resource_filename

html_template = resource_filename('shiptv', 'tmpl/phylocanvas.html')
phylocanvas_metadata_plugin_html = resource_filename('shiptv', 'tmpl/phylocanvas_metadata_plugin.html')
resources = {
    'bootstrap_css': 'https://stackpath.bootstrapcdn.com/bootstrap/4.2.1/css/bootstrap.min.css',
    'jquery_js': 'https://code.jquery.com/jquery-3.3.1.slim.min.js',
    'popper_js': 'https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.6/umd/popper.min.js',
    'bootstrap_js': 'https://stackpath.bootstrapcdn.com/bootstrap/4.2.1/js/bootstrap.min.js',
    'phylocanvas_js': 'https://unpkg.com/phylocanvas@2.8.1/dist/phylocanvas.js',
    'polyfill_js': 'https://unpkg.com/phylocanvas@2.8.1/polyfill.js',
    'chroma_js': 'https://unpkg.com/chroma-js@2.0.2/chroma.js',
    'lodash_js': 'https://unpkg.com/lodash@4.17.15/lodash.min.js',
    'ag_grid_js': 'https://unpkg.com/ag-grid-community/dist/ag-grid-community.min.noStyle.js',
    'ag_grid_css': 'https://unpkg.com/ag-grid-community/dist/styles/ag-grid.css',
    'ag_grid_theme_css': 'https://unpkg.com/ag-grid-community/dist/styles/ag-theme-balham.css',
    'select2_js': 'https://cdnjs.cloudflare.com/ajax/libs/select2/3.5.4/select2.min.js',
    'select2_css': 'https://cdnjs.cloudflare.com/ajax/libs/select2/3.5.4/select2.min.css',
}


def read_lines(filepath: Path) -> List[str]:
    out = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            out.append(line)
    return out


def genbank_source_metadata(rec: SeqRecord) -> Dict[str, str]:
    """Get source feature metadata dictionary for a SeqRecord"""
    return {k: v[0] if v is not None and len(v) == 1 else v
            for k, v in rec.features[0].qualifiers.items()}


def genbank_metadata(genbank: Path) -> pd.DataFrame:
    """Parse genome metadata from Genbank file into Pandas DataFrame.
    """
    id_to_rec = {r.id: r for r in SeqIO.parse(genbank, 'genbank')}
    df_metadata = pd.DataFrame({gid: genbank_source_metadata(rec)
                                for gid, rec in id_to_rec.items()}).transpose()
    if 'isolate' in df_metadata and 'strain' in df_metadata:
        df_metadata['strain'] = df_metadata['isolate'] \
            .combine_first(df_metadata['strain'])
    return df_metadata


def fix_host_metadata(df: pd.DataFrame) -> None:
    cattle_syn = '''
    Bos taurus
    bovine
    cattle (Ankole cow, sentinel herd)
    Ankole cow
    Cattle
    '''.strip().split('\n')
    df.loc[df.host.isin(cattle_syn), 'host'] = 'cattle'
    sheep_syn = '''
    ovine
    '''.strip().split('\n')
    df.loc[df.host.isin(sheep_syn), 'host'] = 'sheep'
    pig_syn = '''
    sus scrofa domesticus
    swine
    porcine'''.strip().split('\n')
    df.loc[df.host.isin(pig_syn), 'host'] = 'pig'


def fix_collection_date(df_metadata: pd.DataFrame):
    dates = pd.to_datetime(df_metadata.collection_date, errors='coerce')
    years = [x.year if not pd.isnull(x) else None for x in dates]
    df_metadata['collection_year'] = years
    df_metadata.collection_date = [str(x).split()[0] if not pd.isnull(x) else None for x in dates]


def fix_country_region(df: pd.DataFrame):
    df['region'] = df.country.str.extract(r'.*:\s*(.*)\s*')
    df['country'] = df.country.str.extract(r'([^:]+)(:\s*.*\s*)?')[0]


def add_user_metadata(df: pd.DataFrame, user_sample_metadata: Path) -> None:
    df_user_metadata = pd.read_csv(user_sample_metadata, sep='\t', index_col=0, dtype='str')
    logging.info(f'type of "{df_user_metadata.index.values[0]}" is {type(df_user_metadata.index.values[0])} '
                 f'| dtype={df_user_metadata.index.dtype}')
    logging.info(f'Read user sample metadata table from '
                 f'"{user_sample_metadata}" with '
                 f'{df_user_metadata.shape[0]} rows and columns: {list(df_user_metadata.columns)}')
    for user_column in df_user_metadata.columns:
        if user_column not in df.columns:
            df[user_column] = None
    for idx, row in df_user_metadata.iterrows():
        if idx not in df.index:
            logging.info(f'idx "{idx}" not in df.index!')
            continue
        original_row = df.loc[idx, :]
        row_dict = original_row[~pd.isnull(original_row)].to_dict()
        row_dict.update(row.to_dict())
        df.loc[idx, :] = pd.Series(row_dict)


def parse_tree(newick: Path, outgroup=None, midpoint_root=False, ladderize=False) -> Tree:
    # Read phylogenetic tree newick file using ete3
    tree: Tree = Phylo.read(newick, 'newick')
    # setting user specified outgroup takes precedence over setting midpoint node as outgroup/tree root
    if outgroup:

        tree.root_with_outgroup(dict(name=outgroup))
    elif midpoint_root:
        tree.root_at_midpoint()
    # ladderize the tree for subjectively cleaner viz
    if ladderize:
        tree.ladderize()
    return tree


def write_html_tree(df_metadata: pd.DataFrame,
                    output_html: Path,
                    tree: Tree) -> None:
    with open(html_template) as fh, \
         open(phylocanvas_metadata_plugin_html) as fh_pmd, \
         open(output_html, 'w') as fout:
        tmpl = jinja2.Template(fh.read())
        logging.info('Retrieving JS and CSS for HTML tree')
        scripts_css = {}
        for k, v in resources.items():
            logging.info(f'Getting HTML resource "{k}" from "{v}"')
            scripts_css[k] = requests.get(v).text
        logging.info('Retrieved JS and CSS for HTML tree')
        from io import StringIO
        sio = StringIO()
        Phylo.write(tree, sio, 'newick')
        fout.write(tmpl.render(newick_string=sio.getvalue().strip(),
                               metadata_json_string=df_metadata.to_json(orient='index'),
                               phylocanvas_metadata_plugin=fh_pmd.read(),
                               **scripts_css))


def parse_leaf_list(leaflist_path: Path) -> Optional[List[str]]:
    leaflist = None
    if leaflist_path:
        leaflist = read_lines(leaflist_path)
        logging.info(f'Read {len(leaflist)} leaf names from "{leaflist_path}"')
    return leaflist


def prune_tree(df_metadata: pd.DataFrame, leaflist: Optional[List[str]], tree: Tree):
    if leaflist:
        n_nodes_before_prune = len(tree.get_terminals())
        leaves_to_keep = set(leaflist)
        for node in tree.get_terminals():
            if node.name not in leaves_to_keep:
                tree.prune(node)
        logging.info(f'Pruned tree to {len(tree.get_terminals())} leaves from {n_nodes_before_prune} leaves.')
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
                            'region',
                            'collection_date',
                            'collection_year',
                            'host',
                            'isolation_source']
    return metadata_columns


def add_user_samples_field(df_metadata: pd.DataFrame, metadata_columns: List[str],
                           tree_leaf_names: List[str]) -> pd.DataFrame:
    set_user_samples = set(tree_leaf_names) - set(df_metadata.index)
    logging.info(f'Found {len(set_user_samples)} user samples. {";".join(set_user_samples)}')
    df_metadata = subset_metadata_table(df_metadata, metadata_columns, tree_leaf_names)
    df_metadata['user_sample'] = None
    df_metadata.loc[set_user_samples, 'user_sample'] = 'Yes'
    logging.info(f'Highlighting {len(set_user_samples)} user samples with '
                 f'new field "user_sample" where user\'s samples have '
                 f'value "Yes"')
    return df_metadata


def subset_metadata_table(df_metadata, metadata_columns, tree_leaf_names):
    df_metadata = df_metadata.reindex(index=tree_leaf_names)
    df_metadata = df_metadata.loc[:, metadata_columns]
    return df_metadata


def try_fix_serotype_metadata(df_metadata: pd.DataFrame) -> None:
    if 'serotype' in df_metadata.columns:
        df_metadata.serotype = df_metadata.serotype.str.replace(r'\s', '')
        logging.info('Stripped internal whitespace from serotype metadata field')


def try_fix_country_metadata(df_metadata: pd.DataFrame) -> None:
    if 'country' in df_metadata.columns:
        fix_country_region(df_metadata)
        logging.info('Split country metadata field into country and region on ":"')


def try_fix_collection_date_metadata(df_metadata: pd.DataFrame) -> None:
    if 'collection_date' in df_metadata.columns:
        fix_collection_date(df_metadata)
        logging.info('Normalized collection_date metadata field to ISO date '
                     'format and extracted year integer into '
                     'collection_year field')


def try_fix_host_metadata(df_metadata: pd.DataFrame) -> None:
    if 'host' in df_metadata.columns:
        before_host_types = df_metadata.host.unique().size
        host_before = set(df_metadata.host.unique())
        fix_host_metadata(df_metadata)
        after_host_types = df_metadata.host.unique().size
        host_after = set(df_metadata.host.unique())
        logging.info(f'Harmonized some host info for cattle, sheep and pig '
                     f'to use more consistent host type categories. '
                     f'Before: {before_host_types} host types. '
                     f'After: {after_host_types} host types. Difference: {host_after ^ host_before}')


def collapse_branches(tree: Tree, collapse_support: float) -> None:
    """Collapse internal branches below support threshold

    Note:
        This function modifies the supplied `tree` object.

    Args:
        tree: Bio.Phylo Tree object
        collapse_support: Support threshold
    """
    node: Clade
    for node in tree.get_nonterminals():
        if node.confidence and node.confidence < collapse_support:
            tree.collapse(node)
