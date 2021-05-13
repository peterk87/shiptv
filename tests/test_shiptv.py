#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `shiptv` package."""

from os.path import exists
from pathlib import Path

import pandas as pd
from Bio import Phylo
from Bio.Phylo.Newick import Tree
from pandas.testing import assert_frame_equal, assert_series_equal
from typer.testing import CliRunner

from shiptv.cli import app
from shiptv.shiptv import collapse_branches, fix_collection_date

runner = CliRunner()

dirpath = Path(__file__).parent
input_ref_genbank = dirpath / 'data/fmdv-5.gb'
input_newick = dirpath / 'data/fmdv-5.newick'
expected_table = dirpath / 'data/expected_table.tsv'


def test_command_line_interface():
    """Test the CLI."""
    # check that test input files exist
    assert input_ref_genbank.exists()
    assert input_newick.exists()
    assert expected_table.exists()
    # CLI tests start here
    result = runner.invoke(app)
    assert result.exit_code != 0
    assert 'Error: Missing option' in result.output
    help_result = runner.invoke(app, ['--help'])
    assert help_result.exit_code == 0
    assert 'Show this message and exit.' in help_result.output
    with runner.isolated_filesystem():
        out_html = 'test.html'
        out_table = 'test.tsv'
        out_newick = 'test.newick'
        test_result = runner.invoke(app, ['-r', str(input_ref_genbank.absolute()),
                                          '-n', str(input_newick.absolute()),
                                          '-N', out_newick,
                                          '-o', out_html,
                                          '-M', out_table])
        assert test_result.exit_code == 0
        assert exists(out_html)
        assert exists(out_table)
        assert exists(out_newick)
        assert open(input_newick).read() != open(out_newick).read()
        df_out = pd.read_table(out_table)
        df_exp = pd.read_table(expected_table)
        print(f'Col diff: {set(df_out.columns) ^ set(df_exp.columns)}')
        assert_frame_equal(df_exp, df_out)

    with runner.isolated_filesystem():
        out_html = 'test.html'
        out_table = 'test.tsv'
        test_result = runner.invoke(app, ['-r', str(input_ref_genbank.absolute()),
                                          '-n', str(input_newick.absolute()),
                                          '-o', out_html,
                                          '-M', out_table,
                                          '-C', 95])
        assert test_result.exit_code == 0
        assert exists(out_html)
        assert exists(out_table)
        assert_frame_equal(pd.read_csv(expected_table, sep='\t'), pd.read_csv(out_table, sep='\t'))


def test_collapse_branches():
    before_tree_ascii = """
  _________________________________________________________________ MK088171.1
 |
_|_______________ MK071699.1
 |
 |      _______________________________________________ MH845413.2
 |_____|
       |                  _____ MH784405.1
       |_________________|
                         |______ MH784404.1
    """.strip()
    after_tree_ascii = """
  _________________________________________________________________ MK088171.1
 |
 |_______________ MK071699.1
_|
 |_____________________________________________________ MH845413.2
 |
 |                        _____ MH784405.1
 |_______________________|
                         |______ MH784404.1
    """.strip()

    tree: Tree = Phylo.read(input_newick, 'newick')
    from io import StringIO
    sio = StringIO()
    Phylo.draw_ascii(tree, sio)
    pre_collapse_tree = sio.getvalue().strip()
    assert pre_collapse_tree == before_tree_ascii
    collapse_branches(tree, 95)
    sio = StringIO()
    Phylo.draw_ascii(tree, sio)
    post_collapse_tree = sio.getvalue().strip()
    assert post_collapse_tree == after_tree_ascii


def test_fix_collection_date():
    df = pd.DataFrame(dict(collection_date=['1994/1995', '2000', 'not-a-date', '2009/03/31']))
    fix_collection_date(df)
    expected_years = pd.Series([None, 2000, None, 2009], name='collection_year')
    assert_series_equal(df.collection_year, expected_years)
