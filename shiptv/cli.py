# -*- coding: utf-8 -*-

"""Console script for shiptv."""
import logging
from pathlib import Path
from sys import version_info
from typing import Optional

from Bio import Phylo
import pandas as pd
import typer
from rich.logging import RichHandler

from shiptv.shiptv import (
    genbank_metadata,
    add_user_metadata,
    parse_tree,
    write_html_tree,
    parse_leaf_list,
    prune_tree,
    reorder_metadata_fields,
    get_metadata_fields,
    add_user_samples_field,
    try_fix_serotype_metadata,
    try_fix_country_metadata,
    try_fix_collection_date_metadata,
    try_fix_host_metadata,
    collapse_branches,
    subset_metadata_table,
)
from shiptv import __version__

app = typer.Typer()


def version_callback(value: bool):
    if value:
        typer.echo(f"shiptv version {__version__}")
        raise typer.Exit()


@app.command(
    epilog=f"shiptv version {__version__}; Python {version_info.major}.{version_info.minor}.{version_info.micro}"
)
def main(
    newick: Path = typer.Option(
        ..., "-n", "--newick", help="Phylogenetic tree Newick file"
    ),
    output_html: Path = typer.Option(
        "shiptv-tree.html", "-o", "--output-html", help="Output HTML tree path"
    ),
    output_newick: Optional[Path] = typer.Option(
        None, "-N", "--output-newick", help="Output Newick file"
    ),
    ref_genomes_genbank: Optional[Path] = typer.Option(
        None,
        "-r",
        "--ref-genomes-genbank",
        help="Reference genome sequences Genbank file",
    ),
    output_metadata_table: Optional[Path] = typer.Option(
        None, "-M", "--output-metadata-table", help="Output metadata table path"
    ),
    leaflist: Path = typer.Option(
        None,
        help="Optional leaf names to select from phylogenetic tree "
        "for pruned tree visualization. One leaf name per line.",
    ),
    genbank_metadata_fields: Path = typer.Option(
        None,
        help="Optional fields to extract from Genbank source "
        "metadata. One field per line.",
    ),
    user_sample_metadata: Path = typer.Option(
        None,
        "-m",
        "--metadata",
        help="Optional tab-delimited metadata for user samples to join with"
        " metadata derived from reference genome sequences Genbank "
        "file. Sample IDs must be in the first column.",
    ),
    metadata_fields_in_order: Path = typer.Option(
        None,
        help="Optional list of fields in order to output in metadata "
        "table and HTML tree visualization. One field per line.",
    ),
    fix_metadata: bool = typer.Option(
        True, help="Try to automatically fix metadata from reference Genbank file."
    ),
    collapse_support: float = typer.Option(
        -1.0,
        "-C",
        "--collapse-support",
        help="Collapse internal branches below specified bootstrap "
        "support value (default -1 for no collapsing)",
    ),
    highlight_user_samples: bool = typer.Option(
        False, help="Highlight user samples with metadata field in tree."
    ),
    outgroup: str = typer.Option(None, help="Tree outgroup taxa"),
    midpoint_root: bool = typer.Option(False, help="Set midpoint root"),
    verbose: bool = typer.Option(False, help="Verbose logs"),
    version: Optional[bool] = typer.Option(
        None,
        callback=version_callback,
        help=f'Print "shiptv version {__version__}" and exit',
    ),
):
    """shiptv - create an HTML tree visualization with optional metadata highlighting

    Typical usage:

    $ shiptv --newick tree.nwk --metadata tree-metadata.tsv --output-html tree.html

    Note: The metadata for reference genomes can be extracted from the specified Genbank file.
    """
    from rich.traceback import install

    install(show_locals=True)

    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.INFO if not verbose else logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )

    tree = parse_tree(newick, outgroup=outgroup, midpoint_root=midpoint_root)
    leaf_names = [x.name for x in tree.get_terminals()]
    if collapse_support != -1:
        logging.info(
            f"Collapsing internal branches with support values less "
            f"than {collapse_support}"
        )
        collapse_branches(tree, collapse_support)
    if output_newick:
        Phylo.write([tree], output_newick, "newick")
    if ref_genomes_genbank:
        df_metadata = genbank_metadata(ref_genomes_genbank)
        logging.info(
            f'Parsed metadata from "{ref_genomes_genbank}" with columns: {list(df_metadata.columns)}'
        )
        if fix_metadata:
            try_fix_serotype_metadata(df_metadata)
            try_fix_country_metadata(df_metadata)
            try_fix_collection_date_metadata(df_metadata)
            try_fix_host_metadata(df_metadata)
        else:
            logging.warning("Not fixing any genome metadata.")

        metadata_fields = get_metadata_fields(genbank_metadata_fields)
        # only use columns present in the reference genome metadata
        metadata_fields = [x for x in metadata_fields if x in list(df_metadata.columns)]
        logging.info(f"Metadata table fields: {metadata_fields}")

    else:
        df_metadata = pd.DataFrame(index=leaf_names)
        logging.info(
            f"df_metadata index dtype = {df_metadata.index.dtype} "
            f"| {df_metadata.index.values[0]} is {type(df_metadata.index.values[0])}"
        )
        metadata_fields = []

    if highlight_user_samples:
        df_metadata = add_user_samples_field(df_metadata, metadata_fields, leaf_names)
    else:
        df_metadata = subset_metadata_table(df_metadata, metadata_fields, leaf_names)
    if leaflist:
        df_metadata = prune_tree(df_metadata, parse_leaf_list(leaflist), tree)
    if user_sample_metadata:
        add_user_metadata(df_metadata, user_sample_metadata)
    df_metadata = reorder_metadata_fields(df_metadata, metadata_fields_in_order)
    if output_metadata_table:
        df_metadata.to_csv(output_metadata_table, sep="\t")
        logging.info(
            f"Wrote tab-delimited genome metadata table with "
            f"{df_metadata.shape[0]} rows and "
            f"{df_metadata.shape[1]} columns to "
            f'"{output_metadata_table}"'
        )
    write_html_tree(df_metadata, output_html, tree)
    logging.info(f'Wrote HTML tree to "{output_html}"')


if __name__ == "__main__":
    app()  # pragma: no cover
