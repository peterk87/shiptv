=================================================
Standalone HTML Interactive Phylogenetic Tree Viz
=================================================


.. image:: https://img.shields.io/pypi/v/shiptv.svg
        :target: https://pypi.python.org/pypi/shiptv

.. image:: https://img.shields.io/travis/peterk87/shiptv.svg
        :target: https://travis-ci.org/peterk87/shiptv

.. image:: https://readthedocs.org/projects/shiptv/badge/?version=latest
        :target: https://shiptv.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Generate a standalone HTML file with an interactive phylogenetic tree using PhyloCanvas


* Free software: Apache Software License 2.0
* Documentation: https://shiptv.readthedocs.io.

*See test shiptv HTML output here:*

- `fmdv-5-shiptv.html`_


**Phylogenetic tree of 5 FMDV genomes**

.. image:: docs/images/fmdv5.png
        :alt: Phylogenetic tree of 5 FMDV genomes

**Phylogenetic tree of IAV HA gene sequences**

.. image:: docs/images/iav-ha-gene-tree-scrn.png
        :alt: Phylogenetic tree of IAV HA gene sequences


Features
--------

* Interactively view your tree in the browser with metadata highlighted beside your tree using PhyloCanvas.
* Automatically retrieve metadata from a GenBank file!
* Visualize your own metadata! Provide a tab-delimited table as input with ``--user-sample-metadata /path/to/your-table.tsv``
* Highlight branches with low support in the browser.
* Collapse branches with low support (e.g. ``-C 95`` for IQ-TREE trees with UFBoot ``-bb 1000`` to collapse branches with less than 95% support). 

Usage
-----

Show help

.. code-block:: bash

    shiptv --help


Help output

.. code-block:: 

    Usage: shiptv [OPTIONS]

      Create HTML tree visualization with metadata.

      The metadata for reference genomes can be extracted from the specified
      Genbank file.

      Any leaf names that are present in the tree but not present in the Genbank
      file are assumed to be user samples and are flagged as such in the
      metadata table as "user_sample"="Yes".

    Options:
      -r, --ref-genomes-genbank FILE  Reference genome sequences Genbank file
      -n, --newick FILE               Phylogenetic tree Newick file  [required]
      -N, --output-newick PATH        Output Newick file
      -o, --output-html PATH          Output HTML tree path  [required]
      -m, --output-metadata-table PATH
                                      Output metadata table path  [required]
      --leaflist PATH                 Optional leaf names to select from
                                      phylogenetic tree for pruned tree
                                      visualization. One leaf name per line.
      --genbank-metadata-fields PATH  Optional fields to extract from Genbank
                                      source metadata. One field per line.
      --user-sample-metadata PATH     Optional tab-delimited metadata for user
                                      samples to join with metadata derived from
                                      reference genome sequences Genbank file.
                                      Sample IDs must be in the first column.
      --metadata-fields-in-order PATH
                                      Optional list of fields in order to output
                                      in metadata table and HTML tree
                                      visualization. One field per line.
      --dont-fix-metadata             Do not automatically fix metadata (only on
                                      Genbank file metadata)
      -C, --collapse-support FLOAT    Collapse internal branches below specified
                                      bootstrap support value (default -1 for no
                                      collapsing)
      --help                          Show this message and exit.


With a reference sequence Genbank file `ref.gb`, a Newick format phylogenetic tree `tree.nwk`, output a `tree.html` standalone HTML interactive phylogenetic tree visualization and a `metadata-table.tsv` tab-delimited table of metadata from `ref.gb`.

.. code-block:: bash

    shiptv -r ref.gb -n tree.nwk -o tree.html -m metadata-table.tsv



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _`fmdv-5-shiptv.html`: docs/data/fmdv-5-shiptv.html
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
