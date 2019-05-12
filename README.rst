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


Features
--------

* TODO

Usage
-----

Show help

.. code-block:: bash

    shiptv --help


Help output

.. code-block:: 

    Usage: shiptv [OPTIONS]

      Create HTML tree visualization with metadata.

      The metadata for reference genomes is extracted from the specified Genbank
      file.

      Any leaf names that are present in the tree but not present in the Genbank
      file are assumed to be user samples and are flagged as such in the
      metadata table as "user_sample"="Yes".

    Options:
      -r, --ref-genomes-genbank TEXT  Reference genome sequences Genbank file
                                      [required]
      -n, --newick TEXT               Phylogenetic tree Newick file  [required]
      -o, --output-html TEXT          Output HTML tree path  [required]
      -m, --output-metadata-table TEXT
                                      Output metadata table path  [required]
      --leaflist TEXT                 Optional leaf names to select from
                                      phylogenetic tree for pruned tree
                                      visualization. One leaf name per line.
      --genbank-metadata-fields TEXT  Optional fields to extract from Genbank
                                      source metadata. One field per line.
      --user-sample-metadata TEXT     Optional tab-delimited metadata for user
                                      samples to join with metadata derived from
                                      reference genome sequences Genbank file.
                                      Sample IDs must be in the first column.
      --metadata-fields-in-order TEXT
                                      Optional list of fields in order to output
                                      in metadata table and HTML tree
                                      visualization. One field per line.
      --dont-fix-metadata             Do not automatically fix metadata
      --help                          Show this message and exit.


With a reference sequence Genbank file `ref.gb`, a Newick format phylogenetic tree `tree.nwk`, output a `tree.html` standalone HTML interactive phylogenetic tree visualization and a `metadata-table.tsv` tab-delimited table of metadata from `ref.gb`.

.. code-block:: bash

    shiptv -r ref.gb -n tree.nwk -o tree.html -m metadata-table.tsv



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
