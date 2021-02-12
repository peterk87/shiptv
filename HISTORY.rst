=======
History
=======

0.4.0 (2021-02-12)
------------------

* Add Select2 metadata field select box to show metadata fields in the order required
* Add Shift+Ctrl+F to toggle tree viz full window mode
* Using BioPython's Phylo module for phylogenetic tree manipulation. Removed ete3.
* Using Requests library to get JS and CSS to embed into HTML output
* Migrated to Typer from Click for CLI
* Added Rich for nicer logging and tracebacks
* Changed some options to be optional; only required options are ``-n/--newick`` and ``-o/--output-html``
* By default don't highlight user samples as a field of metadata in tree viz
* Allow user to optionally specify outgroup taxa name
* Allow user to optionally re-root tree at midpoint node
* Fixed warnings from Pandas when modifying view of dataframe by using ``.loc[]``
* Moved to GitHub Actions CI and PyPI deployment


0.3.0 (2019-10-02)
------------------

* Fix rendering of numeric and ISO date metadata fields
* Reference genomes Genbank now optional
* Updated docs


0.2.0 (2019-06-28)
------------------

* Added low support branch highlighting in tree HTML file
* Added option to collapse low support branches (`-C/--collapse-support`)
* Added option to output modified Newick tree file (`-N/--output-newick`)
* Fixed date/time parsing from Genbank files (#1)

0.1.0 (2019-05-10)
------------------

* First release on PyPI.
