pypette
=======

Pythonic utilities for the analysis of high throughput sequencing data

Pypette is a collection of command line utilities and libraries for analyzing biological data. Current functionality includes:
  * Gene expression quantification
  * Copy number analysis
  * Mutation analysis
  * [Chromosomal rearrangement analysis](doc/StructuralVariation.markdown)
  * FASTA, VCF and SAM file manipulation
  * Other miscellaneous functionality

Pypette has no dependencies beyond Python and [samtools](http://samtools.sourceforge.net/). For performance reasons, it is recommended that the software is run using [PyPy](http://pypy.org/), a high performance implementation of the Python language.
