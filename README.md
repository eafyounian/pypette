pypette
=======

Pythonic utilities for the analysis of high throughput sequencing data

Pypette is a collection of command line utilities and libraries for analyzing biological data. Current functionality includes:
  * [ExpressionAnalysis Gene expression quantification]
  * Copy number analysis
  * [MutationCalling Mutation analysis]
  * [StructuralVariation Chromosomal rearrangement analysis]
  * FASTA, VCF and SAM file manipulation
  * Other miscellaneous functionality

Pypette has no dependencies beyond Python and [http://samtools.sourceforge.net/ samtools]. For performance reasons, it is recommended that the software is run using [http://pypy.org/ PyPy], a high performance implementation of the Python language.
