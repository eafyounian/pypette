Pypette is a collection of command line utilities and libraries for analyzing biological data. Current functionality includes:
  * Gene expression quantification
  * Copy number analysis
  * Mutation and SNP analysis
  * [Chromosomal rearrangement analysis](doc/StructuralVariation.markdown)
  * FASTA, VCF and SAM file manipulation
  * Other miscellaneous functionality

Pypette has no dependencies beyond Python and [samtools](http://samtools.sourceforge.net/). For performance reasons, it is recommended that the software is run using [PyPy](http://pypy.org/), a high performance implementation of the Python language. You can build PyPy from [source](http://pypy.org/download.html#building-from-source), or download [portable binaries](https://github.com/squeaky-pl/portable-pypy).  

Installation
------------

To install Pypette, download the latest release and extract it to a folder. Then run the Makefile and add the subfolder bin/ to your PATH. See below for an example:

    wget https://github.com/annalam/pypette/archive/pypette-0.6.tar.gz
    tar -xzf pypette-0.6.tar.gz
    cd pypette-0.6
    make
    export PATH=/some/folder/pypette-0.6/bin:$PATH
