Introduction
------------

Pypette is a collection of command line utilities and libraries for analyzing biological data. Current functionality includes:
  * Gene expression quantification
  * Copy number analysis
  * Mutation and SNP analysis
  * Chromosomal rearrangement analysis
  * FASTA, VCF and SAM file manipulation
  * Other miscellaneous functionality

Pypette has no dependencies beyond Python and [samtools](https://github.com/samtools/samtools). For performance reasons, it is recommended that the software is run using [PyPy](http://pypy.org/), a high performance implementation of the Python language. You can build PyPy from [source](http://pypy.org/download.html#building-from-source), or download [portable binaries](https://github.com/squeaky-pl/portable-pypy).  

Installation
------------

To install Pypette, download the latest release and extract it to a folder. Then run the Makefile and add the subfolder bin/ to your PATH. See below for an example:

    wget https://github.com/annalam/pypette/archive/pypette-0.6.tar.gz
    tar -xzf pypette-0.6.tar.gz
    cd pypette-0.6
    make
    export PATH=/some/folder/pypette-0.6/bin:$PATH

Examples
--------

#### Chromosomal rearrangements

Identify structural variants based on read pairs and split reads. Split reads must overlap at least 25bp on both sides of the breakpoint:

    breakfast detect -a25 test.bam test
    breakfast detect -a25 control.bam control

Construct a blacklist of false positive regions based on the control sample:

    breakfast blacklist control.sv > blacklist.txt

Only keep structural variants with at least one read pair and three split reads of evidence, or at least 10 split reads (one -r option must be satisfied). Also discard blacklisted structural variants found in the control:

    breakfast filter -r 1-3-0 -r 0-10-0 --blacklist=blacklist.txt test.sv > filtered.sv

Annotate structural variants with information about nearby genes:

    breakfast annotate -b ensembl.bed filtered.sv > annotated.sv
