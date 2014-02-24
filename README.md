Introduction
------------

Pypette is a collection of command line utilities and libraries for analyzing biological data. Current functionality includes:
  * Gene expression quantification
  * Copy number analysis
  * Mutation and SNP analysis
  * Chromosomal rearrangement analysis
  * FASTA, VCF and SAM file manipulation
  * Other miscellaneous functionality

For performance reasons, it is recommended that the software is run using [PyPy](http://pypy.org/), a high performance implementation of the Python language. You can build PyPy from [source](http://pypy.org/download.html#building-from-source), or download [portable binaries](https://github.com/squeaky-pl/portable-pypy).  

Installation
------------

To install Pypette, download the latest release and extract it to a folder. Then run the Makefile and add the subfolder bin/ to your PATH. See below for an example:

    wget https://github.com/annalam/pypette/archive/0.7.tar.gz
    tar -xzf pypette-0.7.tar.gz
    cd pypette-0.7
    make
    export PATH=/some/folder/pypette-0.7/bin:$PATH

Some Pypette functionality requires external software to be installed:
- [samtools](https://github.com/samtools/samtools): required for mutation and structural variant calling
- [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml): required for structural variant calling
- [ANNOVAR](http://www.openbioinformatics.org/annovar/): required for mutation annotation

Examples
--------

#### Chromosomal rearrangements

Create a flattened genome FASTA file for use in split read analysis:

    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
    fasta flatten human_g1k_v37.fasta.gz hg19_flat

Identify structural variants based on read pairs and split reads. Split reads must overlap at least 25bp on both sides of the breakpoint:

    breakfast detect -a25 test.bam hg19_flat test
    breakfast detect -a25 control.bam hg19_flat control

Construct a blacklist of false positive regions based on the control sample:

    breakfast blacklist control.sv > blacklist.txt

Only keep structural variants with at least one read pair and three split reads of evidence, or at least 10 split reads (one -r option must be satisfied). Also discard blacklisted structural variants found in the control:

    breakfast filter -r 1-3-0 -r 0-10-0 --blacklist=blacklist.txt test.sv > filtered.sv

Annotate structural variants with information about nearby genes:

    wget ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.gz
    gtf to gene bed Homo_sapiens.GRCh37.74.gtf.gz > ensembl_genes.bed
    breakfast annotate -b ensembl_genes.bed filtered.sv > annotated.sv
