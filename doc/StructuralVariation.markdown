# Structural variation

## Usage

Identify structural variants based on paired and single end reads:
```
cd genome_alignments
breakfast detect -a25 test.bam test
breakfast detect -a25 control.bam control
```

Construct a blacklist of false positive regions based on a control sample:
```
breakfast blacklist control.sv > blacklist.txt
```

Filter structural variants based on level of evidence:
```
breakfast filter -r 1-3-0 -r 0-10-0 --blacklist=blacklist.txt test.sv > filtered.sv
```

Annotate structural variants with information about nearby genes:
```
breakfast annotate -b ensembl.bed filtered.sv > annotated.sv
```

## Introduction

Breakfast identifies paired end reads where the mates align discordantly. A paired alignment is considered discordant if both mates align to the genome but to separate chromosomes or more than 100 kb apart. Mates with an alignment quality phred value < 20 are discarded from analysis. Next, individual mates that did not initially align to the reference genome are split into short anchors of configurable length (e.g. 25bp). The anchor pairs are then realigned and searched for discordant alignments using the same criteria as with paired end reads. The full sequences corresponding to discordant anchor pairs are compared against the reference genome to identify exact breakpoints and to analyze for sequence homologies. A discordant anchor pair is discarded if the sequence homology between the read and one of the breakpoint flanking sequences was above 70% for the nucleotides matching with the discordant anchor. The exact breakpoint was determined by selecting the breakpoint associated with the lowest amount of nucleotide mismatches.

After identifying discordant pairs from paired end and split reads, the discordant pairs were reoriented so as to always have the pair with the lower chromosome or coordinate first. Discordant pairs were then clustered using a sliding window approach. A cluster of discordant pairs was accepted as a putative structural variant if it contained at least one paired end read and one split read indicating the structural variant.

To filter out false positives, structural variants can also be called in control samples, and all genomic regions within 1 kb of a breakpoint identified in a control sample can be blacklisted.
