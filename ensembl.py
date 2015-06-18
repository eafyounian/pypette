#!/bin/env pypy

"""
Tools for analyzing Ensembl annotations.

Usage:
  ensembl gene bed <gtf_file>
  ensembl transcript bed <gtf_file>
  ensembl exon bed <gtf_file>
  ensembl cds bed <gtf_file>
  ensembl cleanup <gtf_file>

Options:
  -h --help         Show this screen.
"""

from __future__ import print_function
import sys, docopt, re, os
from pypette import zopen, shell, revcomplement, info, error, shell_stdout
from sam import read_sam

human_chr = ['%d' % c for c in range(1, 23)] + ['X', 'Y']
human_chr = set(human_chr + ['chr%s' % c for c in human_chr])
reasonable_gene_types = ['protein_coding', 'processed_transcript',
	'pseudogene', 'lincRNA']


####################
# ENSEMBL GENE BED #
####################

def ensembl_gene_bed(gtf_path):
	gene_id_to_name = {}
	gene_exons = {}
	
	gtf_file = zopen(gtf_path)
	for line in gtf_file:
		if line.startswith('#'): continue
		c = line.rstrip('\n').split('\t')
		if not c[0] in human_chr: continue
		if not c[1] in reasonable_gene_types: continue
		if c[2] != 'exon': continue
		
		chr, start, end, strand = c[0], int(c[3]), int(c[4]), c[6]
		if not chr.startswith('chr'): chr = 'chr' + chr
		
		m = re.search(r'gene_id "(.+?)"', line)
		gene_id = m.group(1)
		
		m = re.search(r'gene_name "(.+?)"', line)
		gene_name = m.group(1)
		
		exons = gene_exons.setdefault(gene_id, [])
		exons.append((chr, strand, pos[0], pos[1]))
		
		gene_id_to_name[gene_id] = gene_name
		
	for gene_id, exons in gene_exons.iteritems():
		if not all(exon[0] == exons[0][0] for exon in exons):
			info('Chromosome confusion detected.')
		if not all(exon[1] == exons[0][1] for exon in exons):
			info('Strand confusion detected.')
		
		start, end = min(ex[2] for ex in exons), max(ex[3] for ex in exons)
		print('%s\t%d\t%d\t%s (%s)\t\t%s' % (exons[0][0], start - 1, end,
			gene_id_to_name[gene_id], gene_id, exons[0][1]))
	
	







##########################
# ENSEMBL TRANSCRIPT BED #
##########################

def ensembl_transcript_bed(gtf_path):
	tx_id_to_gene = {}
	tx_exons = {}
	
	gtf_file = zopen(gtf_path)
	for line in gtf_file:
		if line.startswith('#'): continue
		c = line.rstrip('\n').split('\t')
		if not c[0] in human_chr: continue
		if not c[1] in reasonable_gene_types: continue
		if c[2] != 'exon': continue
		
		chr, start, end, strand = c[0], int(c[3]), int(c[4]), c[6]
		if not chr.startswith('chr'): chr = 'chr' + chr
		
		tx_id = re.search(r'transcript_id "(.+?)"', line).group(1)		
		gene_id = re.search(r'gene_id "(.+?)"', line).group(1)
		gene_name = re.search(r'gene_name "(.+?)"', line).group(1)
		
		exons = tx_exons.setdefault(tx_id, [])
		exons.append((chr, strand, start, end))
		tx_id_to_gene[tx_id] = (gene_id, gene_name)
		
	for tx_id, exons in tx_exons.iteritems():
		start, end = min(ex[2] for ex in exons), max(ex[3] for ex in exons)
		print('%s\t%d\t%d\t%s:%s:%s\t\t%s' % (exons[0][0], start - 1, end,
			tx_id_to_gene[tx_id][0], tx_id_to_gene[tx_id][1], tx_id,
			exons[0][1]))
	








####################
# ENSEMBL EXON BED #
####################

def ensembl_exon_bed(gtf_path):
	for line in zopen(gtf_path):
		if line.startswith('#'): continue
		c = line.rstrip('\n').split('\t')
		if not c[0] in human_chr: continue
		if not c[1] in reasonable_gene_types: continue
		if c[2] != 'exon': continue

		chr, start, end = c[0], int(c[3]), int(c[4])
		if not chr.startswith('chr'): chr = 'chr' + chr
		gene_id = '%s:%s' % (re.search(r'gene_id "(.+?)"', line).group(1),
			re.search(r'gene_name "(.+?)"', line).group(1))
		
		print('%s\t%d\t%s\t%s\t\t%s' % (chr, start-1, end, gene_id, c[6]))










####################
# ENSEMBL CDS BED #
####################

def ensembl_cds_bed(gtf_path):
	for line in zopen(gtf_path):
		if line.startswith('#'): continue
		c = line.rstrip('\n').split('\t')
		if not c[0] in human_chr: continue
		if not c[1] in reasonable_gene_types: continue
		if c[2] != 'CDS': continue

		chr, start, end = c[0], int(c[3]), int(c[4])
		if not chr.startswith('chr'): chr = 'chr' + chr
		gene_id = '%s:%s' % (re.search(r'gene_id "(.+?)"', line).group(1),
			re.search(r'gene_name "(.+?)"', line).group(1))
		
		print('%s\t%d\t%s\t%s\t\t%s' % (chr, start-1, end, gene_id, c[6]))







###################
# ENSEMBL CLEANUP #
###################

def ensembl_cleanup(gtf_path):
	valid_chr = set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
		'11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
		'21', '22', 'X', 'Y', 'MT'])
	for line in zopen(gtf_path):
		chr = line[:line.find('\t')]
		if chr in valid_chr and not 'nonsense_mediated' in line:
			sys.stdout.write('chr')
			sys.stdout.write(line)





#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['gene'] and args['bed']:
		ensembl_gene_bed(args['<gtf_file>'])
	elif args['transcript'] and args['bed']:
		ensembl_transcript_bed(args['<gtf_file>'])
	elif args['exon'] and args['bed']:
		ensembl_exon_bed(args['<gtf_file>'])
	elif args['cds'] and args['bed']:
		ensembl_cds_bed(args['<gtf_file>'])
	elif args['cleanup']:
		ensembl_cleanup(args['<gtf_file>'])



