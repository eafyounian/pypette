#!/bin/env pypy

"""
Tools for the manipulation of GTF files.

Usage:
  gtf cleanup ensembl <gtf_file>
  gtf to gene bed <gtf_file>
  gtf to transcript bed <gtf_file>
  gtf to exon bed <gtf_file>
  gtf to composite exon bed <gtf_file>

Options:
  -h --help         Show this screen.
"""

from __future__ import print_function
import sys, docopt, re, os
from pypette import zopen, shell, revcomplement, info, error, shell_stdout
from sam import read_sam
import numpypy as np




###################
# GTF TO GENE BED #
###################

def gtf_to_gene_bed(gtf_path):
	gene_id_to_name = {}
	gene_exons = {}
	
	gtf_file = zopen(gtf_path)
	for line in gtf_file:
		tokens = line.split('\t')
		if not tokens[2] == 'exon': continue
		if tokens[1] in ['nonsense_mediated_decay', 'retained_intron']: continue
		
		chr = tokens[0]
		if not chr.startswith('chr'):
			chr = 'chr' + chr
			
		pos = (int(tokens[3]), int(tokens[4]))
		strand = tokens[6]
		
		m = re.search('gene_id "(.+?)"', line)
		gene_id = m.group(1)
		
		m = re.search('gene_name "(.+?)"', line)
		gene_name = m.group(1)
		
		exons = gene_exons.setdefault(gene_id, [])
		exons.append((chr, strand, pos[0], pos[1]))
		
		gene_id_to_name[gene_id] = gene_name
		
	for gene_id, exons in gene_exons.iteritems():
		if not all([exon[0] == exons[0][0] for exon in exons]):
			info('Chromosome confusion detected.')
		if not all([exon[1] == exons[0][1] for exon in exons]):
			info('Strand confusion detected.')
		
		pos = (min([exon[2] for exon in exons]),
			max([exon[3] for exon in exons]))
		
		print('%s\t%d\t%d\t%s (%s)\t\t%s' % (exons[0][0], pos[0] - 1, pos[1],
			gene_id_to_name[gene_id], gene_id, exons[0][1]))
	
	gtf_file.close()
	







#########################
# GTF TO TRANSCRIPT BED #
#########################

def gtf_to_transcript_bed(gtf_path):
	tx_id_to_name = {}
	tx_exons = {}
	
	gtf_file = open(gtf_path)
	for line in gtf_file:
		tokens = line.split('\t')
		if not tokens[2] == 'exon': continue
		if tokens[1] in ['nonsense_mediated_decay', 'retained_intron']: continue
		
		chr = tokens[0]
		if not chr.startswith('chr'):
			chr = 'chr' + chr
			
		pos = (int(tokens[3]), int(tokens[4]))
		strand = tokens[6]
		
		m = re.search('transcript_id "(.+?)"', line)
		tx_id = m.group(1)
		
		m = re.search('gene_name "(.+?)"', line)
		gene_name = m.group(1)
		
		exons = tx_exons.setdefault(tx_id, [])
		exons.append((chr, strand, pos[0], pos[1]))
		
		tx_id_to_name[tx_id] = gene_name
		
	for tx_id, exons in tx_exons.iteritems():
		#if not all([exon[0] == exons[0][0] for exon in exons]):
		#if not all([exon[1] == exons[0][1] for exon in exons]):
	#		print('Strand confusion detected.')
		
		pos = (min([exon[2] for exon in exons]),
			max([exon[3] for exon in exons]))
		
		print('%s\t%d\t%d\t%s (%s)\t\t%s' % (exons[0][0], pos[0] - 1, pos[1],
			tx_id_to_name[tx_id], tx_id, exons[0][1]))
	
	gtf_file.close()








###################
# GTF TO EXON BED #
###################

def gtf_to_exon_bed(gtf_path):
	for line in zopen(gtf_path):
		c = line.split('\t')
		if not c[2] == 'exon': continue
		if c[1] in ['nonsense_mediated_decay', 'retained_intron']: continue
		if c[0].startswith('HSCHR'): continue
		if c[0].startswith('GL'): continue
		if c[0].startswith('HG'): continue
		
		gene = re.search('gene_id "(.+?)"', line).group(1)
		chr = 'chr'+c[0] if not c[0].startswith('chr') else c[0]
		print('%s\t%d\t%s\t%s\t\t%s' % (chr, int(c[3])-1, c[4], gene, c[6]))









#############################
# GTF TO COMPOSITE EXON BED #
#############################

def gtf_to_composite_exon_bed(bed_path):
	
	genes = {}
	
	bed_file = open(bed_path)
	for line in bed_file:
		chr, start, end, gene_id = line[:-1].split('\t')[:4]
		start = int(start)
		end = int(end)
		
		gene = genes.setdefault(gene_id, [chr, [(start, end)]])
		if chr != gene[0]:
			print('ERROR: Chromosome mismatch.')
		
		exons = gene[1]
		overlapping = [ex for ex in exons if end >= ex[0] and start <= ex[1]]
		disjoint = [ex for ex in exons if not (end >= ex[0] and start <= ex[1])]
		
		disjoint.append((min([start] + [ex[0] for ex in overlapping]),
			max([end] + [ex[1] for ex in overlapping])))
		gene[1] = disjoint
		
	for gene_id, gene in genes.iteritems():
		exons = gene[1]
		for ex in exons:
			print('%s\t%d\t%d\t%s' % (gene[0], ex[0], ex[1], gene_id))




###################
# CLEANUP ENSEMBL #
###################

def cleanup_ensembl(gtf_path):
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
	if args['to'] and args['gene'] and args['bed']:
		gtf_to_gene_bed(args['<gtf_file>'])
	elif args['to'] and args['transcript'] and args['bed']:
		gtf_to_transcript_bed(args['<gtf_file>'])
	elif args['to'] and args['exon'] and args['bed']:
		gtf_to_exon_bed(args['<gtf_file>'])
	elif args['to'] and args['composite'] and args['exon'] and args['bed']:
		gtf_to_composite_exon_bed(args['<gtf_file>'])
	elif args['cleanup'] and args['ensembl']:
		cleanup_ensembl(args['<gtf_file>'])



