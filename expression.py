#!/bin/env pypy

"""
Tools for the analysis of gene expression data.

Usage:
  expression summarize <expr_file>
  expression ensembl to hugo <expr_file> <gtf_file>
  expression normalize rpk <expr_file>
  expression normalize rpkm <expr_file>
  expression normalize mor <expr_file>
  expression exon outliers <expr_file>
  expression visualize splicing <genes> <fastq_prefix> <out_prefix>
  expression sqlite exons <expr_file> <sqlite_file>
  expression sqlite genes <expr_file>

Options:
  -h --help      Show this screen.
"""

from __future__ import print_function
import sys, docopt, re, os, itertools
from collections import defaultdict, namedtuple
from pypette import info, shell, zopen, mkdir, error, read_flat_seq, revcomplement
import numpypy as np







########################
# EXPRESSION SUMMARIZE #
########################

def summarize(exon_expr_path):
	file = zopen(exon_expr_path)
	header = next(file)
	samples = header.rstrip().split('\t')[4:]
	S = len(samples)
	
	features = {}
	for line in file:
		cols = line.rstrip().split('\t')
		if len(cols) < S: continue
		fname = cols[3]
		expr = np.array([float(x) for x in cols[4:]])
		
		f = features.setdefault(cols[3], [np.zeros(S), 0])
		f[0] += expr
		f[1] += int(cols[2]) - int(cols[1])
	
	print('FEATURE\tLENGTH\t' + '\t'.join(samples))
	for name, f in features.iteritems():
		print('%s\t%d\t%s' % (name, f[1], '\t'.join(str(e) for e in f[0])))






##############################
# EXPRESSION ENSEMBL TO HUGO #
##############################

def ensembl_to_hugo(expr_path, gtf_path):
	translations = {}
	for line in zopen(gtf_path):
		m = re.search('gene_id "(.*?)";.*; gene_name "(.*?)"', line)
		if not m: continue
		translations[m.group(1)] = m.group(2)
	
	file = zopen(expr_path)
	header = next(file)
	headers = header[:-1].split('\t')
	feature_col = headers.index('FEATURE')
	
	sys.stdout.write(header)
	for line in file:
		cols = line[:-1].split('\t')
		cols[feature_col] = translations[cols[feature_col]] + ':' + cols[feature_col]
		print('\t'.join(cols))
	






############################
# EXPRESSION NORMALIZE RPK #
############################

def normalize_rpk(expr_path, bed_path):
	
	features = GenomicFeatures.from_bed(bed_path)
	
	# Construct a dictionary of gene lengths.
	feature_len = defaultdict(int)
	for exon in features:
		feature_len[exon[4]] += exon[3] - exon[2] + 1
	
	# Read the gene expression data in.
	expr_file = zopen(expr_path)
	sys.stdout.write(next(expr_file))
	for line in expr_file:
		cols = line[:-1].split('\t')
		flen = feature_len[cols[0]]
		expr = [float(c) / (flen / 1000.0) for c in cols[1:]]
		print('%s\t%s' % (cols[0], '\t'.join(['%.2f' % e for e in expr])))
		









#############################
# EXPRESSION NORMALIZE RPKM #
#############################

def normalize_rpkm(expr_path, bed_path):
	
	if txome_file.lower().endswith('.bed'):
		features = GenomicFeatures.from_bed(txome_file)
	
	# Construct a dictionary of gene lengths.
	feature_len = defaultdict(int)
	for exon in features:
		feature_len[exon[4]] += exon[3] - exon[2] + 1
	
	# Read the gene expression data in.
	data = np.genfromtxt(expr_path, names=True, dtype=None)
	features = data['FEATURE']
	samples = list(data.dtype.names[1:])
	expr = np.zeros((len(features), len(samples)))
	for s in range(len(samples)):
		expr[:, s] = data[samples[s]]
	
	expr /= np.tile(expr.sum(0) / 1000000, (expr.shape[0], 1))
	
	feature_lens = [float(feature_len[f]) for f in features]
	expr /= np.tile(feature_lens, (len(samples), 1)).T / 1000
	
	print('FEATURE\t%s' % '\t'.join(samples))
	str_format = '\t%.3f' * len(samples)
	for f in range(len(features)):
		sys.stdout.write(features[f])
		print(str_format % tuple(expr[f, :].tolist()))









###########################
# EXPRESSION SQLITE EXONS #
###########################

def sqlite_exons(expr_path, sqlite_path):
	import sqlite3
	db = sqlite3.connect(sqlite_path, isolation_level=None)
	db.execute('PRAGMA synchronous = off;')
	db.execute('''CREATE TABLE exons (
		gene TEXT, exon INTEGER, chrom TEXT, start INTEGER, end INTEGER,
		expression TEXT);''')
	db.execute('CREATE TABLE clinical (key TEXT, value TEXT);')
	
	file = zopen(expr_path)
	header = next(file)
	samples = header[:-1].split('\t')[4:]
	S = len(samples)
	
	features = {}
	
	for line in file:
		cols = line[:-1].split('\t')
		fname = cols[3]
		expr = np.array(float(x) for x in line[:-1].split('\t')[4:])
		
		f = features.setdefault(cols[3], [np.zeros(S), 0])
		f[0] += expr
		f[1] += int(cols[2]) - int(cols[1])
	
	print('FEATURE\tLENGTH\t' + '\t'.join(samples))
	for name, f in features.iteritems():
		print('%s\t%d\t%s' % (name, f[1], '\t'.join(str(e) for e in expr)))
	
	db.close()
	
	
	

#################################
# EXPRESSION VISUALIZE SPLICING #
#################################

def visualize_splicing(genes, fastq_prefix, out_prefix):
	genome_path = '/data/csb/organisms/homo_sapiens/hg19_flat'
	bed_path = '/data/csb/organisms/homo_sapiens/ensembl_68/exons.bed'
	genes = genes.replace(' ', '').split(',')
	min_anchor = 15
	read_len = 90
	trim = read_len - min_anchor
	
	chromosomes = read_flat_seq('/data/csb/organisms/homo_sapiens/hg19_flat')
	
	donors = []
	acceptors = []
	exons = []
	for line in zopen(bed_path):
		cols = line[:-1].split('\t')
		if cols[3] in genes:
			chr = cols[0] if cols[0].startswith('chr') else 'chr'+cols[0]
			chr_seq = chromosomes[chr]
			pos = (int(cols[1])+1, int(cols[2]))
			if cols[5] == '+':
				acceptors.append((chr, '+', pos[0],
					chr_seq[pos[0]-1:pos[0]-1+trim]))
				donors.append((chr, '+', pos[1], chr_seq[pos[1]-trim:pos[1]]))
			elif cols[5] == '-':
				acceptors.append((chr, '-', pos[1],
					revcomplement(chr_seq[pos[1]-trim:pos[1]])))
				donors.append((chr, '-', pos[0],
					revcomplement(chr_seq[pos[0]-1:pos[0]-1+trim])))
			exons.append(pos)
				
	# Remove duplicate acceptors and donors.
	acceptors = list(set(acceptors))
	donors = list(set(donors))
	exons = list(set(exons))
	
	# Calculate the contiguous genomic sequence
	chr = acceptors[0][0]
	if any(a[0] != chr for a in acceptors):
		error('Genes must be in the same chromosome!')
	
	genome_window = (min(a[2] for a in acceptors)-2000,
		max(a[2] for a in acceptors)+2000)
	#contig = chromosomes[chr][genome_window[0]:genome_window[1]]
	
	# Calculate junction sequences
	class Junction(object):
		def __init__(self, name, seq):
			self.name = name
			self.sequence = seq
			self.reads = 0
			self.ratio = 0
	
	junctions = defaultdict(list)   # Group junctions by donor
	for left in donors:
		for right in acceptors:
			name = '%d[%s]_%d[%s]' % (left[2], left[1], right[2], right[1])
			junctions[left].append(Junction(name, left[3] + right[3]))
	print('Generated %d junctions.' % (len(donors) * len(acceptors)))
	
	# Build Bowtie index
	index_fasta_path = '%s_ref.fa' % out_prefix
	index = open(index_fasta_path, 'w')
	#index.write('>contig\n%s\n' % contig)
	for donor in junctions:
		for junc in junctions[donor]:
			index.write('>%s\n%s\n' % (junc.name, junc.sequence))
	index.close()
	shell('/data/csb/tools/bowtie-0.12.9/bowtie-build -q %s %s_index' %
		(index_fasta_path, out_prefix))
	
	# Align reads against junctions and tally junction read counts.
	shell('bowtie -v1 -B1 -p8 %s_index <(gunzip -c %s_1.fq.gz %s_2.fq.gz) '
		'> %s.bowtie' % (out_prefix, fastq_prefix, fastq_prefix, out_prefix))
	junction_by_name = {}
	for donor in junctions:
		for j in junctions[donor]: junction_by_name[j.name] = j
	for line in open('%s.bowtie' % out_prefix):
		cols = line[:-1].split('\t')
		if not '_' in cols[2]: continue
		junction_by_name[cols[2]].reads += 1
	
	# Calculate junction power relative to all outgoing links from donor
	for donor in junctions:
		total = sum(j.reads for j in junctions[donor])
		if total <= 0: continue
		for j in junctions[donor]:
			j.ratio = float(j.reads) / total
			if j.reads > 0:
				print('%s: %.1f%% (%d)' % (j.name, j.ratio*100, j.reads))
		
	# Check which exons actually participate in the mature transcripts
	active_edges = []
	for donor in junctions:
		for j in junctions[donor]:
			if j.ratio < 0.05: continue
			active_edges += [int(x[:-3]) for x in j.name.split('_')]
	
	exons = [[ex[0], ex[1], False] for ex in exons]
	ties = []
	for edge in set(active_edges):
		matches = [ex for ex in exons if edge in ex]
		if len(matches) == 1: matches[0][2] = True  # Unique match, mark active
		if len(matches) > 1: ties.append(matches)
	for tie in ties:
		if not any(ex[2] for ex in tie):
			for ex in tie: ex[2] = True   # If still tied, mark all tied active
	
	# Print exon map
	from svgfig import Rect, Frame, Poly
	rects = [Rect(ex[0], 1, ex[1], 2, stroke='none',
		fill='whitesmoke', stroke_linejoin='miter')
		for ex in exons if not ex[2]]
	rects += [Rect(ex[0], 1, ex[1], 2, stroke='none',
		fill='black', stroke_linejoin='miter') for ex in exons if ex[2]]
	lines = []
	for donor in junctions:
		for j in junctions[donor]:
			start, end = [int(x[:-3]) for x in j.name.split('_')]
			lines.append(Poly([(start,2), ((start+end)/2,3), (end,2)],
				stroke_opacity=j.ratio))
		
	Frame(genome_window[0], genome_window[1], 0, 10, *(rects+lines),
		width=500).SVG().save('%s.svg' % out_prefix)
	
	shell('rm %s_index.* %s_ref.fa' % (out_prefix, out_prefix))



#######################
# COMMAND LINE PARSER #
#######################
	
if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	
	if args['summarize']:
		summarize(args['<expr_file>'])
	elif args['ensembl'] and args['to'] and args['hugo']:
		ensembl_to_hugo(args['<expr_file>'], args['<gtf_file>'])
	elif args['normalize'] and args['rpk']:
		normalize_rpk(args['<expr_file>'])
	elif args['normalize'] and args['rpkm']:
		normalize_rpkm(args['<expression_file>'])
	elif args['normalize'] and args['mor']:
		normalize_mor(args['<expr_file>'])
	elif args['visualize'] and args['splicing']:
		visualize_splicing(args['<genes>'], args['<fastq_prefix>'],
			args['<out_prefix>'])
	elif args['exon'] and args['outliers']:
		exon_outliers(args['<expr_file>'])
	elif args['sqlite'] and args['exons']:
		sqlite_exons(args['<expr_file>'], args['<sqlite_file>'])
	

