#!/bin/env python

"""
Tools for creating Betastasis-compatible files.

Usage:
  betastasis gene expression <expr_file> <sqlite_file>
  betastasis small rna expression <expr_file> <sqlite_file>
  betastasis convert gene expression <sqlite_file>

Options:
  -h --help      Show this screen.
"""

from __future__ import print_function
import sys, docopt, re, os, sqlite3, json
from pypette import info, shell, zopen, mkdir, error


class Archive:
	def __init__(self, path, mode='r'):
		self.mode = mode
		self.db = sqlite3.connect(path, isolation_level=None)
		self.db.execute('PRAGMA synchronous = off;')
		if mode == 'w':
			self.db.execute('CREATE TABLE entries (key TEXT, value TEXT);')
			self.db.execute('BEGIN')

	def add(self, key, value):
		self.db.execute('INSERT INTO entries VALUES (?, ?)',
			(key, json.dumps(value)))

	def close(self):
		if self.mode == 'w': self.db.execute('COMMIT')
		self.db.close()



##############################
# BETASTASIS GENE EXPRESSION #
##############################

def gene_expression(expr_path, sqlite_path):
	file = zopen(expr_path)
	header = next(file)
	samples = header.rstrip().split('\t')[1:]
	genes = []

	archive = Archive(sqlite_path, 'w')
	for n, line in enumerate(file):
		cols = line.rstrip().split('\t')
		gene = cols[0]
		archive.add('gene/' + gene, {
#			'chromosome': cols[0], 'start': int(cols[1]), 'end': int(cols[2]),
			'expression': [float(c) for c in cols[1:]]
		})
		genes.append(gene)
	archive.add('clinical', { 'sample_id': samples })
	archive.add('index', { 'genes': sorted(genes) })
	archive.close()





###################################
# BETASTASIS SMALL RNA EXPRESSION #
###################################

def small_rna_expression(expr_path, sqlite_path):
	file = zopen(expr_path)
	header = next(file)
	samples = header[:-1].split('\t')[1:]
	samples = [re.sub('.bam$', '', s, re.I) for s in samples]
	
	names = []

	archive = Archive(sqlite_path, 'w')
	for n, line in enumerate(file):
		cols = line[:-1].split('\t')
		archive.add('gene/' + cols[0], {
			'expression': [float(c) for c in cols[1:]]
		})
		names.append(cols[0])
	archive.add('clinical', { 'sample_id': samples })
	archive.add('index', { 'genes': sorted(names) })
	archive.close()







######################################
# BETASTASIS CONVERT GENE EXPRESSION #
######################################

def convert_gene_expression(sqlite_path):
	archive = Archive(sqlite_path, 'w')
	genes = []
	for dirpath, _, files in os.walk('expr'):
		for path in files:
			if not path.endswith('.json'): continue
			if path == 'features.json': continue
			file = open(dirpath + '/' + path)
			d = json.load(file)
			file.close()
			archive.add('gene/' + path.replace('.json', ''), {
				'expression': [2**x for x in d['data']]
			})
			genes.append(path.replace('.json', ''))
	clinical = json.load(open('clinical.json'))
	archive.add('clinical', clinical)
	archive.add('index', { 'genes': sorted(genes) })
	archive.close()





	
	


#######################
# COMMAND LINE PARSER #
#######################
	
if __name__ == '__main__':
	args = docopt.docopt(__doc__)

	if args['convert'] and args['gene'] and args['expression']:
		convert_gene_expression(args['<sqlite_file>'])
	elif args['gene'] and args['expression']:
		gene_expression(args['<expr_file>'], args['<sqlite_file>'])
	elif args['small'] and args['rna'] and args['expression']:
		small_rna_expression(args['<expr_file>'], args['<sqlite_file>'])

	

