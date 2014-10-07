#!/bin/env pypy

"""
Toolkit for the manipulation of FASTA and FASTQ sequence data.

Usage:
  fasta filter <fasta> <regexp>
  fasta flatten <fasta> <out_dir>
  fasta split <fasta> <tag_length> <anchors_5p_path> <anchors_3p_path>
  fasta split interleaved <fasta> <tag_length>
  fasta trim <fasta> <trim_length>
  fasta find adapters <fasta> <seed>
  fasta remove adapters <fasta> <adapter>
  fasta color2nuc <fasta>
  fasta rna2dna <fasta>
  fasta dna2rna <fasta>
  fasta check <fasta_1> <fasta_2> <out_fasta_1> <out_fasta_2>
  fasta repair <fasta_1> <fasta_2> <out_fasta_1> <out_fasta_2>
  fasta from sra <sra_file>

Options:
  -h --help      Show this screen.
"""

from __future__ import print_function
import sys, subprocess, docopt, re, os
from collections import defaultdict
from pypette import zopen, shell, read_fasta, revcomplement, GenomicFeatures
from pypette import info, error









################
# FASTA FILTER #
################

def fasta_filter(fasta_path, regexp):
	fasta = open(fasta_path)
	keep = False
	for line in fasta:
		if line[0] in ('>', '@'):
			keep = re.search(regexp, line) != None
		if keep:
			sys.stdout.write(line)





#################
# FASTA FLATTEN #
#################

def fasta_flatten(fasta_path, output_dir):
	fasta = zopen(fasta_path)
	flat_file = None
	for line in fasta:
		if line[0] in ('>'):
			if flat_file: flat_file.close()
			flat_file = open(output_dir + '/' + line[1:].strip() + '.seq', 'w')
		else:
			flat_file.write(line.strip())
	flat_file.close()









###############
# FASTA SPLIT #
###############

def fasta_split(fasta_path, tag_length, anchors_5p_path, anchors_3p_path):
	fasta = zopen(fasta_path, 'r')
	tags_5p = zopen(anchors_5p_path, 'w')
	tags_3p = zopen(anchors_3p_path, 'w')
	
	R = 0
	
	for line in fasta:
		if line[0] == '+':
			next(fasta)
			continue
		if line[0] in '>@#': continue
		
		R += 1
		if len(line) <= 2 * tag_length: continue
		tags_5p.write('>%d_%s/1\n%s\n' % (R, line[:-1], line[0:tag_length]))
		tags_3p.write('>%d_%s/2\n%s\n' % (R, line[:-1], line[-tag_length-1:-1]))
	
	tags_5p.close()
	tags_3p.close()






###########################
# FASTA SPLIT INTERLEAVED #
###########################

def fasta_split_interleaved(fasta_path, tag_length):
	fasta = zopen(fasta_path, 'r')
	
	R = 0
	for line in fasta:
		if line[0] == '+':
			next(fasta)
			continue
		if line[0] in '>@#': continue
		
		R += 1
		if len(line) <= 2 * tag_length: continue
		print('>%d_%s/1\n%s' % (R, line[:-1], line[0:tag_length]))
		print('>%d_%s/2\n%s' % (R, line[:-1], line[-tag_length-1:-1]))








##############
# FASTA TRIM #
##############

def fasta_trim(fasta_path, trim_len):
	fasta = zopen(fasta_path)
	for line in fasta:
		if line[0] in '@+>':
			sys.stdout.write(line)
			line = next(fasta)
			if len(line) - 1 > trim_len:
				print(line[:trim_len])
			else:
				sys.stdout.write(line)








#######################
# FASTA FIND ADAPTERS #
#######################

def fasta_find_adapters(fasta_path, seed):
	suffixes = []
	seed_len = len(seed)
	
	fasta = open(fasta_path)
	for line in fasta:
		if line[0] not in '>@': continue
		seq = next(fasta)[:-1]
		pos = seq.find(seed)
		if pos != -1:
			suffixes.append(seq[pos+seed_len:])
	
	suffix_count = Counter(suffixes)
	for suffix, count in suffix_count.most_common(20):
		print('%s\t%d' % (suffix, count))









#########################
# FASTA REMOVE ADAPTERS #
#########################

def fasta_remove_adapters(fasta_path, adapter):
	# Convert the adapter into a regular expression
	if len(adapter) < 5: error('Adapter sequence is too short.')
	adapter_re = adapter[:5]
	for base in adapter[5:]:
		adapter_re += '(?:' + base
	adapter_re += (len(adapter) - 5) * ')?'
	adapter_re = re.compile(adapter_re)

	info('Adapter regular expression: %s' % adapter_re)

	fasta = zopen(fasta_path)
	for line in fasta:
		if line[0] == '#':
			sys.stdout.write(line)
		elif line[0] == '>':
			sys.stdout.write(line)
			seq = next(fasta)[:-1]
			m = adapter_re.search(seq)
			if m: seq = seq[:m.start()]
			print(seq)
		elif line[0] == '@':
			sys.stdout.write(line)
			seq = next(fasta)[:-1]
			m = adapter_re.search(seq)
			trim_len = m.start() if m else len(seq)
			print(seq[:trim_len])
			sys.stdout.write(next(fasta))
			print(next(fasta)[:trim_len])








#################
# FASTA RNA2DNA #
#################

def fasta_rna_to_dna(fasta_path):
	fasta = zopen(fasta_path)
	for line in fasta:
		if line[0] in '#>@+':
			sys.stdout.write(line)
		else:
			sys.stdout.write(line.upper().replace('U', 'T'))
	fasta.close()





###############
# FASTA CHECK #
###############

def fasta_check(fasta_1_path, fasta_2_path, out_1_path, out_2_path):
	fasta_1 = zopen(fasta_1_path)
	fasta_2 = zopen(fasta_2_path)
	
	out_1 = zopen(out_1_path, 'w')
	out_2 = zopen(out_2_path, 'w')
	
	bad_out_1 = zopen('bad.' + out_1_path, 'w')
	bad_out_2 = zopen('bad.' + out_2_path, 'w')
	
	while 1:
		
		discard = False
		
		if fasta_1:
			while 1:
				line = fasta_1.readline()
				if line == '':
					fasta_1.close()
					fasta_1 = None
					break
				
				if not line[0] in '>@': continue
				
				header_1 = line[:-1].replace('/1', '')
				seq_1 = fasta_1.readline()[:-1]
				qual_1 = None
				if header_1[0] == '@':
					while line[0] != '+': line = fasta_1.readline()
					qual_1 = fasta_1.readline()[:-1]
					if len(seq_1) != len(qual_1):
						discard = True
				
				break
				
		if fasta_2:
			while 1:
				line = fasta_2.readline()
				if line == '':
					fasta_2.close()
					fasta_2 = None
					break
				
				if not line[0] in '>@': continue
					
				header_2 = line[:-1].replace('/2', '')
				seq_2 = fasta_2.readline()[:-1]
				if header_2[0] == '@':
					while line[0] != '+': line = fasta_2.readline()
					qual_2 = fasta_2.readline()[:-1]
					if len(seq_2) != len(qual_2):
						discard = True
				
				break
	
		if fasta_1 == None and fasta_2 == None: break
	
		if (fasta_1 == None) ^ (fasta_2 == None):
			info('File terminated abruptly.')
			break
		
		if header_1 != header_2: discard = True
		
		if discard:
			if qual_1 and qual_2:
				bad_out_1.write('%s/1\n%s\n+\n%s\n' % (header_1, seq_1, qual_1))
				bad_out_2.write('%s/2\n%s\n+\n%s\n' % (header_2, seq_2, qual_2))
			else:
				bad_out_1.write('%s/1\n%s\n' % (header_1, seq_1))
				bad_out_2.write('%s/2\n%s\n' % (header_2, seq_2))
		else:
			if qual_1 and qual_2:
				out_1.write('%s/1\n%s\n+\n%s\n' % (header_1, seq_1, qual_1))
				out_2.write('%s/2\n%s\n+\n%s\n' % (header_2, seq_2, qual_2))
			else:
				out_1.write('%s/1\n%s\n' % (header_1, seq_1))
				out_2.write('%s/2\n%s\n' % (header_2, seq_2))

	out_1.close()
	out_2.close()
		
	bad_out_1.close()
	bad_out_2.close()
		











################
# FASTA REPAIR #
################

def fasta_repair(fasta_1_path, fasta_2_path, out_1_path, out_2_path):
	fasta_1 = zopen(fasta_1_path)
	fasta_2 = zopen(fasta_2_path)
	
	out_1 = zopen(out_1_path, 'w')
	out_2 = zopen(out_2_path, 'w')
	
	orphans_1 = {}
	orphans_2 = {}
	
	while not (fasta_1 == None and fasta_2 == None):
		
		if fasta_1:
			while 1:
				line = fasta_1.readline()
				if line == '':
					fasta_1.close()
					fasta_1 = None
					break
				
				if not line[0] in '>@': continue
				
				header = line[:-1].replace('/1', '')
				seq = fasta_1.readline()[:-1]
				qual = None
				if header[0] == '@':
					while line[0] != '+': line = fasta_1.readline()
					qual = fasta_1.readline()[:-1]
					
				# Check that there are as many quality values as nucleotides.
				if len(seq) != len(qual):
					info('Read %s/1 discarded due to corrupted qualities.' %
						header[:-1])
					break
				
				read = orphans_2.get(header)
				if read:
					del orphans_2[header]
					if qual and read[1]:
						out_1.write('%s/1\n%s\n+\n%s\n' %
							(header, seq, qual))
						out_2.write('%s/2\n%s\n+\n%s\n' % 
							(header, read[0], read[1]))
					else:
						out_1.write('%s/1\n%s\n' % (header, seq))
						out_2.write('%s/2\n%s\n' % (header, read[0]))
				else:
					orphans_1[header] = (seq, qual)
					
				break
				
		if fasta_2:
			while 1:
				line = fasta_2.readline()
				if line == '':
					fasta_2.close()
					fasta_2 = None
					break
				
				if not line[0] in '>@': continue
					
				header = line[:-1].replace('/2', '')
				seq = fasta_2.readline()[:-1]
				qual = None
				if header[0] == '@':
					while line[0] != '+': line = fasta_2.readline()
					qual = fasta_2.readline()[:-1]
				
				# Check that there are as many quality values as nucleotides.
				if len(seq) != len(qual):
					info('Read %s/2 discarded due to corrupted qualities.' %
						header[:-1])
					break
				
				read = orphans_1.get(header)
				if read:
					del orphans_1[header]
					if qual and read[1]:
						out_1.write('%s/1\n%s\n+\n%s\n' %
							(header, read[0], read[1]))
						out_2.write('%s/2\n%s\n+\n%s\n' % 
							(header, seq, qual))
					else:
						out_1.write('%s/1\n%s\n' % (header, read[0]))
						out_2.write('%s/2\n%s\n' % (header, seq))
				else:
					orphans_2[header] = (seq, qual)
				
				break

			
	out_1.close()
	out_2.close()



			
###################
# FASTA COLOR2NUC #
###################
			
def fasta_color2nuc(color_seq):
	color_seq = color_seq.upper()
	
	prev_nuc = color_seq[0]
	if prev_nuc not in 'ACGTacgt':
		print('ERROR: Colorspace sequence must begin with a nucleotide.')
		exit(-1)
	
	decode = { 'A': 'ACGT', 'C': 'CATG', 'G': 'GTAC', 'T': 'TGCA' }
	
	nuc_seq = ''
	for k in range(1, len(color_seq)):
		prev_nuc = decode[prev_nuc][int(color_seq[k])]
		nuc_seq += prev_nuc
	
	return nuc_seq

		
		






##################
# FASTA FROM SRA #
##################

def fasta_from_sra(sra_path):
	shell('~/tools/sratoolkit*/fastq-dump --split-3 --gzip %s' % sra_path)







	

#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['filter']:
		fasta_filter(args['<fasta>'], args['<regexp>'])
	elif args['flatten']:
		fasta_flatten(args['<fasta>'], args['<out_dir>'])
	elif args['split'] and args['interleaved']:
		fasta_split_interleaved(args['<fasta>'], int(args['<tag_length>']))
	elif args['split']:
		fasta_split(args['<fasta>'], int(args['<tag_length>']),
			args['<anchors_5p_path>'], args['<anchors_3p_path>'])
	elif args['trim']:
		fasta_trim(args['<fasta>'], int(args['<trim_length>']))
	elif args['repair']:
		fasta_repair(args['<fasta_1>'], args['<fasta_2>'],
			args['<out_fasta_1>'], args['<out_fasta_2>'])
	elif args['check']:
		fasta_check(args['<fasta_1>'], args['<fasta_2>'],
			args['<out_fasta_1>'], args['<out_fasta_2>'])
	elif args['rna2dna']:
		fasta_rna_to_dna(args['<fasta>'])
	elif args['adapters'] and args['find']:
		fasta_find_adapters(args['<fasta>'])
	elif args['remove'] and args['adapters']:
		fasta_remove_adapters(args['<fasta>'], args['<adapter>'])
	elif args['from'] and args['sra']:
		fasta_from_sra(args['<sra_file>'])

	
