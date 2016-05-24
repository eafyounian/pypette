#!/bin/env pypy

"""
Calculate small RNA expression based on high throughput sequencing data.

Usage:
  smallrna expression <srna_fasta> <reads>...
  smallrna expression bgi <srna_fasta> <count_files>...
  smallrna parse mirbase <mirbase_gff>

Options:
  -h --help      Show this screen.
  --adapter=SEQ  Adapter sequence [default: ''].
"""

from __future__ import print_function
import sys, subprocess, docopt, re, datetime, pwd, os
from collections import Counter, defaultdict
from pypette import read_fasta, error, info, zopen

def write_fasta(seq_dict, fasta_path):
	fasta = open(fasta_path, 'w')
	for k, v in seq_dict.iteritems():
		fasta.write('>%s\n%s\n' % (k, v))
	fasta.close()
	






def smallrna_expression(read_paths, srna_reference_path):
	
	S = len(read_paths)
	min_other_reads = 100
	
	# Read the FASTA file containing small RNA reference sequences and
	# construct a new FASTA file that includes potential isoforms and
	# variants of these small RNA sequences.
	info('Constructing database of reference small RNA sequences...')
	seq_names = defaultdict(lambda: '')
	counts = defaultdict(lambda: [0]*S)

	for name, seq in read_fasta(srna_reference_path).iteritems():
		name = re.sub(' MIMAT.*', '', name)
		seq = seq.upper().replace('U', 'T')
		seq_names[seq] = name
		seq_names[seq[:-1]] = name + '-1'
		seq_names[seq + 'A'] = name + '+A'
		seq_names[seq + 'C'] = name + '+C'
		seq_names[seq + 'G'] = name + '+G'
		seq_names[seq + 'T'] = name + '+T'
	
	info('Counting reads aligning to small RNA sequences...')
	for s, read_path in enumerate(read_paths):
		fasta = zopen(read_path)
		for line in fasta:
			if not line or line[0] == '#': continue
			if line[0] in '>@':
				seq = next(fasta)[:-1]
				counts[seq][s] += 1

	counts = { seq: count for seq, count in counts.iteritems()
		if seq in seq_names or sum(count >= min_other_reads) >= 2 }
	
	print('NAME\tSEQUENCE\t%s' % '\t'.join(read_paths))
	for seq in sorted(counts.iterkeys(), key=lambda x: seq_names[x]):
		sys.stdout.write('%s\t%s' % (seq_names[seq], seq))
		for x in counts[seq]: sys.stdout.write('\t%d' % x)
		sys.stdout.write('\n')




def smallrna_expression_bgi(count_paths, srna_reference_path):
	
	S = len(count_paths)
	min_other_reads = 100
	
	# Read the FASTA file containing small RNA reference sequences and
	# construct a new FASTA file that includes potential isoforms and
	# variants of these small RNA sequences.
	seq_names = defaultdict(lambda: '')
	counts = defaultdict(lambda: [0]*S)
	for name, seq in read_fasta(srna_reference_path).iteritems():
		name = re.sub(' MIMAT.*', '', name)
		seq = seq.upper().replace('U', 'T')
		seq_names[seq] = name
		seq_names[seq[:-1]] = name + '-1'
		seq_names[seq + 'A'] = name + '+A'
		seq_names[seq + 'C'] = name + '+C'
		seq_names[seq + 'G'] = name + '+G'
		seq_names[seq + 'T'] = name + '+T'
	
	for s, count_path in enumerate(count_paths):
		for line in open(count_path):
			tokens = line.strip().split('\t')
			count = int(tokens[1])
			if tokens[2] in seq_names or count >= min_other_reads:
				counts[tokens[2]][s] += count
	
	print('NAME\tSEQUENCE\t%s' % '\t'.join(count_paths))
	for seq in sorted(counts.iterkeys(), key=lambda x: seq_names[x]):
		if seq_names[seq] == '' and sum(c >= min_other_reads for c in counts[seq]) < 2:
			continue
		sys.stdout.write('%s\t%s' % (seq_names[seq], seq))
		for x in counts[seq]: sys.stdout.write('\t%d' % x)
		sys.stdout.write('\n')





def small_rna_expression(reads_path, srna_reference_path, adapter=''):
	
	tmp_path = ''
	if not tmp_path:
		tmp_path = '/data/tmp/%s/small_rna_expression_%s' % (
			pwd.getpwuid(os.getuid()).pw_name,
			datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))
		os.makedirs(tmp_path)
	
	# Read the FASTA file containing small RNA reference sequences and
	# construct a new FASTA file that includes potential isoforms and
	# variants of these small RNA sequences.
	small_rna = {}
	srna_ref = open(srna_reference_path)
	for line in srna_ref:
		if not line[0] == '>': continue
		srna_name = line[1:-1]
		srna_name = re.sub(' .*', '', srna_name)
		srna_seq = next(srna_ref)[:-1].upper().replace('U', 'T')
		
		small_rna[srna_name] = srna_seq + adapter
		small_rna[srna_name + '-1'] = srna_seq[:-1] + adapter
		small_rna[srna_name + '+A'] = srna_seq + 'A' + adapter
		small_rna[srna_name + '+C'] = srna_seq + 'C' + adapter
		small_rna[srna_name + '+G'] = srna_seq + 'G' + adapter
		small_rna[srna_name + '+T'] = srna_seq + 'T' + adapter
		
	srna_isoforms_fasta = tmp_path + '/small_rna_isoforms.fa'
	srna_isoforms_index = tmp_path + '/small_rna_isoforms'
	write_fasta(small_rna, srna_isoforms_fasta)
	
	colorspace = True
	
	subprocess.check_call('bowtie-build -q %s %s %s' % 
		(' -C' if colorspace else '', srna_isoforms_fasta, srna_isoforms_index),
		shell=True)
	
	alignments_path = tmp_path + '/alignments.bowtie'
	# FIXME: Take out the trimming part
	subprocess.check_call(
		'bowtie %s -f -p8 --best --trim3 8 --norc -v2 -k1 -m1 %s %s '
		'--suppress=1,2,7,8 > %s' % 
		(' -C' if colorspace else '', srna_isoforms_index, reads_path,
		alignments_path), shell=True)
		
	srna_reads = Counter({key: 0 for key in small_rna})
	
	alignments = open(alignments_path)
	for line in alignments:
		tokens = line[:-1].split('\t')
		
		# Check that the microRNA does not have too much 5' end variation.
		if int(tokens[1]) > 3: continue
			
		srna_reads[tokens[0]] += 1
	
	srna_expr_path = 'srna_expr.tsv'
	srna_expr = open(srna_expr_path, 'w')
	srna_expr.write('NAME\tNUM_READS\n')
	for key in sorted(srna_reads.iterkeys()):
		srna_expr.write('%s\t%d\n' % (key, srna_reads[key]))
	
	
	




##########################
# SMALLRNA PARSE MIRBASE #
##########################

def smallrna_parse_mirbase(mirbase_gff_path):
	info('Printing mature microRNA loci in BED format.')
	
	mirna_name_re = re.compile(r';Name=([^\s;]+)')
	
	for line in open(mirbase_gff_path):
		if line[0] == '#': continue
		tokens = line[:-1].split('\t')
		if tokens[2] != 'miRNA': continue
			
		print('%s\t%d\t%d\t%s' % (tokens[0], int(tokens[3])-1, int(tokens[4]),
			mirna_name_re.search(tokens[8]).group(1)))
	
	
		
	
	
if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['expression'] and args['bgi']:
		smallrna_expression_bgi(args['<count_files>'],args['<srna_fasta>'])
	elif args['expression']:
		smallrna_expression(args['<reads>'], args['<srna_fasta>'])
	elif args['parse'] and args['mirbase']:
		smallrna_parse_mirbase(args['<mirbase_gff>'])

