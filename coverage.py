#!/bin/env pypy

"""
Tools for copy number analysis and visualization.

Usage:
  coverage grid <genome_file> <window_size> [-s N]
  coverage cds <bam_file> <gtf_file>
  coverage telomere <bam_file>
  coverage downsample <wig_file> <fold>
  coverage sum <wig_files>...
  coverage unbias logratios <wig_file>
  coverage median filter <wig_file> <window_size>
  coverage format igv <wig_file>

Options:
  -h --help       Show this screen.
  -s --step=N     Step size for window placement [default: window size / 2].
"""

from __future__ import print_function
import sys, docopt, re, os, math, tempfile, shutil, subprocess
from pypette import zopen, shell, revcomplement, info, error, shell_stdout
from pypette import Object
from sam import read_sam, ref_sequence_sizes

class WigTrack: pass

def read_fixed_wig(wig_path):
	tracks = {}
	chr = ''
	for line in zopen(wig_path):
		if line.startswith('fixedStep'):
			if chr:
				track.values = values[0:N]	 # Remove preallocated space
			track = Object()
			values = np.zeros(1000000)
			m = re.search(r'chrom=(\w+)', line)
			chr = m.group(1)
			m = re.search(r'start=(\d+)', line)
			track.start = int(m.group(1))
			m = re.search(r'step=(\d+)', line)
			track.step = int(m.group(1))
			m = re.search(r'span=(\d+)', line)
			track.span = int(m.group(1)) if m else -1
			N = 0
			tracks[chr] = track
			continue
		
		if chr:
			values[N] = float(line)
			N += 1
	
	if chr: track.values = values[0:N]	 # Remove preallocated space
	return tracks
		



def parse_wig_header(line):
	m = re.search(r'chrom=(\w+)', line)
	chr = m.group(1)
	m = re.search(r'start=(\d+)', line)
	start = int(m.group(1))
	m = re.search(r'step=(\d+)', line)
	step = int(m.group(1))
	m = re.search(r'span=(\d+)', line)
	span = int(m.group(1)) if m else -1
	return (chr, start, step, span)





	


#################
# COVERAGE GRID #
#################

def coverage_grid(genome_path, winsize, step):
	for line in zopen(genome_path):
		if not line.strip(): continue
		c = line.rstrip('\n').split('\t')
		chr, chr_len = c[0], int(c[1])
		start = 1
		while start + winsize < chr_len:
			print('%s\t%d\t%d' % (chr, start - 1, start + winsize - 1))
			start += step









################
# COVERAGE CDS #
################

def coverage_cds(bam_path, gtf_path):
	
	chr_sizes = ref_sequence_sizes(bam_path)
	
	info('Constructing a map of coding regions...')
	coding = {}
	for chr, size in chr_sizes.iteritems():
		coding[chr] = [False] * size
	for line in zopen(gtf_path):
		if line.startswith('#'): continue
		cols = line.split('\t')
		if cols[2] != 'CDS': continue
		if len(cols[0]) > 5: continue   # Ignore chromosomes other than chrXX
		if not cols[0] in coding: continue
		coding[cols[0]][int(cols[3])-1:int(cols[4])] = True
		
	info('Calculating a coverage histogram...')
	coverage_hist = [0] * 200
	chr = ''
	pos = 0
	for line in shell_stdout('bedtools genomecov -d -split -ibam %s' % bam_path):
		cols = line.split('\t')
		if cols[0] != chr:
			chr = cols[0]
			cds = coding[chr]
			pos = int(cols[1])-2
			info('%s...' % chr)
		pos += 1
		if cds[pos]:
			coverage_hist[min(int(cols[2]), len(coverage_hist)-1)] += 1
			
	print('Coverage histogram:')
	print('===================')
	for cov in range(0, len(coverage_hist)):
		print('%d: %d' % (cov, coverage_hist[cov]))
	





#####################
# COVERAGE TELOMERE #
#####################

def coverage_telomere(bam_path):
	
	# Method is based on "Assessing telomeric DNA content in pediatric cancers
	# using whole-genome sequencing data" by Parker et al.
	telo_seq = 'TTAGGG' * 4
	rev_telo_seq = revcomplement(telo_seq)
	
	telo_count = 0
	for al in read_sam(bam_path, 'u'):
		if telo_seq in al[9] or rev_telo_seq in al[9]: telo_count += 1
	
	print('%s\t%d' % (bam_path, telo_count))

	
		
	
		
		
#######################
# COVERAGE DOWNSAMPLE #
#######################
		
def coverage_downsample(wig_path, fold):
	wig, step = read_fixed_wig(wig_path)
	mid = (fold-1)/2
	for chr, data in wig.iteritems():
		for k in xrange(data.size / fold):
			#data[k] = np.median(data[k*fold:k*fold+fold])    # FIXME
			data[k] = sorted(data[k*fold:k*fold+fold])[mid]
		data = data[0:data.size/fold]
		print('fixedStep chrom=%s start=1 step=%d span=%d' %
			(chr, step*fold, step*fold))
		for v in data: print('%.2f' % v)





################
# COVERAGE SUM #
################
		
def coverage_sum(wig_paths):
	wigs = [read_fixed_wig(p) for p in wig_paths]
	for chr in wigs[0]:
		total = wigs[0][chr].values.copy()
		for wig in wigs[1:]: total += wig[chr].values
		span = wigs[0][chr].span
		print('fixedStep chrom=%s start=%d step=%d%s' % (chr,
			wigs[0][chr].start, wigs[0][chr].step,
			(' span=%d' % span) if span > 0 else ''))
		for v in total: print(v)








#############################
# COVERAGE UNBIAS LOGRATIOS #
#############################
		
def coverage_unbias_logratios(wig_path):
	wig = read_fixed_wig(wig_path)

	# Find mode for each chromosome
	modes = {}
	for chr in wig:
		bins = np.arange(-4.975, 4.976, .05)
		hist = np.zeros(len(bins))
		for val in wig[chr].values:
			if not -5 <= val < 5: continue
			hist[int((val + 5) // 0.05)] += 1
		modes[chr] = bins[np.argmax(hist)]

	mode = sorted(modes.values())[len(modes)/2]

	for chr in wig:
		values = wig[chr].values
		values -= mode
		span = wig[chr].span
		print('fixedStep chrom=%s start=%d step=%d%s' % (chr,
			wig[chr].start, wig[chr].step,
			(' span=%d' % span) if span > 0 else ''))
		for v in values: print(v)








##########################
# COVERAGE MEDIAN FILTER #
##########################
		
def coverage_median_filter(wig_path, win_size):
	wig = read_fixed_wig(wig_path)

	for chr in wig:
		orig = wig[chr].values


		bins = np.arange(-4.975, 4.976, .05)
		hist = np.zeros(len(bins))
		for val in wig[chr].values:
			if not -5 <= val < 5: continue
			hist[int((val + 5) // 0.05)] += 1
		modes[chr] = bins[np.argmax(hist)]

	mode = sorted(modes.values())[len(modes)/2]

	for chr in wig:
		values = wig[chr].values
		values -= mode
		span = wig[chr].span
		print('fixedStep chrom=%s start=%d step=%d%s' % (chr,
			wig[chr].start, wig[chr].step,
			(' span=%d' % span) if span > 0 else ''))
		for v in values: print(v)




	

#######################
# COVERAGE FORMAT IGV #
#######################

def coverage_format_igv(wig_path):
	for line in zopen(wig_path):
		if line == 'nan\n': print('NaN')
		elif line == '-inf\n': print('-1000')
		elif line == 'inf\n': print('1000')
		else: sys.stdout.write(line)





#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['grid']:
		wsize = int(args['<window_size>'])
		step = wsize / 2
		if args['--step'].isdigit(): step = int(args['--step'])
		coverage_grid(args['<genome_file>'], wsize, step=step)
	elif args['cds']:
		coverage_cds(args['<bam_file>'], args['<gtf_file>'])
	elif args['telomere']:
		coverage_telomere(args['<bam_file>'])
	elif args['downsample']:
		coverage_downsample(args['<wig_file>'], int(args['<fold>']))
	elif args['sum']:
		coverage_sum(args['<wig_files>'])
	elif args['unbias'] and args['logratios']:
		coverage_unbias_logratios(args['<wig_file>'])
	elif args['median'] and args['filter']:
		coverage_median_filter(args['<wig_file>'], args['<window_size>'])
	elif args['format'] and args['igv']:
		coverage_format_igv(args['<wig_file>'])

