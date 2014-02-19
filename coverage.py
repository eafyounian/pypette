#!/bin/env pypy

"""
Tools for copy number analysis and visualization.

Usage:
  coverage tiled <bam_file> <window_size> [-S|-P|-M] [-s N] [-q N]
  coverage cds <bam_file> <gtf_file>
  coverage telomere <bam_file>
  coverage logratio <test_wig> <ref_wig> <min_ref>
  coverage downsample <wig_file> <fold>
  coverage format igv <wig_file>

Options:
  -h --help         Show this screen.
  -q --quality=N    Minimum alignment quality [default: 10].
  -s --step=N       Step size for window placement [default: window size].
  -S --single       Use all reads for coverage calculation, not just paired.
  -P --plus         Calculate coverage only for the plus strand.
  -M --minus        Calculate coverage only for the minus strand.
"""

from __future__ import print_function
import sys, docopt, re, os
from pypette import zopen, shell, revcomplement, info, error, shell_stdout
from sam import read_sam, ref_sequence_sizes
import numpypy as np





def read_fixed_wig(wig_path):
	header_re = re.compile('chrom=(\w+) .* step=(\d+)')
	tracks = {}
	chr = ''
	for line in zopen(wig_path):
		if line.startswith('fixedStep'):
			if chr: tracks[chr] = cov[0:N]  # Remove preallocated space
			m = header_re.search(line)
			chr = m.group(1)
			step = int(m.group(2))
			cov = np.zeros(1000000)
			N = 0
			continue
		
		if not chr: continue
		#if cov.size < N: cov.resize(cov.size * 2)   # Preallocate more space
		
		cov[N] = float(line)
		N += 1
	
	if chr: tracks[chr] = cov[0:N]   # Don't forget the last chromosome
	return (tracks, step)
		
	
		
		
		
		



	


##################
# COVERAGE TILED #
##################

def coverage_tiled(bam_path, window_size, quality, mode, step,
	max_frag_len=1000):
	
	chr_sizes = ref_sequence_sizes(bam_path)
	chr_cov = { chr: np.zeros(int(size / step) + 1, np.uint32)
		for chr, size in chr_sizes.iteritems() }
	
	win_overlap = window_size - step
	
	empty = np.array([])
	cov = empty
	chr = ''
	
	if mode == 'paired':
		# Only count concordant paired end reads.
		for al in read_sam(bam_path, 'A1', min_quality=quality):
			pos = int(al[3]); mpos = int(al[7])
			if al[6] != '=' or abs(pos - mpos) > max_frag_len: continue

			# Discard spliced and clipped reads.
			if 'N' in al[5] or 'S' in al[5]: continue

			start = (min(pos, mpos)-1 - win_overlap) / step
			stop = (max(pos, mpos)+len(al[9])-1 + win_overlap) / step
			if al[2] != chr:
				chr = al[2]
				cov = chr_cov.get(al[2], empty)
			if not cov.size: continue
			start = max(start, 0)
			stop = min(stop, cov.size-1)
			cov[start:stop+1] += 1
		
	else:
		sam_flags = {'single': 'a', 'plus': 'a+', 'minus': 'a-'}
			
		# Count all individual reads.
		for al in read_sam(bam_path, sam_flags[mode], min_quality=quality):
			pos = int(al[3])

			# Discard spliced and clipped reads.
			if 'N' in al[5] or 'S' in al[5]: continue

			start = (pos-1 - win_overlap) / step
			stop = (pos+len(al[9])-1 + win_overlap) / step
			if al[2] != chr:
				chr = al[2]
				cov = chr_cov.get(al[2], empty)
			if not cov.size: continue
			start = max(start, 0)
			stop = min(stop, cov.size-1)
			cov[start:stop+1] += 1

	
	for chr in chr_cov:
		print('fixedStep chrom=%s start=1 step=%d' % (chr, step))
		for x in chr_cov[chr]: print(x)









################
# COVERAGE CDS #
################

def coverage_cds(bam_path, gtf_path):
	
	chr_sizes = ref_sequence_sizes(bam_path)
	
	info('Constructing a map of coding regions...')
	coding = {}
	for chr, size in chr_sizes.iteritems():
		coding[chr] = np.zeros(size, np.bool_)
	for line in open(gtf_path):
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




#####################
# COVERAGE LOGRATIO #
#####################

def coverage_logratio(test_wig_path, ref_wig_path, min_ref=1):
	if min_ref <= 0:
		error('<min_ref> must be a positive number.')
		
	test, step = read_fixed_wig(test_wig_path)
	ref, step = read_fixed_wig(ref_wig_path)
	
	for chr in test:
		info(chr)
		#pos = np.where(valid) * step + 1
		logratios = np.log2(test[chr] / ref[chr])
		logratios[ref[chr] < min_ref] = np.nan
		print('fixedStep chrom=%s start=1 step=%d span=%d' % (chr, step, step))
		for v in logratios: print('%.2f' % v)
	
		
	
		
		
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
	if args['tiled']:
		wsize = int(args['<window_size>'])
		step = wsize if args['--step'] == 'window size' else int(args['--step'])
		
		mode = 'paired'
		if args['--single']: mode = 'single'
		if args['--plus']: mode = 'plus'
		if args['--minus']: mode = 'minus'
		
		coverage_tiled(args['<bam_file>'], wsize,
			quality=int(args['--quality']), mode=mode, step=step)
	elif args['cds']:
		coverage_cds(args['<bam_file>'], args['<gtf_file>'])
	elif args['telomere']:
		coverage_telomere(args['<bam_file>'])
	elif args['logratio']:
		coverage_logratio(args['<test_wig>'], args['<ref_wig>'],
			min_ref=float(args['<min_ref>']))
	elif args['downsample']:
		coverage_downsample(args['<wig_file>'], int(args['<fold>']))
	elif args['format'] and args['igv']:
		coverage_format_igv(args['<wig_file>'])

