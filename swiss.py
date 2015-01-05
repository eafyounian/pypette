#!/bin/env pypy

"""
Tools for miscellaneous tasks.

Usage:
  swiss rename <tsv_file>
  swiss link <tsv_file>
  swiss swap <file_1> <file_2>
  swiss xls2tsv <xls_file>
  swiss igv <tsv_file> <data_col>
  swiss download sra <sra_study>
  swiss wig2tsv <wig>
  swiss split wig [-H] <wig_file> <out_prefix>
  swiss ega checksum <tsv_file>
  swiss draw karyotype
  swiss annotate <input_file> <bed_file>
  swiss count segments <min> <max> <seg_files>...

Options:
  -h --help             Show this screen.
  -H --no-header        Discard the WIG headers.
"""

from __future__ import print_function
import sys, docopt, re, os, xlrd, csv, itertools
from collections import defaultdict
from datetime import datetime
from pypette import zopen, shell, read_fasta, revcomplement, GenomicFeatures
from pypette import info, error, Object



################
# SWISS RENAME #
################

def swiss_rename(tsv_path):
	
	for line in open(tsv_path, 'U'):
		line = line.replace('\n', '')
		tokens = line.split('\t')
		if len(tokens) != 2: continue
			
		(source, dest) = tokens 
		if not os.path.exists(source):
			info('Source file %s does not exist.' % source)
			continue
		if os.path.exists(dest):
			info('Destination file %s exists. Will not overwrite.' % dest)
			continue
		
		os.rename(source, dest)
		
		



################
# SWISS LINK #
################

def swiss_link(tsv_path):
	for line in open(tsv_path, 'U'):
		line = line.replace('\n', '')
		tokens = line.split('\t')
		if len(tokens) != 2: continue
			
		(source, dest) = tokens 
		if not os.path.exists(source):
			info('Source file %s does not exist.' % source)
			continue
		if os.path.lexists(dest):
			info('Destination file %s exists. Will not overwrite.' % dest)
			continue
		
		os.symlink(source, dest)






		
		
##############
# SWISS SWAP #
##############

def swiss_swap(file_1, file_2):
	tmp_file = file_2 + '.swaptmp'
	os.rename(file_1, tmp_file)
	os.rename(file_2, file_1)
	os.rename(tmp_file, file_2)







#################
# SWISS XLS2TSV #
#################

def swiss_xls2tsv(xls_file):
	book = xlrd.open_workbook(xls_file)
	sh = book.sheet_by_index(0)

	def read_xls_cell(cell):
		if cell.ctype == 2:
			text = unicode(cell.value)
			return text[:-2] if text[-2:] == '.0' else text
		if cell.ctype == 3: return datetime(*xlrd.xldate_as_tuple(
			cell.value, book.datemode)).date().isoformat()
		else:
			return cell.value.replace('\t', '    ').replace('\n', '\\n')

	for r in range(sh.nrows):
		#print([read_xls_cell(c) for c in sh.row(r)])
		line = '\t'.join(read_xls_cell(c) for c in sh.row(r))
		print(line.encode('UTF-8'))





#############
# SWISS IGV #
#############

def swiss_igv(tsv_path, data_col, one_based=True):
	tsv_file = zopen(tsv_path)
	headers = next(tsv_file)[:-1].split('\t')

	chrom_col = [i for i, h in enumerate(headers[:data_col])
		if re.match('chrom', h, re.I)]
	if len(chrom_col) != 1: error('Cannot find chromosome column.')
	chrom_col = chrom_col[0]

	pos_col = [i for i, h in enumerate(headers[:data_col])
		if re.match('pos', h, re.I)]
	if len(pos_col) != 1: error('Cannot find position column.')
	pos_col = pos_col[0]

	print('CHROMOSOME\tSTART\tEND\tFEATURE\t' + '\t'.join(headers[data_col:]))
	for line in tsv_file:
		tokens = line[:-1].split('\t')
		chr = tokens[chrom_col]
		pos = int(tokens[pos_col])
		if one_based: pos -= 1
		sys.stdout.write('%s\t%d\t%d\t-\t' % (chr, pos, pos+1))
		print('\t'.join(tokens[data_col:]))






######################
# SWISS DOWNLOAD SRA #
######################

def swiss_download_sra(sra_study):
	if not sra_study.startswith('SRP'):
		error('SRA study identifier must begin with "SRP".')
	
	shell('/data/csb/tools/ncftp-3.2.5/bin/ncftpget -R -v '
		'ftp-trace.ncbi.nlm.nih.gov ./ '
		'/sra/sra-instant/reads/ByStudy/sra/SRP/%s/%s' %
		(sra_study[:6], sra_study))
	



	
	
	
	
#################
# SWISS WIG2TSV #
#################

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
	'11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
	'21', '22', 'X', 'Y', 'M']
chr_index = {}
for idx, chr in enumerate(chromosomes):
	chr_index[chr] = idx + 1
	chr_index['chr'+chr] = idx + 1

def swiss_wig2tsv(wig_path):
	fixed_re = re.compile('fixedStep chrom=(\w+) start=(\d+) step=(\d+)')
	chr = ''
	pos = 0
	step = 0
	for line in zopen(wig_path):
		m = fixed_re.match(line)
		if m:
			chr = chr_index[m.group(1)]
			step = int(m.group(3))
			pos = int(m.group(2))
			continue
		print('%d\t%d\t%f' % (chr, pos, float(line)))
		pos += step








###################
# SWISS SPLIT WIG #
###################

def swiss_split_wig(wig_path, out_prefix, no_header=False):
	chrom_re = re.compile('chrom=(\w+)')
	chr = ''
	for line in zopen(wig_path):
		if line.startswith(('fixedS', 'variableS')):
			if chr: out.close()
			chr = chrom_re.search(line).group(1)
			out = zopen('%s_%s.wig.gz' % (out_prefix, chr), 'w')
			if not no_header: out.write(line)
			continue
		
		if chr: out.write(line)








######################
# SWISS EGA CHECKSUM #
######################

def swiss_ega_checksum(tsv_path):
	file = open(tsv_path, 'U')
	headers = file.next()[:-1].split('\t')
	print('\t'.join(headers))

	def checksum(path):
		with open(path) as f: return f.read().split()[0]
	
	if 'forward_file_name' in headers:
		fw_path = headers.index('forward_file_name')
		rev_path = headers.index('reverse_file_name')
		fw_gpg_checksum = headers.index('forward_file_md5')
		fw_plain_checksum = headers.index('forward_file_unencrypted_md5')
		rev_gpg_checksum = headers.index('reverse_file_md5')
		rev_plain_checksum = headers.index('reverse_file_unencrypted_md5')

		for line in file:
			tokens = line[:-1].split('\t')
			tokens[fw_gpg_checksum] = checksum(tokens[fw_path] + '.gpg.md5')
			tokens[fw_plain_checksum] = checksum(tokens[fw_path] + '.md5')
			tokens[rev_gpg_checksum] = checksum(tokens[rev_path] + '.gpg.md5')
			tokens[rev_plain_checksum] = checksum(tokens[rev_path] + '.md5')
			print('\t'.join(tokens))
			
	else:
		path = headers.index('file_name')
		gpg_checksum = headers.index('file_md5')
		plain_checksum = headers.index('file_unencrypted_md5')
		
		for line in file:
			tokens = line[:-1].split('\t')
			tokens[gpg_checksum] = checksum(tokens[path] + '.gpg.md5')
			tokens[plain_checksum] = checksum(tokens[path] + '.md5')
			print('\t'.join(tokens))




			
		



########################
# SWISS DRAW KARYOTYPE #
########################

def swiss_draw_karyotype():
	from svgfig import Rect, Fig, window
	chr = ''
	chr_y = 0
	clen = 0
	rects = []
	shades = {
		'gneg': '#ffffff', #(0,0%,100%)',
		'gpos25': '#c0c0c0',
		'gpos50': '#808080',
		'gpos75': '#404040',
		'gpos100': '#000000',
		'acen': '#ff0000',
		'gvar': '#ffffff',
		'stalk': '#ffffff',
	}

	def draw_border():
		rects.append(Rect(-w, chr_y, w, chr_y + clen, fill=None,
			stroke_width='0.05pt'))

	for line in open('/data/csb/organisms/homo_sapiens/cytoband.txt'):
		c = line.rstrip().split('\t')
		pos = (float(c[1]) / 1000000, float(c[2]) / 1000000)
		if c[0] != chr:
			if chr: draw_border()
			chr = c[0]
			chr_y += clen + 20
			clen = 0
		clen = max(clen, pos[1])
		w = 0.01 if c[4] == 'acen' else 0.02
		rects.append(Rect(-w, chr_y + pos[0], w, chr_y + pos[1],
			stroke='none', fill=shades[c[4]], stroke_linejoin='miter'))

	draw_border()

	Fig(*rects, trans=window(0, 20, 0, chr_y + clen, width=500)).SVG() \
		.save('/home/annalam/karyotype.svg')









##################
# SWISS ANNOTATE #
##################

def distance(pos, region):
	return max([0, region[0] - pos, pos - region[1]])

def swiss_annotate(input_path, bed_path):
	features = []
	for line in zopen(bed_path):
		c = line.rstrip().split('\t')
		features.append((c[0], c[5], (int(c[1]), int(c[2])), c[3]))

	for line in zopen(input_path):
		t = line.rstrip().split('\t')
		chr, pos = None, None

		m = re.match(r'(chr.+):(\d+)', t[0])
		if m:
			chr = m.group(1)
			pos = int(m.group(2))

		if chr == None:
			sys.stdout.write(line)
			continue
		
		nearby = []
		for f in features:
			if f[0] != chr: continue
			dist = distance(pos, f[2])
			if dist < 50000:
				nearby.append((re.sub(' \(ENSG.*?\)', '', f[3]), dist))

		nearby.sort(key=lambda x: x[1])
		sys.stdout.write(line[:-1])
		sys.stdout.write('\t')
		print(','.join(['%s (%d)' % f for f in nearby]))






########################
# SWISS COUNT SEGMENTS #
########################

def swiss_count_segments(seg_paths, range):

	samples = seg_paths

	# Read all segments into memory
	segments = defaultdict(lambda: defaultdict(list))
	for seg_path in seg_paths:
		seg_file = open(seg_path)
		line = next(seg_file)   # Discard headers
		sample_segs = segments[seg_path]
		for line in seg_file:
			cols = line.split('\t')
			sample_segs[cols[1]].append(
				(int(cols[2]), int(cols[3]), float(cols[4])))
		seg_file.close()

	# For each chromosome, find the union of segment breakpoints
	union = {}
	chromosomes = segments[segments.keys()[0]].keys()
	for chr in chromosomes:
		breakpoints = []
		for sample in segments:
			breakpoints += (seg[0] for seg in segments[sample][chr])
			breakpoints.append(max(seg[1] for seg in segments[sample][chr]))
		union[chr] = sorted(set(breakpoints))
		union[chr] = zip(union[chr][:-1], union[chr][1:])

	# Now collect segment values across all samples
	values = defaultdict(lambda: defaultdict(list))
	for sample in segments:
		for chr in chromosomes:
			val = values[sample][chr]
			segs = segments[sample][chr]
			for useg in union[chr]:
				match = [seg for seg in segs if
					useg[0] >= seg[0] and useg[1] <= seg[1]]
				if len(match) != 1:
					val.append(float('NaN'))
					continue
					#print([sample, chr, useg])
					#print(segs)
					#error('No encompassing segment found.')
				match = match[0]
				val.append(match[2])

	for chr in chromosomes:
		for n, useg in enumerate(union[chr]):
			matches = sum(range[0] <= values[sample][chr][n] <= range[1]
				for sample in values)
			print('%s\t%d\t%d\t%d' % (chr, useg[0], useg[1], matches))











#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__, options_first=True)
	if args['rename']:
		swiss_rename(args['<tsv_file>'])
	elif args['link']:
		swiss_link(args['<tsv_file>'])
	elif args['swap']:
		swiss_swap(args['<file_1>'], args['<file_2>'])
	elif args['xls2tsv']:
		swiss_xls2tsv(args['<xls_file>'])
	elif args['wig2tsv']:
		swiss_wig2tsv(args['<wig>'])
	elif args['igv']:
		swiss_igv(args['<tsv_file>'], int(args['<data_col>'])-1)
	elif args['download'] and args['sra']:
		swiss_download_sra(args['<sra_study>'])
	elif args['split'] and args['wig']:
		swiss_split_wig(args['<wig_file>'], args['<out_prefix>'],
			no_header=args['--no-header'])
	elif args['ega'] and args['checksum']:
		swiss_ega_checksum(args['<tsv_file>'])
	elif args['draw'] and args['karyotype']:
		swiss_draw_karyotype()
	elif args['annotate']:
		swiss_annotate(args['<input_file>'], args['<bed_file>'])
	elif args['count'] and args['segments']:
		swiss_count_segments(args['<seg_files>'],
			(float(args['<min>']), float(args['<max>'])))

