#!/bin/env pypy

"""
Tools for the manipulation of SAM and BAM files.

Usage:
  sam reads [-r] <bam_file> <out_prefix>
  sam unaligned reads <bam_file>
  sam compact <bam_file>
  sam discordant pairs [-q N] <bam_file> <min_distance_kb>
  sam extend fragments <bam_file> [--fragment-length=N]
  sam read length <bam_file>
  sam pileup each <vcf_file> <bam_files>... [--quality=N]
  sam pileup <region> <bam_files>... [--quality=N]
  sam count <bed_file> <bam_file>
  sam merge counts <bed_file> <count_files>...
  sam fragment lengths <bam_file>
  sam translate flags <flags>
  sam statistics <bam_files>...
  sam mismatch frequency [-q N] <bam_files>...

Options:
  -h --help             Show this screen.
  -r --raw              Output in raw sequence format.
  --fragment-length=N   Length to which single end reads will be extended
                        [default: 0].
  -p --parallel=N       Number of parallel processes to use [default: 8].
  -P --partition=PART   SLURM partition to run jobs on [default: local].
  -q --quality=N        Minimum alignment quality [default: 15].
  -s --strand=M         Strand specific read counts [default: none].
"""

from __future__ import print_function
import sys, subprocess, docopt, re, os
from pypette import zopen, shell, info, error, shell_stdout


def read_sam(sam_path, mode='', min_quality=0):
	view_options = ''
	flag_on = 0x0
	flag_off = 0x900       # Ignore secondary and supplementary alignments
	if 'a' in mode: flag_off |= 0x4                   # Aligned
	if 'A' in mode: flag_on |= 0x1; flag_off |= 0xc   # Aligned read pairs
	if 'u' in mode: flag_on |= 0x4                    # Unaligned
	if '1' in mode: flag_on |= 0x40                   # First mates
	if '2' in mode: flag_on |= 0x80                   # Second mates
	if '+' in mode: flag_off |= 0x10                  # Plus strand only
	if '-' in mode: flag_on |= 0x10                   # Minus strand only
	if not 'D' in mode: flag_off |= 0x400             # Flagged duplicates
	
	view_options += '-f 0x%x -F 0x%x ' % (flag_on, flag_off)
	
	if min_quality > 0: view_options += '-q%d ' % min_quality
	
	out = shell_stdout('samtools view %s %s' % (view_options, sam_path))
	for line in out:
		yield line.split('\t')


# This check is used in conjunction with samtools bam2fq, because some BAM
# files include mate suffixes in read IDs and some do not. This affects the
# use of the "-n" flag of samtools bam2fq.
def has_mate_suffixes(sam_path):
	first_ids = []
	for al in read_sam(sam_path):
		if len(first_ids) < 10: first_ids.append(al[0])
	return all(id.endswith('/1') or id.endswith('/2') for id in first_ids)


def ref_sequence_sizes(sam_path):
	out = shell_stdout('samtools view -H %s' % sam_path)
	chr_sizes = {}
	for line in out:
		m = re.match('@SQ\tSN:(\w+)\tLN:(\d+)', line)
		if m:
			chr_sizes[m.group(1)] = int(m.group(2))
	return chr_sizes




#############
# SAM READS #
#############

def sam_reads(bam_path, out_prefix):
	fastq_1 = zopen('%s_1.fq.gz' % out_prefix, 'w')
	fastq_2 = zopen('%s_2.fq.gz' % out_prefix, 'w')
	fastq = zopen('%s.fq.gz' % out_prefix, 'w')
	
	reads_1 = {}
	reads_2 = {}
		
	# FIXME: We assume that each read only has one alignment in the BAM file.
	for al in read_sam(bam_path):
		flags = int(al[1])		
		if flags & 0x40:
			rname = al[0][:-2] if al[0].endswith('/1') else al[0]
			mate = reads_2.pop(rname, None)
			if mate:
				fastq_1.write('@%s/1\n%s\n+\n%s\n' % (rname, al[9], al[10]))
				fastq_2.write('@%s/2\n%s\n+\n%s\n' % (rname, mate[0], mate[1]))
			else:
				reads_1[rname] = (al[9], al[10])
		elif flags & 0x80:
			rname = al[0][:-2] if al[0].endswith('/2') else al[0]
			mate = reads_1.pop(rname, None)
			if mate:
				fastq_1.write('@%s/1\n%s\n+\n%s\n' % (rname, mate[0], mate[1]))
				fastq_2.write('@%s/2\n%s\n+\n%s\n' % (rname, al[9], al[10]))
			else:
				reads_2[rname] = (al[9], al[10])	
		else:
			fastq.write('@%s\n%s\n+\n%s\n' % (al[0], al[9], al[10]))

	info('Found %d orphan first mates.' % len(reads_1))
	for read_id in reads_1.keys()[:5]: info('- Example: %s' % read_id)

	info('Found %d orphan second mates.' % len(reads_2))
	for read_id in reads_2.keys()[:5]: info('- Example: %s' % read_id)

	if len(reads_1) > 0:
		for rname, read in reads_1.iteritems():
			fastq.write('@%s\n%s\n+\n%s\n' % (rname, read[0], read[1]))
	
	if len(reads_2) > 0:
		for rname, read in reads_2.iteritems():
			fastq.write('@%s\n%s\n+\n%s\n' % (rname, read[0], read[1]))

	fastq_1.close()
	fastq_2.close()
	fastq.close()








#################
# SAM READS RAW #
#################

def sam_reads_raw(bam_path, out_prefix):
	out_1 = zopen('%s_1.reads.gz' % out_prefix, 'w')
	out_2 = zopen('%s_2.reads.gz' % out_prefix, 'w')
	out = zopen('%s.reads.gz' % out_prefix, 'w')
	
	reads_1 = {}
	reads_2 = {}

	# The "samtools bam2fq" command does not output supplementary or
	# secondary alignments. Each read only has one primary alignment.
	options = '-n' if has_mate_suffixes(bam_path) else ''
	bam2fq = shell_stdout('samtools bam2fq %s %s' % (options, bam_path))
	for line in bam2fq:
		if line[0] != '@': error('Invalid bam2fq output.')
		line = line[:-1]
		if line.endswith('/1'):
			segname = line[1:-2]
			mate = reads_2.pop(segname, None)
			if mate:
				out_1.write(next(bam2fq))
				out_2.write('%s\n' % mate)
			else:
				reads_1[segname] = next(bam2fq)[:-1]
		elif line.endswith('/2'):
			segname = line[1:-2]
			mate = reads_1.pop(segname, None)
			if mate:
				out_1.write('%s\n' % mate)
				out_2.write(next(bam2fq))
			else:
				reads_2[segname] = next(bam2fq)[:-1]
		else:
			out.write('%s\n' % next(bam2fq)[:-1])

		# Skip per-base qualities. They can start with '@'.
		next(bam2fq)
		next(bam2fq)

	info('Found %d orphan first mates.' % len(reads_1))
	for read_id in reads_1.keys()[:5]: info('- Example: %s' % read_id)

	info('Found %d orphan second mates.' % len(reads_2))
	for read_id in reads_2.keys()[:5]: info('- Example: %s' % read_id)

	if len(reads_1) > 0:
		for read in reads_1.itervalues(): out.write('%s\n' % read)
	
	if len(reads_2) > 0:
		for read in reads_2.itervalues(): out.write('%s\n' % read)

	out_1.close()
	out_2.close()
	out.close()











#######################
# SAM UNALIGNED READS #
#######################

def sam_unaligned_reads(bam_path):
	# The "samtools bam2fq" command does not output supplementary or
	# secondary alignments. Each read has max 1 primary alignment.
	options = '-n' if has_mate_suffixes(bam_path) else ''
	bam2fq = shell_stdout('samtools view -u -f 0x4 -F 0x900 %s | samtools bam2fq %s -' % (bam_path, options))
	for line in bam2fq:
		if line[0] != '@': error('Invalid bam2fq output.')
		sys.stdout.write('>')
		sys.stdout.write(line[1:])
		sys.stdout.write(next(bam2fq))

		# Skip per-base qualities. They can start with '@'.
		next(bam2fq)
		next(bam2fq)






###############
# SAM COMPACT #
###############

def sam_compact(bam_path):
	shell('samtools view -H %s' % bam_path)
	short_id = {}
	id_counter = 1
	for al in read_sam(bam_path, 'D'):
		rid = al[0]
		if rid[-2] == '/': rid = rid[:-2]
		new_id = short_id.get(rid)
		if new_id:
			del short_id[rid]
		else:
			new_id = str(id_counter)
			id_counter += 1
			short_id[rid] = new_id
		al[0] = new_id
		al[10] = '*'
		sys.stdout.write('\t'.join(al))





########################
# SAM DISCORDANT PAIRS #
########################

def sam_discordant_pairs(bam_path, min_distance, min_mapq=15):
	for al in read_sam(bam_path, 'A', min_quality=min_mapq):
		if al[6] == '=' and abs(int(al[7]) - int(al[3])) <= min_distance:
			continue
		if al[6] == '*': continue    # Sometimes unknown chromosomes show up
		if 'M' in al[2] or 'M' in al[6]: continue    # Discard mitochondrial
		sys.stdout.write('\t'.join(al))






########################
# SAM EXTEND FRAGMENTS #
########################

def sam_extend_fragments(bam_path, extend_len=0):
	frag_id = ''
	frag_pieces = 0
	
	out = shell_stdout('samtools sort -on %s %s | bedtools bamtobed -i stdin' %
		(bam_path, bam_path))
	for line in out:
		tokens = line.split('\t')
			
		read_id = tokens[3]
		if read_id[-2] == '/': read_id = read_id[:-2]
		
		if read_id == frag_id:
			frag_pieces += 1
			if tokens[0] == frag_chr:
				frag_start = min(frag_start, int(tokens[1]))
				frag_end = max(frag_end, int(tokens[2]))
			else:
				frag_id = ''    # Fragment is interchromosomal, discard it
		
		else:
			if frag_pieces == 1 and extend_len > 0:
				if frag_strand == '+':
					frag_end = frag_start + extend_len - 1
				else:
					frag_start = frag_end - extend_len + 1
			
			if frag_pieces > 0:
				print('%s\t%d\t%d' % (frag_chr, frag_start - 1, frag_end))
			
			frag_id = read_id
			frag_chr = tokens[0]
			frag_strand = tokens[5]
			frag_start = int(tokens[1])
			frag_end = int(tokens[2])
			frag_pieces = 1
	
			




###################
# SAM READ LENGTH #
###################

def read_length(bam_path):
	read_lens = []
	for al in read_sam(bam_path, 'a'):
		if len(read_lens) >= 100: break
		read_lens.append(len(al[9]))

	if len(set(read_lens)) > 1:
		error('SAM file contains reads of varying length.')
	else:
		return read_lens[0]





##############
# SAM PILEUP #
##############
	
def sam_pileup(region, bam_paths, min_al_quality=0):
	
	# Check the file paths here to ensure a nicer error message if files are
	# missing.
	missing = [path for path in bam_paths if not os.path.isfile(path)]
	for path in missing:
		info('WARNING: File %s was not found.' % path)
	bam_paths = [path for path in bam_paths if os.path.isfile(path)]
	if not bam_paths: return
	
	chr, region = region.replace(' ', '').split(':')
	region = [int(x) for x in region.split('-')]
	if len(region) == 1: region *= 2

	dev_null = open('/dev/null', 'a')
	indel_rx = re.compile('(\w[+-]\d+)?(\w+)(?![+-])')

	for pos in range(region[0], region[1]+1):
		if region[0] != region[1]: print('Pileup for %s:%d:' % (chr, pos))
		
		for bam in bam_paths:
			line = subprocess.check_output(
				'samtools mpileup -A -B -q%d -r %s:%d-%d %s' %
				(min_al_quality, chr, pos, pos, bam), shell=True,
				stderr=dev_null)
			sample_name = re.sub(r'(.*/)?(.*).bam', r'\2', bam)
			if not line:
				print('%s\t' % sample_name)
			else:
				tokens = line[:-1].split('\t')
				bases = re.sub(r'\^.', '', tokens[4]).upper()
				bases = re.sub(r'[$<>]', '', bases)
				
				# Parse the pileup string for indels
				indel_tokens = indel_rx.findall(bases)
				bases = ''.join([m[1][int(m[0][2:]):] if m[0] else m[1]
					for m in indel_tokens])
				indels = [m[0][:2] + m[1][:int(m[0][2:])]
					for m in indel_tokens if m[0]]
				
				bases = ''.join(sorted(bases))
				if bases: bases += ' '
				print('%s\t%s%s' % (sample_name, bases, ' '.join(indels)))

	dev_null.close()







###################
# SAM PILEUP EACH #
###################
	
def sam_pileup_each(vcf_path, bam_paths, min_al_quality=0):
	vcf = zopen(vcf_path)
	#for line in vcf:
	#	if line.startswith('CHROM'): break
	
	#headers = line.strip().split('\t')
	#ref_allele_col = headers.index(re.findall('ref\w*', line, re.I)[0])
	#alt_allele_col = headers.index(re.findall('alt\w*', line, re.I)[0])
	#nearby_genes_col = []
	#if re.search('nearby\w*', line, re.I):
	#	nearby_genes_col = headers.index(re.findall('nearby\w*', line, re.I)[0])
	
	for line in vcf:
		if not line.startswith('chr'): continue
		tokens = line.split('\t')
		#print('Mutation %s:%s:%s>%s (%s):' % (tokens[0], tokens[1],
		#	tokens[ref_allele_col], tokens[alt_allele_col],
		#	tokens[nearby_genes_col] if nearby_genes_col else ''))
		sys.stdout.write(line)
		sam_pileup('%s:%s' % (tokens[0], tokens[1]), bam_paths,
			min_al_quality=min_al_quality)
		print()










#############
# SAM COUNT #
#############

def sam_count(bam_path, bed_path,
	genome_path='~/organisms/homo_sapiens/hg19.chrom.sizes'):
	shell('bedtools coverage -split -sorted -counts -g %s -a %s -b %s | cut -f4' %
		(genome_path, bed_path, bam_path))









####################
# SAM MERGE COUNTS #
####################

def sam_merge_counts(bed_path, count_paths):
	samples = [p.replace('.tsv', '') for p in count_paths]
	bed_file = open(bed_path)
	cols = next(bed_file).rstrip('\n').split('\t')
	header = ['CHROMOSOME', 'START', 'END']
	if len(cols) >= 4: header.append('FEATURE')
	while len(header) < len(cols): header.append('')
	header += samples
	print('\t'.join(header))
	for line in shell_stdout('paste %s %s' % (bed_path,' '.join(count_paths))):
		sys.stdout.write(line)






	
	
	



########################
# SAM FRAGMENT LENGTHS #
########################

def sam_fragment_lengths(bam_path):
	max_size = 10000
	hist = [0] * (max_size + 1)
	for al in read_sam(bam_path, 'A1'):
		if al[6] != '=': continue
		if 'N' in al[5] or 'S' in al[5]: continue   # Spliced and clipped
		
		fragsize = abs(int(al[7]) - int(al[3])) + len(al[9])
		fragsize = min(fragsize, max_size)
		hist[fragsize] += 1
	
	for fragsize, count in enumerate(hist):
		print('%d\t%d' % (fragsize, count))

		







	
#######################
# SAM TRANSLATE FLAGS #
#######################

def sam_flags(flags):
	flags = int(flags)
	print('- Paired end' if flags & 0x1 else '- Single end')
	if flags & 0x1:
		if flags & 0x40:
			print('  * Read is first mate')
		if flags & 0x80:
			print('  * Read is second mate')
		
		print('  * Mate is unmapped' if flags & 0x8 else '  * Mate is mapped')
		print('  * Mate aligned to reverse complement of reference'
			if flags & 0x20 else '  * Mate aligned to forward reference')
		
	print('- Read is unmapped' if flags & 0x4 else '- Read is mapped')
	print('- Aligned to reverse complement of reference' if flags & 0x10 else
		'- Aligned to forward reference')
	if flags & 0x100:
		print('- Secondary alignment')
	if flags & 0x200:
		print('- Did not pass quality controls')
	if flags & 0x400:
		print('- PCR or optical duplicate')
	if flags & 0x800:
		print('- Supplementary alignment')
	









##################
# SAM STATISTICS #
##################

def sam_statistics(bam_paths):
	samples = [re.sub('\.bam$', '', s) for s in bam_paths]
	print('SAMPLE\tTOTAL READS\tALIGNED READS\tALIGNED READS WITH ALIGNED MATE\tALIGNED READS WITH CONCORDANT MATE\tMITOCHONDRIAL')

	for bam_path in bam_paths:
		total = -1
		aligned = -1
		aligned_with_aligned_mate = -1
		aligned_with_concordant_mate = -1
		mitochondrial = -1

		for line in shell_stdout('samtools flagstat %s' % bam_path):
			m = re.search(r'(\d+) \+ (\d+) in total', line)
			if m: total = int(m.group(1)) + int(m.group(2))
			
			m = re.search(r'(\d+) \+ (\d+) mapped', line)
			if m: aligned = int(m.group(1)) + int(m.group(2))

			m = re.search(r'(\d+) \+ (\d+) properly paired', line)
			if m: aligned_with_concordant_mate = \
				int(m.group(1)) + int(m.group(2))

			m = re.search(r'(\d+) \+ (\d+) with itself and mate mapped', line)
			if m: aligned_with_aligned_mate = int(m.group(1)) + int(m.group(2))

		# Count the number of reads aligned to mitochondrial DNA
		for line in shell_stdout('samtools view -c %s chrM' % bam_path):
			mitochondrial = int(line)
			break

		print('%s\t%d\t%d (%.1f%%)\t%d (%.1f%%)\t%d (%.1f%%)\t%d (%.1f%%)' % (
			re.sub('\.bam$', '', bam_path), total,
			aligned, float(aligned) / total * 100,
			aligned_with_aligned_mate,
			float(aligned_with_aligned_mate) / total * 100,
			aligned_with_concordant_mate,
			float(aligned_with_concordant_mate) / total * 100,
			mitochondrial, float(mitochondrial) / total * 100)) 








##########################
# SAM MISMATCH FREQUENCY #
##########################

def sam_mismatch_frequency(bam_paths, min_mapq):
	print('SAMPLE\tTOTAL\tMISMATCH\tFRACTION')
	for bam_path in bam_paths:
		total_mm = 0
		total_nucs = 0
		xm_col = 0
		for al in read_sam(bam_path, 'a', min_quality=min_mapq):
			if not al[xm_col].startswith('XM:i:'):
				xm_col = [c.startswith('XM:i:') for c in al].index(True)
			total_mm += int(al[xm_col][5:])
			total_nucs += len(al[9])
		print('%s\t%d\t%d\t%.3f%%' % (bam_path.replace('.bam', ''), total_nucs,
			total_mm, float(total_mm) / total_nucs * 100))









#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['unaligned'] and args['reads']:
		sam_unaligned_reads(args['<bam_file>'])
	elif args['reads'] and args['--raw']:
		sam_reads_raw(args['<bam_file>'], args['<out_prefix>'])
	elif args['reads']:
		sam_reads(args['<bam_file>'], args['<out_prefix>'])
	elif args['compact']:
		sam_compact(args['<bam_file>'])
	elif args['discordant'] and args['pairs']:
		sam_discordant_pairs(args['<bam_file>'],
			int(args['<min_distance_kb>']) * 1000,
			min_mapq=int(args['--quality']))
	elif args['extend'] and args['fragments']:
		sam_extend_fragments(args['<bam_file>'],
			extend_len=int(args['--fragment-length']))
	elif args['read'] and args['length']:
		read_len = read_length(args['<bam_file>'])
		if not read_len: error('Could not determine read length.')
		else: print('%d' % read_len)
	elif args['pileup'] and args['each']:
		sam_pileup_each(args['<vcf_file>'], args['<bam_files>'],
			min_al_quality=int(args['--quality']))
	elif args['pileup']:
		sam_pileup(args['<region>'], args['<bam_files>'],
			min_al_quality=int(args['--quality']))
	elif args['count']:
		sam_count(args['<bam_file>'], args['<bed_file>'])
	elif args['counts'] and args['merge']:
		sam_merge_counts(args['<bed_file>'], args['<count_files>'])
	elif args['fragment'] and args['lengths']:
		sam_fragment_lengths(args['<bam_file>'])
	elif args['flags']:
		sam_flags(args['<flags>'])
	elif args['statistics']:
		sam_statistics(args['<bam_files>'])
	elif args['mismatch'] and args['frequency']:
		sam_mismatch_frequency(args['<bam_files>'],
			min_mapq=int(args['--quality']))

