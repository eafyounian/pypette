#!/bin/env pypy

'''
BreakFast is a toolkit for detecting chromosomal rearrangements
based on whole genome sequencing data.

Usage:
breakfast detect <bam_file> <genome> <out_prefix> [-a N] [-f N] [-q N]
	[-d N] [-O orientation]
breakfast detect specific [-A] <bam_file> <donors> <acceptors> <genome>
	<out_prefix>
breakfast filter <sv_file> [-r P-S-A]... [--blacklist=PATH]
breakfast annotate <sv_file> <bed_file>
breakfast blacklist [--freq-above=FREQ] <sv_files>...
breakfast visualize <sv_file>
breakfast tabulate rearranged genes <sv_files>...
breakfast tabulate fusions <sv_files>...
breakfast statistics <sv_files>...
breakfast filter by region <sv_file> <region>
breakfast filter by distance <min_distance> <sv_file>
breakfast align junction <reads>

Options:
  -a, --anchor-len=N      Anchor length for split read analysis. When zero,
                          split reads are not used [default: 0].
  -f, --max-frag-len=N    Maximum fragment length [default: 5000].
  -q, --min-mapq=N        Minimum mapping quality to consider [default: 15].
  -d, --min-distance=N    Min kb distance between breakpoints [default: 10].
  -O, --orientation=OR    Read pair orientation produced by sequencer. Either
                          'fr' (converging), 'rf' (diverging) or 'ff'
                          [default: fr].
  -A, --all-reads         Use all reads for rearrangement detection, not just
                          unaligned reads.

  --blacklist <list>      Path to a file containing blacklisted regions.
  -r, --min-reads=P-S-A   Minimum number of spanning reads required to accept
                          a breakpoint. Specified in the format P-S-A, where
                          P=paired, S=split, A=either. For example, -r 1-2-0
                          would require at least one mate pair and two split
                          reads of evidence [default: 0-0-0].
							
  --freq-above=FREQ       Minimum frequency at which a variant must be
                          present among the control samples to be
                          considered a false positive [default: 0].
'''

from __future__ import print_function
import sys, re, docopt, itertools, os, sam, gc
from collections import defaultdict
from pypette import info, error, shell_stdout, shell, zopen, Object
from pypette import read_flat_seq, revcomplement, point_region_distance
from pypette import regions_from_bed, read_fasta

class Rearrangement(object):
	__slots__ = ('chr', 'strand', 'pos', 'mchr', 'mstrand', 'mpos', 'reads')
	
	def __init__(self, chr, strand, pos, mchr, mstrand, mpos, read):
		self.chr = chr
		self.strand = strand
		self.pos = (pos, pos)
		self.mchr = mchr
		self.mstrand = mstrand
		self.mpos = (mpos, mpos)
		self.reads = [read]
		
	def id(self):
		return '%s:%s:%d <-> %s:%s:%d' % \
			(self.chr, self.strand, self.pos,
			self.mchr, self.mstrand, self.mpos)


# Columns 1-4 describe the breakpoint with the lower coordinate.
# Columns 6-9 describe the breakpoint with the higher coordinate.
# Columns 11-13 describe the evidence for the breakpoint.
# If the breakpoint location is only based on spanning fragments, the
# position represents the nucleotide where the mate closest to the
# breakpoint ends.
sv_file_header = (
	'CHROM\tSTRAND\tPOSITION\tNEARBY_FEATURES\t\t'
	'CHROM\tSTRAND\tPOSITION\tNEARBY_FEATURES\t\t'
	'NUM_SPANNING_FRAGMENTS\tNUM_SPANNING_MATES\tSPANNING_MATE_SEQUENCES')







####################
# BREAKFAST DETECT #
####################

def detect_discordant_pairs(sam_path, out_prefix, min_rearrangement_size,
	min_mapq, orientation):
	
	out = zopen(out_prefix + '.discordant_pairs.tsv.gz', 'w')
	N = 0

	sort_tmp_dir = os.path.dirname(out_prefix)
	if not sort_tmp_dir: sort_tmp_dir = './'
	
	# Go through all the first mates and look for discordant pairs.
	info('Searching for discordant read pairs...')
	prev = ['']
	for line in shell_stdout(
		'sam discordant pairs -q%d %s %d | sort -k1,1 -T %s' %
		(min_mapq, sam_path, min_rearrangement_size / 1000, sort_tmp_dir)):
		
		al = line.split('\t')
		if len(al) < 9: continue
		
		# Discard spliced and clipped reads.
		# FIXME: Add support for spliced RNA-seq reads.
		if 'N' in al[5] or 'S' in al[5]: continue
		
		if al[0][-2] == '/': al[0] = al[0][:-2]   # Remove /1 or /2 suffix
		
		if al[0] != prev[0]:
			prev = al
			continue
		
		flags = int(al[1])
		chr = al[2]
		mchr = prev[2]
		strand = '-' if flags & 0x10 else '+'
		mstrand = '-' if flags & 0x20 else '+'
		pos = int(al[3])
		mpos = int(prev[3])
		rlen = len(al[9])
		mrlen = len(prev[9])

		if not chr.startswith('chr'): chr = 'chr' + chr
		if not mchr.startswith('chr'): mchr = 'chr' + mchr

		if orientation == 'fr':
			# Reorient pairs so that the first mate is always upstream.
			if chr > mchr or (chr == mchr and pos > mpos):
				chr, mchr = mchr, chr
				pos, mpos = mpos, pos
				rlen, mrlen = mrlen, rlen
				strand, mstrand = mstrand, strand

			# Convert to forward-forward orientation (flip second mate).
			mstrand = '-' if mstrand == '+' else '+'

		elif orientation == 'rf':
			# Reorient pairs so that the first mate is always upstream.
			if chr > mchr or (chr == mchr and pos > mpos):
				chr, mchr = mchr, chr
				pos, mpos = mpos, pos
				rlen, mrlen = mrlen, rlen
				strand, mstrand = mstrand, strand

			# Convert to forward-forward orientation (flip first mate).
			strand = '-' if strand == '+' else '+'

		elif orientation == 'ff':
			# Reorient pairs so that the first mate is always upstream.
			# If mates are swapped, both mates must be reversed.
			if chr > mchr or (chr == mchr and pos > mpos):
				chr, mchr = mchr, chr
				pos, mpos = mpos, pos
				rlen, mrlen = mrlen, rlen
				strand, mstrand = '+' if mstrand == '-' else '-', \
					'+' if strand == '-' else '-'

		else:
			error('Unsupported read orientation detected.')
		
		# Make positions represent read starts.
		if strand == '-': pos += rlen - 1
		if mstrand == '-': mpos += mrlen - 1

		# Each discordant mate pair is represented as a 7-tuple
		# (chr_1, strand_1, pos_1, chr_2, strand_2, pos_2, None).
		# The None at the end signifies that this is a mate pair.
		# Positions are 1-based and represent read starts.
		out.write('%s\t%s\t%d\t%s\t%s\t%d\t-\n' % (
			chr, strand, pos, mchr, mstrand, mpos))
		N += 1
	
	out.close()
	info('Found %d discordant mate pairs.' % N)
	
	




def detect_discordant_mates(sam_path, genome_path, out_prefix, anchor_len,
	min_rearrangement_size):
	
	out = zopen(out_prefix + '.discordant_singles.tsv.gz', 'w')
	N = 0
	
	info('Splitting unaligned reads into %d bp anchors and aligning against '
		'the genome...' % anchor_len)
	
	# IMPORTANT: Only one thread can be used, otherwise alignment order is not
	# guaranteed and the loop below will fail.
	anchor_alignments = shell_stdout(
		'sam unaligned reads %s | fasta split interleaved - %d | '
		'bowtie -f -p1 -v0 -m1 -B1 --suppress 5,6,7,8 %s -'
		% (sam_path, anchor_len, genome_path))
	
	# FIXME: Would be nice to do this one chromosome at a time, so we wouldn't
	# need to use 3 gigabytes of memory.
	chromosomes = read_flat_seq(genome_path)
	for chr in list(chromosomes.keys()):
		if not chr.startswith('chr'):
			chromosomes['chr' + chr] = chromosomes.pop(chr)
	
	prev = ['']
	for line in anchor_alignments:
		
		al = line.split('\t')
		if al[0][-2] == '/': al[0] = al[0][:-2]
		
		if al[0] != prev[0]:
			prev = al
			continue
		
		chr = prev[2]
		mchr = al[2]
		strand = prev[1]
		mstrand = al[1]
		pos = int(prev[3])
		mpos = int(al[3])
		seq = prev[0][prev[0].find('_')+1:]
		full_len = len(seq)

		if not chr.startswith('chr'): chr = 'chr' + chr
		if not mchr.startswith('chr'): mchr = 'chr' + mchr

		# Ignore anchor pairs where the anchors are too close.
		if chr == mchr and abs(pos - mpos) < min_rearrangement_size:
			continue
			
		# Ignore rearrangements involving mitochondrial DNA.
		if 'M' in chr or 'M' in mchr: continue
			
		# Reorient the pairs so the first anchor is always upstream.
		# If mates are swapped, both mates must be reverse-complemented.
		if chr > mchr or (chr == mchr and pos > mpos):
			chr, mchr = mchr, chr
			pos, mpos = mpos, pos
			strand, mstrand = '+' if mstrand == '-' else '-', \
				'+' if strand == '-' else '-'
			seq = revcomplement(seq)
		
		# Extract the flanking sequences from the chromosome sequences.
		# The range calculations are a bit complex. It's easier to understand
		# them if you first add one to all indices to convert to 1-based
		# genomic coordinates ("pos" and "mpos" are 1-based).
		if strand == '+':
			left_grch = chromosomes[chr][pos-1:pos+full_len-1]
		else:
			left_grch = revcomplement(chromosomes[chr]
				[pos+anchor_len-full_len-1:pos+anchor_len-1])
		
		if mstrand == '+':
			right_grch = chromosomes[mchr][
				mpos+anchor_len-full_len-1:mpos+anchor_len-1]
		else:
			right_grch = revcomplement(chromosomes[mchr]
				[mpos-1:mpos+full_len-1])
		
		# If the read is at the very edge of a chromosome, ignore it.
		if len(left_grch) < full_len or len(right_grch) < full_len:
			continue
		
		#print('-------------------')
		#print([chr, strand, pos, mchr, mstrand, mpos])
		#print(seq)
		#print(left_grch)
		#print(right_grch)
			
		# Check that the read sequence is not too homologous on either side
		# of the breakpoint.
		left_match = float(sum([seq[i] == left_grch[i]
			for i in range(full_len - anchor_len, full_len)])) / anchor_len
		right_match = float(sum([seq[i] == right_grch[i]
			for i in range(anchor_len)])) / anchor_len

		max_homology = 0.7
		if left_match >= max_homology or right_match >= max_homology: continue
		
		# Identify the breakpoint location that minimizes the number of
		# nucleotide mismatches between the read and the breakpoint flanks.
		potential_breakpoints = range(anchor_len, full_len - anchor_len + 1)
		mismatches = [0] * len(potential_breakpoints)
		for k, br in enumerate(potential_breakpoints):
			grch_chimera = left_grch[:br] + right_grch[br:]
			mismatches[k] = sum([seq[i] != grch_chimera[i]
				for i in range(full_len)])
			
		# The best breakpoint placement cannot have more than N mismatches.
		least_mismatches = min(mismatches)
		if least_mismatches > 2: continue
		
		# "br" represent the number of nucleotides in the read
		# before the breakpoint, counting from the 5' end of the read.
		# If there is microhomology, we pick the first breakpoint.
		br = potential_breakpoints[mismatches.index(least_mismatches)]
		
		# Now that we know the exact fusion breakpoint, we mark mismatches
		# with a lower case nucleotide and augment the read
		# sequence with a | symbol to denote the junction.
		grch_chimera = left_grch[:br] + right_grch[br:]
		seq = ''.join([nuc if grch_chimera[k] == nuc else nuc.lower()
			for k, nuc in enumerate(seq)])
		seq = seq[:br] + '|' + seq[br:]
		
		# Make positions represent read starts.
		if strand == '-': pos += anchor_len - 1
		if mstrand == '-': mpos += anchor_len - 1
		
		# Each discordant anchor pair is represented as a 7-tuple
		# (chr_1, strand_1, pos_1, chr_2, strand_2, pos_2, sequence).
		# Positions are 1-based and represent read starts.
		out.write('%s\t%s\t%d\t%s\t%s\t%d\t%s\n' % (
			chr, strand, pos, mchr, mstrand, mpos, seq))
		N += 1
		
	info('Found %d discordant anchor pairs.' % N)
	out.close()
	




def detect_rearrangements(sam_path, genome_path, out_prefix, anchor_len,
	min_rearrangement_size, min_mapq, orientation, max_frag_len, 
	discard_pcr_duplicates=True):
	
	if not os.path.exists(sam_path):
		error('File %s does not exist.' % sam_path)
	
	detect_discordant_pairs(sam_path, out_prefix,
		min_rearrangement_size=min_rearrangement_size, min_mapq=min_mapq,
		orientation=orientation)
	
	# Execute split read analysis if the user has specified an anchor length.
	if anchor_len > 0:
		detect_discordant_mates(sam_path, genome_path, out_prefix, anchor_len,
			min_rearrangement_size=min_rearrangement_size)
	
	info('Sorting discordant pairs by chromosomal position...')
	sort_inputs = '<(gunzip -c %s.discordant_pairs.tsv.gz)' % out_prefix
	if anchor_len > 0:
		sort_inputs +=' <(gunzip -c %s.discordant_singles.tsv.gz)' % out_prefix

	sort_tmp_dir = os.path.dirname(out_prefix)
	if not sort_tmp_dir: sort_tmp_dir = './'

	shell('sort -k1,1 -k3,3n -T %s %s | gzip -c > %s.sorted_pairs.tsv.gz' %
		(sort_tmp_dir, sort_inputs, out_prefix))
	
	def print_rearrangement(out, r, discard_pcr_duplicates):
		if discard_pcr_duplicates:
			r.reads = list(set(r.reads))
			
		if len(r.reads) < 2: return 0
		
		pos = r.pos[1] if r.strand == '+' else r.pos[0]
		mpos = r.mpos[1] if r.mstrand == '+' else r.mpos[0]
		out.write('%s\t%s\t%d\t\t\t%s\t%s\t%d\t\t\t%d\t%d\t%s\n' % (
			r.chr, r.strand, pos, r.mchr, r.mstrand, mpos,
			sum([read[2] == None for read in r.reads]),
			sum([read[2] != None for read in r.reads]),
			';'.join([read[2] for read in r.reads if read[2] != None])))
		return 1

	info('Identifying rearrangements based on clusters of discordant reads...')
	
	out = open('%s.sv' % out_prefix, 'w')
	out.write(sv_file_header + '\n')
	
	N = 0
	rearrangements = []
	for line in zopen('%s.sorted_pairs.tsv.gz' % out_prefix):
		al = line[:-1].split('\t')
		
		chr = al[0]
		strand = al[1]
		pos = int(al[2])
		mchr = al[3]
		mstrand = al[4]
		mpos = int(al[5])
		seq = None if al[6] == '-' else al[6]
		
		# Rearrangements that are too far need not be considered in the future
		far = (r for r in rearrangements if pos - r.pos[1] > max_frag_len)
		for r in far:
			N += print_rearrangement(out, r, discard_pcr_duplicates)
			
		rearrangements = [r for r in rearrangements if
			pos - r.pos[1] <= max_frag_len]
		
		# Check if we already have a rearrangement that matches the new pair.
		# We don't check the distance for the first mate because we already
		# know from above the rearrangements near it. In comparing the strands,
		# we assume that fragments are not oriented.
		matches = [r for r in rearrangements if 
			point_region_distance(mpos, r.mpos) <= max_frag_len and
			chr == r.chr and mchr == r.mchr and
			strand == r.strand and mstrand == r.mstrand]
		
		read = (pos, mpos, seq)
		if matches:
			for match in matches:
				# Only the 3' boundary of the left mate can extend because
				# the discordant pairs were ordered by left mate position.
				match.pos = (match.pos[0], max(match.pos[1], pos))
				match.mpos = (min(match.mpos[0], mpos),
					max(match.mpos[1], mpos))
				match.reads.append(read)
				
		else:
			# No suitable rearrangements, create a new one.
			rearrangements.append(Rearrangement(
				chr, strand, pos, mchr, mstrand, mpos, read))
	
	for r in rearrangements:
		N += print_rearrangement(out, r, discard_pcr_duplicates)
			
	info('Found %d rearrangements with at least 2 reads of evidence.' % N)
	





#############################
# BREAKFAST DETECT SPECIFIC #
#############################

def detect_specific(bam_path, donors_path, acceptors_path, genome_path,
	out_prefix, all_reads):

	read_len = sam.read_length(bam_path)
	info('Using read length %d bp...' % read_len)

	flank_len = read_len - 10
	chromosomes = read_fasta(genome_path)

	donor_exons = regions_from_bed(donors_path)
	donors = []
	for ex in donor_exons:
		chr = ex[0] if ex[0].startswith('chr') else 'chr'+ex[0]
		chr_seq = chromosomes[chr]
		if ex[1] == '+':
			donors.append((chr, '+', ex[3], chr_seq[ex[3]-flank_len:ex[3]]))
		elif ex[1] == '-':
			donors.append((chr, '-', ex[2],
				revcomplement(chr_seq[ex[2]-1:ex[2]-1+flank_len])))

	acceptor_exons = regions_from_bed(acceptors_path)
	acceptors = []
	for ex in acceptor_exons:
		chr = ex[0] if ex[0].startswith('chr') else 'chr'+ex[0]
		chr_seq = chromosomes[chr]
		if ex[1] == '+':
			acceptors.append((chr, '+', ex[2],
				chr_seq[ex[2]-1:ex[2]-1+flank_len]))
		elif ex[1] == '-':
			acceptors.append((chr, '-', ex[3],
				revcomplement(chr_seq[ex[3]-flank_len:ex[3]])))
				
	del chromosomes    # Release 3 GB of memory
	gc.collect()

	# Remove duplicate acceptors and donors.
	acceptors = list(set(acceptors))
	donors = list(set(donors))
	
	# Calculate junction sequences
	junctions = {}
	for left in donors:
		for right in acceptors:
			name = '%s:%s:%d_%s:%s:%d' % (left[:3] + right[:3])
			junctions[name] = Object(sequence=left[3]+right[3], reads=[])
	info('Generated %d junctions.' % len(junctions))
	
	# Build Bowtie index
	info('Constructing junction FASTA file...')
	index_fasta_path = out_prefix + '_ref.fa'
	index = open(index_fasta_path, 'w')
	for name, junction in junctions.iteritems():
		index.write('>%s\n%s\n' % (name, junction.sequence))
	index.close()
	info('Constructing Bowtie index...')
	shell('bowtie-build -q %s %s_index' % (index_fasta_path, out_prefix))
	
	# Align reads against junctions and tally junction read counts.
	if all_reads:
		info('Aligning all reads against index...')
		reads_command = 'sam reads %s' % bam_path
	else:
		info('Aligning unaligned reads against index...')
		reads_command = 'sam unaligned reads %s' % bam_path

	for line in shell_stdout('bowtie -f -v1 -B1 %s_index <(%s)'
		% (out_prefix, reads_command)):
		cols = line.rstrip().split('\t')
		junctions[cols[2]].reads.append(cols[4])
	
	shell('rm %s_index.* %s_ref.fa' % (out_prefix, out_prefix))

	out_file = open(out_prefix + '.tsv', 'w')
	out_file.write('5\' breakpoint\t3\' breakpoint\tNum reads\tSequences\n')
	for name, j in junctions.iteritems():
		if not j.reads: continue
		flanks = name.split('_')
		out_file.write('%s\t%s\t%d\t' % (flanks[0], flanks[1], len(j.reads)))
		#out_file.write(';'.join(j.reads))
		out_file.write('\n')
	out_file.close()








####################
# BREAKFAST FILTER #
####################

def sv_locus_identifiers(chr, pos, resolution=5000):
	bins = int(round(float(pos) / resolution))
	bins = [x*resolution for x in range(bins-1, bins+2)]
	return ['%s:%d' % (chr, x) for x in bins]

def filter_variants(sv_path, min_reads, blacklist_path=None):
	
	read_rules = [r.split('-') for r in min_reads]
	for k, r in enumerate(read_rules):
		if len(r) != 3:
			error('Invalid minimum read rule %s specified.' % min_reads[k])
	
	blacklist = set()
	if blacklist_path:
		blacklist = set([x.rstrip('\n') for x in open(blacklist_path)])
	
	sv_file = open(sv_path)
	sys.stdout.write(next(sv_file))   # Header
	for line in sv_file:
		tokens = line.rstrip('\n').split('\t')
		
		valid = [int(tokens[10]) >= int(rule[0]) and
			int(tokens[11]) >= int(rule[1]) and
			int(tokens[10]) + int(tokens[11]) >= int(rule[2])
			for rule in read_rules]
		if not any(valid): continue

		chrom = tokens[0]
		pos = int(tokens[2])
		loci_1 = set(sv_locus_identifiers(chrom, pos))
		
		chrom = tokens[5]
		pos = int(tokens[7])
		loci_2 = set(sv_locus_identifiers(chrom, pos))

		# We discard a rearrangement if *both* endpoints are located
		# in blacklisted regions.
		if loci_1.isdisjoint(blacklist) or loci_2.isdisjoint(blacklist):
			sys.stdout.write(line)
		
	sv_file.close()
	
	
	
		
	
	
	
	
######################
# BREAKFAST ANNOTATE #
######################
	
def distance_to_gene(sv_pos, gene_pos):
	return max([0, gene_pos[0] - sv_pos, sv_pos - gene_pos[1]])



def annotate_variants(sv_path, bed_path):
	features = []
	bed_file = open(bed_path)
	for line in bed_file:
		c = line.rstrip().split('\t')
		features.append((c[0], c[5], (int(c[1]), int(c[2])), c[3]))
	
	print(sv_file_header)
	sv_file = open(sv_path)
	for line in sv_file:
		if not line.startswith('chr'): continue
		
		tokens = line[:-1].split('\t')
		chr_1 = tokens[0]
		strand_1 = tokens[1]
		pos_1 = int(tokens[2])
		chr_2 = tokens[5]
		strand_2 = tokens[6]
		pos_2 = int(tokens[7])
		
		nearby_features_1 = [(re.sub(' \(ENSG.*?\)', '', f[3]),
			distance_to_gene(pos_1, f[2]))
			for f in features if f[0] == chr_1]
		nearby_features_2 = [(re.sub(' \(ENSG.*?\)', '', f[3]),
			distance_to_gene(pos_2, f[2]))
			for f in features if f[0] == chr_2]
		
		nearby_features_1 = [f for f in nearby_features_1 if f[1] < 100000]
		nearby_features_2 = [f for f in nearby_features_2 if f[1] < 100000]
		
		nearby_features_1.sort(key=lambda x: x[1])
		nearby_features_2.sort(key=lambda x: x[1])
		
		tokens[3] = ', '.join(['%s (%d)' % f for f in nearby_features_1])
		tokens[8] = ', '.join(['%s (%d)' % f for f in nearby_features_2])

		print('%s' % '\t'.join(tokens))
		
	sv_file.close()

	
	
	
	
	
#######################
# BREAKFAST BLACKLIST #
#######################
	
def natural_sorted(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    return sorted(l, key = alphanum_key)



def generate_blacklist(sv_files, min_frequency=0):
	S = len(sv_files)
	
	sample_variants = [[] for s in range(S)]
	
	for s, sv_file in enumerate(sv_files):
		for line in zopen(sv_file):
			if not line.startswith('chr'): continue
			tokens = line.rstrip('\n').split('\t')
			chrom = tokens[0]
			pos = int(tokens[2])
			sample_variants[s] += sv_locus_identifiers(chrom, pos)
			chrom = tokens[5]
			pos = int(tokens[7])
			sample_variants[s] += sv_locus_identifiers(chrom, pos)
	
	sample_variants = [set(loci) for loci in sample_variants]
	
	# Initial blacklist consists of all structural variants found in the
	# control samples.
	blacklist = set()
	for loci in sample_variants:
		blacklist = blacklist.union(loci)
	blacklist = natural_sorted(blacklist)
	
	# For each blacklisted locus, calculate how many percentage of the
	# control samples contain structural variants involving that locus.
	frequency = [0] * len(blacklist)
	for k, bad_variant in enumerate(blacklist):
		bad_in_sample = [bad_variant in loci for loci in sample_variants]
		#print(bad_in_sample, file=sys.stderr)
		frequency[k] = float(sum(bad_in_sample)) / len(bad_in_sample)
	
	blacklist = [x for k, x in enumerate(blacklist) if
		frequency[k] >= min_frequency]
	
	for locus in blacklist:
		print(locus)


	

	




###############################
# CREATE CIRCOS VISUALIZATION #
###############################

def visualize_circos(sv_path):
	for line in open(sv_path):
		if not line.startswith('chr'): continue
		tokens = line[:-1].split('\t')
		
		chr_1 = tokens[0].replace('chr', '')
		chr_2 = tokens[5].replace('chr', '')
		pos_1 = int(tokens[2])
		pos_2 = int(tokens[7])
		
		print('%s:%d::\t%s:%d::' % (chr_1, pos_1, chr_2, pos_2))

	



#######################################
# BREAKFAST TABULATE REARRANGED GENES #
#######################################

def tabulate_rearranged_genes(sv_paths):
	evidence = defaultdict(lambda: [[0,0] for s in sv_paths])
	nearby_genes = defaultdict(list)
	
	sample_names = [re.sub('\.sv$', '', path, re.I) for path in sv_paths]
	
	for sv_path, sample_name in zip(sv_paths, sample_names):
		sv_file = open(sv_path)
		for line in sv_file:
			if not line.startswith('chr'): continue
			
			tokens = line.rstrip().split('\t')
			nearby_1 = re.finditer(r'(\S+) \((\d+)\)', tokens[3])
			nearby_2 = re.finditer(r'(\S+) \((\d+)\)', tokens[8])
			nearby_1 = [m.group(1) for m in nearby_1
				if abs(int(m.group(2))) < 20000]
			nearby_2 = [m.group(1) for m in nearby_2
				if abs(int(m.group(2))) < 20000]
			
			S = sample_names.index(sample_name)
			if nearby_1:
				reads = evidence[nearby_1[0]][S]
				reads[0] += int(tokens[10]); reads[1] += int(tokens[11])
				nearby_genes[nearby_1[0]] += nearby_1[1:]
			if nearby_2:
				reads = evidence[nearby_2[0]][S]
				reads[0] += int(tokens[10]); reads[1] += int(tokens[11])
				nearby_genes[nearby_2[0]] += nearby_2[1:]
					
		sv_file.close()

	print('Gene\tNearby genes\tTotal positive\t%s' % '\t'.join(sample_names))
	for gene, reads in sorted(evidence.iteritems(),
		key=lambda x: sum([r[0]+r[1] > 0 for r in x[1]]), reverse=True):
		sys.stdout.write('%s\t%s\t%d' % (gene,
			', '.join(set(nearby_genes[gene])),
			sum([r[0]+r[1] > 0 for r in reads])))
		for r in reads:
			if not r[0] + r[1] > 0:
				sys.stdout.write('\t')
			else:
				sys.stdout.write('\t%d+%d' % (r[0], r[1]))
		sys.stdout.write('\n')


		
		
##############################
# BREAKFAST TABULATE FUSIONS #
##############################
		
def tabulate_fusions(sv_paths):
	pair_evidence = defaultdict(lambda: [[0,0] for s in sv_paths])
	pair_nearby_genes = defaultdict(list)
	
	sample_names = [re.sub('\.sv$', '', path, re.I) for path in sv_paths]
	
	for sv_path, sample_name in zip(sv_paths, sample_names):
		sv_file = open(sv_path)
		for line in sv_file:
			if not line.startswith('chr'): continue
			
			tokens = line[:-1].split('\t')
			nearby_1 = re.finditer(r'(\S+) \((\d+)\)', tokens[3])
			nearby_2 = re.finditer(r'(\S+) \((\d+)\)', tokens[8])
			nearby_1 = [m.group(1) for m in nearby_1
				if abs(int(m.group(2))) < 20000]
			nearby_2 = [m.group(1) for m in nearby_2
				if abs(int(m.group(2))) < 20000]
			
			# Discard intragenic rearrangements and alternative splicing.
			if set(nearby_1).intersection(nearby_2): continue
			if not nearby_1 or not nearby_2: continue
			
			pair = (nearby_1[0], nearby_2[0])
			evidence = pair_evidence[pair][sample_names.index(sample_name)]
			evidence[0] += int(tokens[10])
			evidence[1] += int(tokens[11])
			pair_nearby_genes[pair] += nearby_1[1:] + nearby_2[1:]
					
		sv_file.close()

	print('Gene pair\tNearby genes\tTotal positive\t%s' %
		'\t'.join(sample_names))
	for pair, evidence in sorted(pair_evidence.iteritems(),
		key=lambda x: sum([r[0] + r[1] > 0 for r in x[1]]), reverse=True):
		sys.stdout.write('%s-%s\t%s\t%d' % (pair[0], pair[1],
			', '.join(set(pair_nearby_genes[pair])),
			sum([r[0] + r[1] > 0 for r in evidence])))
		for r in evidence:
			if not r[0] + r[1] > 0:
				sys.stdout.write('\t')
			else:
				sys.stdout.write('\t%d+%d' % (r[0], r[1]))
		sys.stdout.write('\n')








######################################
# BREAKFAST TABULATE FUSION VARIANTS #
######################################
		
def tabulate_fusion_variants(sv_paths):
	variants = defaultdict(lambda: Object({
		'evidence': [[0,0] for s in sv_paths], 'nearby_genes': []}))
	sample_names = [re.sub('\.sv$', '', path, re.I) for path in sv_paths]
	
	for sv_path, sample_name in zip(sv_paths, sample_names):
		sv_file = open(sv_path)
		for line in sv_file:
			if not line.startswith('chr'): continue
			
			tokens = line[:-1].split('\t')
			nearby_1 = re.finditer(r'(\S+) \((\d+)\)', tokens[3])
			nearby_2 = re.finditer(r'(\S+) \((\d+)\)', tokens[8])
			nearby_1 = [m.group(1) for m in nearby_1
				if abs(int(m.group(2))) < 20000]
			nearby_2 = [m.group(1) for m in nearby_2
				if abs(int(m.group(2))) < 20000]
			
			# Discard intragenic rearrangements and alternative splicing.
			if set(nearby_1).intersection(nearby_2): continue
			if not nearby_1 or not nearby_2: continue
			
			bp = tuple(tokens[0:3] + tokens[5:8])
			v = variants[bp]
			v.breakpoints = bp
			v.genes = (nearby_1[0], nearby_2[0])
			evidence = v.evidence[sample_names.index(sample_name)]
			evidence[0] += int(tokens[10])
			evidence[1] += int(tokens[11])
			v.nearby_genes += nearby_1[1:] + nearby_2[1:]
					
		sv_file.close()

	print('Gene pair\tNearby genes\t5\' breakpoint\t3\' breakpoint\t' +
		'Total positive\t%s' % '\t'.join(sample_names))
	for variant in sorted(variants.values(),
		key=lambda v: sum(r[0]+r[1] > 0 for r in v.evidence), reverse=True):
		bp = variant.breakpoints
		sys.stdout.write('%s-%s\t%s\t%s:%s:%s\t%s:%s:%s\t%d' % (
			variant.genes[0], variant.genes[1],
			', '.join(set(variant.nearby_genes)),
			bp[0], bp[1], bp[2], bp[3], bp[4], bp[5],
			sum([r[0] + r[1] > 0 for r in variant.evidence])))
		for r in variant.evidence:
			if not r[0] + r[1] > 0:
				sys.stdout.write('\t')
			else:
				sys.stdout.write('\t%d+%d' % (r[0], r[1]))
		sys.stdout.write('\n')






########################
# BREAKFAST STATISTICS #
########################

def calculate_statistics(sv_paths):
	sample_names = [re.sub('\.sv$', '', sv_path, flags=re.I)
		for sv_path in sv_paths]

	print('Sample\tRearrangement count')
	for sv_path, sample_name in zip(sv_paths, sample_names):
		sv_count = len([1 for line in open(sv_path) if line.startswith('chr')])
		print('%s\t%d' % (sample_name, sv_count))






############################
# BREAKFAST ALIGN JUNCTION #
############################
		
def align_junction(reads):
	reads = reads.strip().split(';')
	reads = zip(reads, (seq.find('|') for seq in reads))
	reads.sort(key=lambda x: x[1])
	longest = reads[-1][1]
	for read in reads:
		print('%s%s' % (' ' * (longest - read[1]), read[0]))


			
			

################################
# BREAKFAST FILTER BY DISTANCE #
################################

def filter_distance(sv_path, min_distance):
	for line in zopen(sv_path):
		if not line.startswith('chr'):
			sys.stdout.write(line)
			continue
		
		tokens = line[:-1].split('\t')
		
		if tokens[0] != tokens[5] or abs(
			int(tokens[2]) - int(tokens[7])) >= min_distance:
			sys.stdout.write(line)






##############################
# BREAKFAST FILTER BY REGION #
##############################

def filter_by_region(sv_path, region):
	m = re.match(r'(chr.+): *(\d+) *- *(\d+)', region.strip())
	if not m: error('Invalid region specified.')

	chr = m.group(1)
	start = int(m.group(2))
	end = int(m.group(3))

	for line in zopen(sv_path):
		if not line.startswith('chr'):
			sys.stdout.write(line)
			continue
		c = line.rstrip().split('\t')

		if not chr in (c[0], c[5]): continue
		if (start <= int(c[2]) <= end) or (start <= int(c[7]) <= end):
			sys.stdout.write(line)






#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	
	if args['detect'] and args['specific']:
		detect_specific(args['<bam_file>'], args['<donors>'],
			args['<acceptors>'], args['<genome>'], args['<out_prefix>'], all_reads=args['--all-reads'])
	elif args['detect']:
		detect_rearrangements(args['<bam_file>'], args['<genome>'],
			args['<out_prefix>'],
			anchor_len=int(args['--anchor-len']),
			min_mapq=int(args['--min-mapq']),
			min_rearrangement_size=int(args['--min-distance'])*1000,
			orientation=args['--orientation'],
			max_frag_len=int(args['--max-frag-len']))
	elif args['filter'] and args['distance']:
		filter_distance(args['<sv_file>'], int(args['<min_distance>'])*1000)
	elif args['filter'] and args['region']:
		filter_by_region(args['<sv_file>'], args['<region>'])
	elif args['filter']:
		filter_variants(args['<sv_file>'],
			min_reads=args['--min-reads'],
			blacklist_path=args['--blacklist'])
	elif args['annotate']:
		annotate_variants(args['<sv_file>'], args['<bed_file>'])
	elif args['blacklist']:
		generate_blacklist(args['<sv_files>'],
			min_frequency=float(args['--freq-above']))
	elif args['visualize']:
		visualize_circos(args['<sv_file>'])
	elif args['tabulate'] and args['rearranged'] and args['genes']:
		tabulate_rearranged_genes(args['<sv_files>'])
	elif args['tabulate'] and args['fusions']:
		tabulate_fusions(args['<sv_files>'])
	elif args['tabulate'] and args['fusion'] and args['variants']:
		tabulate_fusion_variants(args['<sv_files>'])
	elif args['statistics']:
		calculate_statistics(args['<sv_files>'])
	elif args['align'] and args['junction']:
		align_junction(args['<reads>'])
		


