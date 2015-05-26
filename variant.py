#!/bin/env pypy

"""
Tools for calling short nucleotide variants.

Usage:
  variant call <genome_fasta> <bam_files>... [-r REGION] [--ref=N:R]
      [--hetz=N:R] [--homz=N:R] [-q N] [-Q SAMPLES] [--keep-all]
  variant recall <vcf_file> [--ref=N:R] [--hetz=N:R] [--homz=N:R]
  variant somatic <vcf_file> <tumor,normal>...
  variant discard if in controls <vcf_file> <control_samples>...
  variant discard if in <N> controls <vcf_file> <control_samples>...
  variant discard by position <vcf_file> <pos_file>
  variant discard shallow <vcf_file> <min_coverage>
  variant nonsynonymous <vcf_file>
  variant discard 1000g <vcf_file>
  variant merge <vcf_files>...
  variant annotate <vcf_file>
  variant keep samples <vcf_file> <regex>
  variant discard samples <vcf_file> <regex>
  variant conservation <vcf_file>
  variant plot evidence <vcf_file>
  variant statistics <vcf_file>
  variant signature <vcf_file> <genome_fasta>
  variant top variants <vcf_file>
  variant top mutated regions <vcf_file> <region_size>
  variant list alt samples <vcf_file>
  variant heterozygous bases <vcf_file> <pos_file>
  variant heterozygous concordance <vcf_file> <pos_file> <test> <ref>
  variant allele fractions <vcf_file> <pos_file>

Options:
  -h --help         Show this screen
  -r <region>       Restrict analysis to chromosomal region
  -q N              Minimum mapping quality score [default: 10]
  -Q SAMPLES        Samples for which mapping quality is ignored [default: ]
  --ref=N:R         Minimum evidence for homozygous reference [default: 8:0.9]
  --hetz=N:R        Minimum evidence for heterozygous [default: 4:0.25]
  --homz=N:R        Minimum evidence for homozygous alt [default: 4:0.8]
  --keep-all        Show sites even if they are all homozygous reference
"""

from __future__ import print_function
import sys, subprocess, docopt, re, os, string, math, itertools
import numpy as np
from collections import defaultdict
from pypette import zopen, shell, shell_stdout, shell_stdinout
from pypette import info, error, natural_sorted, revcomplement, read_fasta


gt_symbols = ['', '0/0', '0/1', '1/1']


##############
# BENCHMARKS #
##############

#$ time samtools mpileup -B -q20 -f ~/organisms/homo_sapiens/hg19.fa -r chr1:1-5000000 19888.bam 19893.bam > /dev/null
#real	0m24.901s
#user	0m24.406s
#sys	0m0.458s

#$ time samtools mpileup -uB -q20 -f ~/organisms/homo_sapiens/hg19.fa -r chr1:1-5000000 19888.bam 19893.bam | bcftools view -vcg -p 0.1 - > /dev/null
#real	0m57.550s
#user	0m59.659s
#sys	0m0.906s

#$ time bash -c 'samtools mpileup -B -q20 -f ~/organisms/homo_sapiens/hg19.fa -r chr1:1-5000000 19888.bam 19893.bam | ~/tools/pypette/compiled/spileup > /dev/null'
#real	0m24.477s
#user	0m31.647s
#sys	0m0.659s

#$ time snv call ~/organisms/homo_sapiens/hg19.fa 19888.bam 19893.bam > /dev/null
#real	0m38.491s
#user	0m34.363s
#sys	0m0.342s




def simple_pileup(bam_paths, genome_path, min_mapq=10, min_alt_alleles=3,
	region=None):
	
	helper_dir = os.path.dirname(os.path.realpath(__file__)) + '/compiled'
	
	options = []
	if region:
		options.append('%s %s' % ('-l' if region.endswith('.bed') else '-r', region))
	
	# samtools mpileup will automatically ignore alignments flagged as
	# duplicates
	cmd = 'samtools mpileup -A -sB %s -q0 -f %s %s 2> /dev/null | %s/spileup %d %d' % (' '.join(options), genome_path,
		' '.join(bam_paths), helper_dir, min_alt_alleles, min_mapq)
	info('Pre-filtering mutations with the following command:\n%s' % cmd)
	return shell_stdout(cmd)


def call_genotypes(reads, total_reads, options):
	# 0 = unknown, 1 = ref, 2 = hetz, 3 = homz
	ref = (total_reads - reads >= options.min_ref_reads) & \
		((total_reads - reads) / total_reads >= options.min_ref_ratio)
	hetz = (reads >= options.min_hetz_reads) & \
		(reads / total_reads >= options.min_hetz_ratio)
	homz = (reads >= options.min_homz_reads) & \
		(reads / total_reads >= options.min_homz_ratio)

	return ref + 2*hetz + homz





################
# VARIANT CALL #
################

def variant_call(bam_paths, genome_path, options):
	
	if not os.path.exists(genome_path):
		error('Could not find genome FASTA file %s.' % genome_path)

	if options.region:
		for bam_path in bam_paths:
			if not os.path.exists(bam_path + '.bai'):
				error('No index found for BAM file %s.' % bam_path)
	
	samples = [os.path.basename(p).replace('.bam', '') for p in bam_paths]
	print('CHROM\tPOSITION\tREF\tALT\t%s' % '\t'.join(samples))

	ignore_mapq = [False] * len(samples)
	if options.ignore_mapq:
		for s, sample in enumerate(samples):
			if re.search(options.ignore_mapq, sample) != None:
				ignore_mapq[s] = True
				info('Ignoring mapping quality for sample %s.' % sample)
	
	for line in simple_pileup(bam_paths, genome_path,
		min_mapq=options.min_mapq,
		min_alt_alleles=(0 if options.keep_all else options.min_hetz_reads),
		region=options.region):

		tokens = line[:-1].split('\t')
		if len(tokens) < 3: error('Invalid spileup line:\n%s' % line)
		if tokens[2] == 'N': continue
		pileups = [p.split(' ') for p in tokens[3:]]

		total_reads = np.zeros(len(samples))
		allele_reads = defaultdict(lambda: np.zeros(len(samples)))

		for s, pileup in enumerate(pileups):
			if len(pileup) < 3: continue
			for a in range(0, len(pileup), 3):
				count = int(pileup[a+1]) + \
					(int(pileup[a+2]) if ignore_mapq[s] else 0)
				total_reads[s] += count
				if pileup[a] != '.': allele_reads[pileup[a]][s] = count		

		# Call genotypes for each allele.
		for alt, reads in allele_reads.iteritems():
			genotypes = call_genotypes(reads, total_reads, options)
			if not options.keep_all and all(genotypes < 2): continue
			
			gtypes = ('%s:%d:%d' % (gt_symbols[g], reads[s], total_reads[s])
				for s, g in enumerate(genotypes))

			# Reformat indels in VCF4 format
			ref = tokens[2]
			if len(alt) >= 2:
				if alt[1] == '+':    # Insertion
					alt = (ref if alt[0] == '.' else alt[0]) + alt[2:]
				elif alt[1] == '-':  # Deletion
					ref += alt[2:]
					alt = (ref[0] if alt[0] == '.' else alt[0])

			print('%s\t%s\t%s\t%s\t%s' % (tokens[0], tokens[1], ref,
				alt.upper(), '\t'.join(gtypes)))


			
			
			
##################
# VARIANT RECALL #
##################

def variant_recall(vcf_path, options):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break
	sys.stdout.write(line)

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	for line in vcf_file:
		cols = line.rstrip().split('\t')
		gt_reads = [gt.split(':')[1:] for gt in cols[sample_col:]]
		reads = np.array([float(gt[0]) for gt in gt_reads])
		total_reads = np.array([float(gt[1]) for gt in gt_reads])
		
		genotypes = call_genotypes(reads, total_reads, options)
		if all(genotypes < 2): continue
			
		gtypes = ('%s:%d:%d' % (gt_symbols[g], reads[s], total_reads[s])
			for s, g in enumerate(genotypes))

		print('%s\t%s' % ('\t'.join(cols[:sample_col]), '\t'.join(gtypes)))




	

	

####################
# VARIANT ANNOTATE #
####################

def variant_annotate(vcf_path):
	format_annovar(vcf_path, 'anno_tmp.vcf')
	shell('table_annovar.pl anno_tmp.vcf ~/tools/annovar-090513/humandb '
		'-buildver hg19 --remove --otherinfo --outfile annotated '
		'-operation g,f,f,f '
		'-protocol refGene,cosmic64,1000g2012feb_ALL,esp6500si_all')
	
	anno = open('annotated.hg19_multianno.txt')
	out = zopen('annotated.vcf.gz', 'w')
	anno.next()
	line = anno.next()
	headers = ['CHROM', 'POSITION', 'REF', 'ALT', 'FUNCTION', 'NEARBY_GENES',
		'EXONIC_FUNCTION', 'AA_CHANGE', 'COSMIC', '1000G', 'ESP6500']
	headers += line.rstrip('\n').split('\t')[12:]
	out.write('\t'.join(headers) + '\n')
	for line in anno:
		tokens = line.rstrip('\n').split('\t')
		out.write('\t'.join(tokens[0:2] + tokens[3:]))
		out.write('\n')
	out.close()

	os.remove('anno_tmp.vcf')
	os.remove('annotated.hg19_multianno.txt')
	os.remove('annotated.invalid_input')
	os.remove('annotated.refGene.invalid_input')
	
	
def format_annovar(vcf_path, out_path):
	out = open(out_path, 'w')
	
	for line in zopen(vcf_path):
		if line.startswith(('CHROM', '#')):
			headers = line[:-1].split('\t')
			headers[1] = 'START'
			headers.insert(2, 'END')
			out.write('\t'.join(headers) + '\n')
			continue
			
		cols = line[:-1].split('\t')
		cols.insert(2, cols[1])    # Add end coordinate
		ref = cols[3]; alt = cols[4]

		if len(ref) == 1 and len(alt) > 1 and ref[0] == alt[0]:
			# Simple insertion
			ref = '-'
			alt = alt[1:]
			cols[1] = str(int(cols[1])+1)
			cols[2] = cols[1]
		elif len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
			# Simple deletion
			ref = ref[1:]
			alt = '-'
			cols[1] = str(int(cols[1])+1)
			cols[2] = str(int(cols[1]) + len(ref) - 1)
		elif len(ref) > 1 or len(alt) > 1:
			# Block substitution
			cols[2] = str(int(cols[1]) + len(ref) - 1)

		cols[3] = ref; cols[4] = alt
		out.write('\t'.join(cols))
		out.write('\n')
		
	out.close()






##################
# VARIANT FILTER #
##################

def variant_filter(vcf_path, nonsynonymous, no_1000g):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break
	sys.stdout.write(line)

	headers = line[:-1].split('\t')
	
	if nonsynonymous and not 'EXONIC_FUNCTION' in headers:
		error('Cannot find exonic function column.')
	if no_1000g and not '1000G' in headers:
		error('Cannot find 1000 Genomes column.')
		
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	col_1000g = headers.index('1000G')
	col_exonic_func = headers.index('EXONIC_FUNCTION')
	
	for line in vcf_file:
		cols = line[:-1].split('\t')

		if nonsynonymous:
			if not cols[col_exonic_func].startswith(
				('nonsynonymous', 'frameshift', 'stopgain', 'stoploss')):
				continue

		if no_1000g:
			if cols[col_1000g]: continue
				
		sys.stdout.write(line)



###################
# VARIANT SOMATIC #
###################

def somatic(vcf_path, sample_pairs):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	# Convert sample pair names into index 2-tuples.
	sample_pairs = [pair.split(',') for pair in sample_pairs]
	if not all(len(pair) == 2 for pair in sample_pairs):
		info([pair for pair in sample_pairs if len(pair) != 2])
		error('Test and control samples must be in "test,control" format.')
	for pair in sample_pairs:
		if not pair[0] in samples:
			error('Test sample %s was not found in VCF file.' % pair[0])
		if not pair[1] in samples:
			error('Control sample %s was not found in VCF file.' % pair[1])
	sample_pairs = [(samples.index(pair[0]), samples.index(pair[1]))
		for pair in sample_pairs]
	
	sys.stdout.write(line)
	
	for line in vcf_file:
		cols = line.rstrip().split('\t')
		gt_cols = cols[sample_col:]
		
		genotypes = [gt_symbols.index(g[:g.find(':')]) for g in gt_cols]
		
		somatic = [genotypes[pair[0]] >= 2 and genotypes[pair[1]] == 1
			for pair in sample_pairs]
		if not any(somatic): continue
			
		sys.stdout.write(line)







##################################
# VARIANT DISCARD IF IN CONTROLS #
##################################

def discard_if_in_controls(vcf_path, control_samples, threshold):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	control = np.array([any(re.search(rx, s) for rx in control_samples)
		for s in headers[sample_col:]])
	if not any(control): error('No control samples found.')

	info('Using these %d control samples:' % sum(control))
	for s, c in zip(headers[sample_col:], control):
		if c: info('- %s' % s)

	sys.stdout.write(line)
	for line in vcf_file:
		cols = line.rstrip().split('\t')[sample_col:]
		gt = np.array([gt_symbols.index(c[:c.find(':')]) for c in cols])
		if sum(control & (gt > 1)) >= threshold: continue
		sys.stdout.write(line)







###################
# DISCARD SHALLOW #
###################

def discard_shallow(vcf_path, min_coverage):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	sys.stdout.write(line)
	for line in vcf_file:
		cols = line.rstrip().split('\t')[sample_col:]
		total_reads = sum(int(c.split(':')[2]) for c in cols)
		if float(total_reads) / len(cols) < min_coverage: continue
		sys.stdout.write(line)






###############################
# VARIANT DISCARD BY POSITION #
###############################

def variant_discard_by_position(vcf_path, pos_path):
	info('Reading list of blacklisted positions...')
	pos_file = zopen(pos_path)
	blacklist = []
	for line in pos_file:
		cols = line.rstrip().split('\t')
		if len(cols) < 2: continue
		chr = cols[0][3:] if cols[0].startswith('chr') else cols[0]
		blacklist.append(chr + ':' + cols[1])
	blacklist = set(blacklist)

	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	sys.stdout.write(line)
	for line in vcf_file:
		cols = line.rstrip().split('\t')
		chr = cols[0][3:] if cols[0].startswith('chr') else cols[0]
		if not chr + ':' + cols[1] in blacklist:
			sys.stdout.write(line)





#################
# VARIANT MERGE #
#################

def variant_merge(vcf_paths):
	sort_in, sort_out = shell_stdinout('sort -k2,2 -k3,3n -k4,4 -k5,5')
	cons_headers = []    # Consensus headers
	vcf_samples = []     # Sample names of each VCF
	for vcf_index, vcf_path in enumerate(vcf_paths):
		info('Merging VCF file %s...' % vcf_path)
		vcf = zopen(vcf_path)
		for line in vcf:
			if not line.startswith('#'): break
		headers = line.rstrip('\n').split('\t')
		gtype_col = (4 if not 'ESP6500' in headers else
			headers.index('ESP6500') + 1)
		if not cons_headers: cons_headers = headers[:gtype_col]
		if cons_headers != headers[:gtype_col]: error('Header mismatch!')
		vcf_samples.append(headers[gtype_col:])
		for line in vcf:
			sort_in.write('%d\t%s' % (vcf_index, line))
	sort_in.close()

	print('\t'.join(cons_headers + sum(vcf_samples, [])))
	vcf_sample_counts = [len(samples) for samples in vcf_samples]
	S = sum(vcf_sample_counts)
	vcf_sample_col = [sum(vcf_sample_counts[0:k])
		for k in range(len(vcf_samples))]

	prev = None
	calls = [':0:0'] * S
	for line in sort_out:
		cols = line.rstrip('\n').split('\t')
		vcf_index = int(cols[0])
		call_col = vcf_sample_col[vcf_index]
		if prev != cols[1:5]:
			if prev != None:
				print('\t'.join(prev + calls))
			prev = cols[1:5]
			calls = [':0:0'] * S
		calls[call_col:call_col+vcf_sample_counts[vcf_index]] = \
			cols[gtype_col+1:]

	print('\t'.join(prev + calls))    # Handle the last line







	




########################
# VARIANT CONSERVATION #
########################

def variant_conservation(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ALT') + 1
	samples = headers[sample_col:]
	chr = None

	Sd2 = math.ceil(len(samples) / 2.0)

	for line in vcf_file:
		cols = line.rstrip().split('\t')
		if cols[0] != chr:
			chr = cols[0]
			print('variableStep chrom=%s' % chr)

		genotypes = np.array([gt_symbols.index(gt[:gt.find(':')])
			for gt in cols[sample_col:]])
		if any(genotypes == 0): continue

		is_alt = (genotypes >= 2)
		conserved = max(sum(is_alt), sum(1 - is_alt))
		conserved = (conserved - Sd2) / (len(samples) - Sd2)  # -> [0,1]
		print('%s\t%.2f' % (cols[1], conserved))
				
				
		
#########################
# VARIANT PLOT EVIDENCE #
#########################

def variant_plot_evidence(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ALT') + 1
	samples = headers[sample_col:]
	
	#pipe = shell_stdin('gnuplot -persist')
	#pipe.write('set terminal svg\n')
	#pipe.write('set output "out.svg"\n')
	#pipe.write('plot "-" with points\n')
	
	max_total = 100
	ratio_bin = 0.05
	bins = list(np.arange(0, 1, ratio_bin)) + [1]
	hist = np.zeros((max_total+1, len(bins)))
	
	N = 0
	for line in vcf_file:
		cols = line[:-1].split('\t')
		
		gt_cols = [gt.split(':') for gt in cols[sample_col:]]
		genotypes = [(gt_symbols.index(gt[0]), int(gt[1]), int(gt[2]))
			for gt in gt_cols]
		
		for gt in genotypes:
			ratio = float(gt[1]) / gt[2]
			total = gt[2]
			if total > max_total: continue
			bin = int(ratio / ratio_bin)
			hist[total, bin] += 1
			
		#if N > 10000: break
	
	print('TOTAL\t%s' % '\t'.join([str(c) for c in bins]))
	for total in range(hist.shape[0]):
		print('%d\t%s' % (total, '\t'.join(
			[str(c) for c in list(hist[total, :])])))
			
	#pipe.close()
			
			






######################
# VARIANT STATISTICS #
######################

def variant_statistics(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	nearby_gene_col = headers.index('NEARBY_GENES') \
		if 'NEARBY_GENES' in headers else None
	
	mutations_per_sample = np.zeros(len(samples))
	mutations_per_chr = defaultdict(lambda: np.zeros(len(samples)))
	mutations_per_gene = defaultdict(lambda: np.zeros(len(samples)))
	
	for line in vcf_file:
		cols = line[:-1].split('\t')
		gtypes = [gt.split(':')[0] for gt in cols[sample_col:]]
		gtypes = np.array([gt_symbols.index(gt) for gt in gtypes])
		mutations_per_sample += (gtypes > 1)
		mutations_per_chr[cols[0]] += (gtypes > 1)

		if nearby_gene_col:
			for nearby in cols[nearby_gene_col].split(','):
				mutations_per_gene[nearby] += (gtypes > 1)
	
	print('Sample mutation counts:')
	for s, sample_name in enumerate(samples):
		print('%s: %d' % (sample_name, mutations_per_sample[s]))

	print('Mutations per chromosome:')
	chrs = natural_sorted(mutations_per_chr.keys())
	print('SAMPLE\t%s' % '\t'.join(chrs))
	for s, sample_name in enumerate(samples):
		total = sum(mutations_per_chr[chr][s] for chr in chrs)
		if total == 0: continue
		sys.stdout.write(sample_name)
		for chr in chrs:
			sys.stdout.write('\t%d (%.1f)' % (mutations_per_chr[chr][s],
				float(mutations_per_chr[chr][s]) / total * 100))
		sys.stdout.write('\n')

	print('Top mutated genes:')
	top_genes = sorted(mutations_per_gene.iteritems(),
		key=lambda x: sum(x[1] > 0), reverse=True)
	for top in top_genes[0:100]:
		mut_samples = sum(top[1] > 0)
		if mut_samples < 2: continue
		print('%s\t%d samples' % (top[0], mut_samples))
	





#####################
# VARIANT SIGNATURE #
#####################

def variant_signature(vcf_path, genome_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	chromosomes = read_fasta(genome_path)

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]

	substitutions = []
	for ref in 'CT':
		for alt in ('AGT' if ref == 'C' else 'ACG'):
			for pre in 'ACGT':
				for post in 'ACGT':
					substitutions.append(pre+ref+post+'>'+pre+alt+post)
	sub_count = np.zeros((len(substitutions), len(samples)))

	for line in vcf_file:
		cols = line[:-1].split('\t')
		if not cols[2] in 'ACGT' or not cols[3] in 'ACGT': continue
		chr = chromosomes[cols[0]]; pos = int(cols[1])
		if chr[pos-1] != cols[2]: error('Reference mismatch!')
		ref = chr[pos-2:pos+1]
		alt = ref[0] + cols[3] + ref[2]
		if ref[1] in 'AG':
			ref = revcomplement(ref)
			alt = revcomplement(alt)

		for s, gt in enumerate(cols[sample_col:]):
			if gt_symbols.index(gt.split(':')[0]) > 1:
				sub_count[substitutions.index(ref + '>' + alt), s] += 1

	print('SUBSTITUTION\t%s' % '\t'.join(samples))
	for sub in substitutions:
		sys.stdout.write(sub)
		for count in sub_count[substitutions.index(sub), :]:
			sys.stdout.write('\t%d' % count)
		sys.stdout.write('\n')








################
# TOP VARIANTS #
################

def top_variants(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		sys.stdout.write(line)
		if not line.startswith('#'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	variants = []
	for line in vcf_file:
		cols = line[:-1].split('\t')
		gtypes = [gt.split(':')[0] for gt in cols[sample_col:]]
		variants.append((line, sum(gt_symbols.index(gt) > 1 for gt in gtypes)))

	variants = sorted(variants, key=lambda x: int(x[1]), reverse=True)
	for var in variants: sys.stdout.write(var[0])






############################
# VARIANT LIST ALT SAMPLES #
############################

def list_alt_samples(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		sys.stdout.write(line)
		if not line.startswith('#'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	for line in vcf_file:
		cols = line.rstrip('\n').split('\t')
		gtypes = [gt.split(':')[0] for gt in cols[sample_col:]]
		sys.stdout.write('\t'.join(cols[:sample_col]))
		for s, gt in enumerate(gtypes):
			if gt_symbols.index(gt) > 1: sys.stdout.write('\t%s' % samples[s])
		print()







###############################
# VARIANT TOP MUTATED REGIONS #
###############################

def variant_top_mutated_regions(vcf_path, region_size):
	if region_size % 2: error('Region size must be divisible by two.')
	step = region_size / 2

	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line.rstrip('\n').split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]

	# Construct chromosome map
	chr_sizes = defaultdict(int)
	for line in vcf_file:
		cols = line.rstrip('\n').split('\t')
		chr_sizes[cols[0]] = max(chr_sizes[cols[0]], int(cols[1]))
	vcf_file.close()

	mutated = {}		# Which samples are mutated in each bin
	variant_pos = {}	# Position of variant in bin, -1 if various
	for chr in chr_sizes:
		mutated[chr] = np.zeros((chr_sizes[chr] / step + 1, len(samples)),
			dtype=np.bool)
		variant_pos[chr] = np.zeros(chr_sizes[chr] / step + 1, dtype=np.int32)

	# Reopen VCF file (might be compressed), identify columns
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	# Tally mutated samples in each region
	print('Tallying mutated samples...')
	for line in vcf_file:
		cols = line.rstrip('\n').split('\t')
		pos = int(cols[1])
		bin = (pos - 1) / step

		vpos = variant_pos[cols[0]]
		vpos[bin] = -1 if vpos[bin] > 0 and vpos[bin] != pos else pos
		if bin > 0:
			vpos[bin-1] = -1 if vpos[bin-1] > 0 and vpos[bin-1] != pos else pos

		mut = mutated[cols[0]]
		for s, gt in enumerate(cols[sample_col:]):
			if gt_symbols.index(gt.split(':')[0]) <= 1: continue
			mut[bin, s] = True
			if bin > 0: mut[bin-1, s] = True

	# Convert mutation bitmasks into counts
	print('Convert to counts...')
	for chr in mutated:
		mutated[chr] = mutated[chr].sum(axis=1)

	# Print regions in descending order starting with highest recurrence
	print('Find maximum...')
	highest = 0
	for chr in mutated:
		highest = max(highest, max(mutated[chr]))

	print('Top regions with two or more mutated sites:')
	for n in range(highest, 1, -1):
		for chr in mutated:
			mut = mutated[chr]
			vpos = variant_pos[chr]
			for bin in range(len(mut)):
				if mut[bin] != n or vpos[bin] != -1: continue
				print('%s:%d-%d\t%d samples' % (
					chr, bin*step+1, bin*step+region_size, n))














########################
# VARIANT KEEP SAMPLES #
########################

def variant_keep_samples(vcf_path, regex, discard=False):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	keep_col = [c < sample_col or (re.search(regex, sample) != None) != discard
		for c, sample in enumerate(headers)]
	
	print('\t'.join(c for c, keep in zip(headers, keep_col) if keep))
	for line in vcf_file:
		cols = line[:-1].split('\t')
		print('\t'.join(c for c, keep in zip(cols, keep_col) if keep))






##############################
# VARIANT HETEROZYGOUS BASES #
##############################

def variant_heterozygous_bases(vcf_path, kgenomes_path):
	is_snp = np.zeros(300*1000*1000, np.bool_)
	for line in zopen(kgenomes_path):
		pos = int(line[:-1].split('\t')[1])
		is_snp[pos] = True

	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	print('\t'.join(headers[0:4]))
	for line in vcf_file:
		cols = line[:-1].split('\t')
		if not is_snp[int(cols[1])]: continue

		gt_cols = cols[sample_col:]
		genotypes = [gt_symbols.index(gt[:gt.find(':')]) for gt in gt_cols]
		total_reads = [float(gt.split(':')[2]) for gt in gt_cols]
		if not any(g == 2 and r >= 15 for g, r in zip(genotypes, total_reads)):
			continue

		print('\t'.join(cols[0:4]))





############################
# VARIANT ALLELE FRACTIONS #
############################

def variant_allele_fractions(vcf_path, pos_path):
	snps = set()
	for line in zopen(pos_path):
		pos = ':'.join(line[:-1].split('\t')[0:2])
		if not pos.startswith('chr'): pos = 'chr'+pos
		snps.add(pos)

	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	sys.stdout.write(line)
	for line in vcf_file:
		cols = line[:-1].split('\t')
		if not ':'.join(cols[0:2]) in snps: continue

		reads = [gt.split(':')[1:3] for gt in cols[sample_col:]]
		total_reads = np.array([float(r[1]) for r in reads])
		reads = np.array([float(r[0]) for r in reads])
		frac = reads / total_reads
		#frac = [f if f <= 0.5 else 1.0 - f for f in frac]
		cols[sample_col:] = ['%.2f' % f for f in frac]
		print('\t'.join(cols))








####################################
# VARIANT HETEROZYGOUS CONCORDANCE #
####################################

def variant_heterozygous_concordance(vcf_path, kgenomes_path, test_rx, ref_rx):
	is_snp = np.zeros(300*1000*1000, np.bool_)
	for line in zopen(kgenomes_path):
		pos = int(line[:-1].split('\t')[1])
		is_snp[pos] = True

	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	test_col = [i for i, h in headers if
		re.search(test_rx, h) and i >= sample_col]
	ref_col = [i for i, h in headers if
		re.search(ref_rx, h) and i >= sample_col]
	if len(test_col) != 1: error('Test sample not found.')
	if len(ref_col) != 1: error('Reference sample not found.')

	total_hetz_in_ref = 0
	total_concordant = 0
	for line in vcf_file:
		cols = line[:-1].split('\t')
		if not is_snp[int(cols[1])]: continue

		test = cols[test_col]
		test_gt = gt_symbols.index(test[:test.find(':')])
		ref = cols[ref_col]
		ref_gt = gt_symbols.index(ref[:ref.find(':')])
		
		if ref_gt == 2:
			total_hetz_in_ref += 1
			total_concordant += (test_gt == 2)

	print('Concordance was %.1f%% (%d / %d).' % (float(total_concordant) / total_hetz_in_ref * 100, total_concordant, total_hetz_in_ref))








#######################
# COMMAND LINE PARSER #
#######################

class Options(object):
	def __init__(self, args):
		self.region = args['-r']

		self.min_mapq = int(args['-q'])
		self.ignore_mapq = args['-Q']
		self.keep_all = args['--keep-all']
		
		ref = args['--ref'].split(':')
		self.min_ref_reads = int(ref[0])
		self.min_ref_ratio = float(ref[1])
		
		hetz = args['--hetz'].split(':')
		self.min_hetz_reads = int(hetz[0])
		self.min_hetz_ratio = float(hetz[1])
		
		homz = args['--homz'].split(':')
		self.min_homz_reads = int(homz[0])
		self.min_homz_ratio = float(homz[1])

		

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['call']:
		options = Options(args)
		variant_call(args['<bam_files>'], args['<genome_fasta>'], options)
	elif args['recall']:
		variant_recall(args['<vcf_file>'], Options(args))
	elif args['annotate']:
		variant_annotate(args['<vcf_file>'])
	elif args['somatic']:
		somatic(args['<vcf_file>'], args['<tumor,normal>'])
	elif args['discard'] and args['in'] and args['controls']:
		discard_if_in_controls(args['<vcf_file>'], args['<control_samples>'],
			threshold=int(args['<N>']) if args['<N>'] else 1)
	elif args['discard'] and args['shallow']:
		discard_shallow(args['<vcf_file>'], float(args['<min_coverage>']))
	elif args['nonsynonymous']:
		variant_filter(args['<vcf_file>'], nonsynonymous=True, no_1000g=False)
	elif args['discard'] and args['1000g']:
		variant_filter(args['<vcf_file>'], nonsynonymous=False, no_1000g=True)
	elif args['discard'] and args['position']:
		variant_discard_by_position(args['<vcf_file>'], args['<pos_file>'])
	elif args['merge']:
		variant_merge(args['<vcf_files>'])
	elif args['keep'] and args['samples']:
		variant_keep_samples(args['<vcf_file>'], args['<regex>'])
	elif args['discard'] and args['samples']:
		variant_keep_samples(args['<vcf_file>'], args['<regex>'], True)
	elif args['plot'] and args['evidence']:
		variant_plot_evidence(args['<vcf_file>'])
	elif args['statistics']:
		variant_statistics(args['<vcf_file>'])
	elif args['signature']:
		variant_signature(args['<vcf_file>'], args['<genome_fasta>'])
	elif args['conservation']:
		variant_conservation(args['<vcf_file>'])
	elif args['top'] and args['variants']:
		top_variants(args['<vcf_file>'])
	elif args['top'] and args['mutated'] and args['regions']:
		variant_top_mutated_regions(args['<vcf_file>'],
			int(args['<region_size>']))
	elif args['list'] and args['alt'] and args['samples']:
		list_alt_samples(args['<vcf_file>'])
	elif args['heterozygous'] and args['bases']:
		variant_heterozygous_bases(args['<vcf_file>'], args['<pos_file>'])
	elif args['allele'] and args['fractions']:
		variant_allele_fractions(args['<vcf_file>'], args['<pos_file>'])
		


