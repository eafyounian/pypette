#!/bin/env pypy

"""
Tools for calling short nucleotide variants.

Usage:
  variant call <genome_fasta> <bam_files>... [-r REGION] [--ref=N:R]
      [--hetz=N:R] [--homz=N:R] [-q N] [-Q SAMPLES] [--keep-all]
  variant recall <vcf_file> [--ref=N:R] [--hetz=N:R] [--homz=N:R]
  variant filter [--nonsynonymous] [--no-1000g] [--min-refs=N] <vcf_file>
  variant filter with paired controls <vcf_file> <tumor,normal>...
  variant filter with interleaved controls <vcf_file>
  variant filter with controls <vcf_file> <control_samples>...
  variant filter by best evidence <vcf_file> <N:R>
  variant filter contingent <vcf_file>
  variant annotate <vcf_file>
  variant rank <vcf_file>
  variant keep samples <vcf_file> <samples>...
  variant discard samples <vcf_file> <samples>...
  variant conservation <vcf_file>
  variant plot evidence <vcf_file>
  variant statistics <vcf_file>
  variant heterozygous bases <vcf_file> <pos_file>
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

  --min-refs=N      Discard site if less than N are reference [default: 0]
  --nonsynonymous   Only keep non-synonymous mutations
  --no-1000g        Discard mutations found in 1000 Genomes
"""

from __future__ import print_function
import sys, subprocess, docopt, re, os, string
import numpypy as np
from collections import defaultdict
from pypette import zopen, shell, shell_stdin, shell_stdout, argsort, temp_dir
from pypette import info, error


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
	if re.match(r'.+\.bed', region, re.I):
		options.append('-l %s' % region)
	else:
		options.append('-r %s' % region)
		
	return shell_stdout(
		'samtools mpileup -sB %s -q0 -f %s %s 2> /dev/null | %s/spileup %d %d'
		% (' '.join(options), genome_path, ' '.join(bam_paths),
		helper_dir, min_alt_alleles, min_mapq))


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
		for allele, reads in allele_reads.iteritems():
			genotypes = call_genotypes(reads, total_reads, options)
			if not options.keep_all and all(genotypes < 2): continue
			
			gtypes = ('%s:%d:%d' % (gt_symbols[g], reads[s], total_reads[s])
				for s, g in enumerate(genotypes))

			print('%s\t%s\t%s\t%s\t%s' % (tokens[0], tokens[1], tokens[2],
				allele.upper(), '\t'.join(gtypes)))


			
			
			
##################
# VARIANT RECALL #
##################

def variant_recall(vcf_path, options):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break
	sys.stdout.write(line)

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	for line in vcf_file:
		cols = line[:-1].split('\t')
		gt_reads = [gt.split(':')[1:] for gt in cols[sample_col:]]
		reads = np.array([float(gt[0]) for gt in gt_reads])
		total_reads = np.array([float(gt[1]) for gt in gt_reads])
		
		genotypes = call_genotypes(reads, total_reads, options)
		if all(genotypes < 2): continue
			
		gtypes = ('%s:%d:%d' % (gt_symbols[g], reads[s], total_reads[s])
			for s, g in enumerate(genotypes))

		print('%s\t%s' % ('\t'.join(cols[:sample_col]), '\t'.join(gtypes)))







######################
# VARIANT BACKGROUND #
######################

def variant_background(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break
	sys.stdout.write(line)

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	for line in vcf_file:
		cols = line[:-1].split('\t')[sample_col:]
		genotypes = np.array([gt_symbols.index(g[:g.find(':')]) for g in cols])
		gt_reads = [gt.split(':')[1:] for gt in cols]
		reads = np.array([float(gt[0]) for gt in gt_reads])
		total_reads = np.array([float(gt[1]) for gt in gt_reads])
		
		nonmut = (genotypes < 2) & (total_reads >= 8) 
		nonmut_frac = sum(reads[genotypes < 2]) / \
			sum(total_reads[genotypes < 2])

		line = line[:-1]
		line += '\t%.1f' % (nonmut_frac*100)
		print(line)



	



	

####################
# VARIANT ANNOTATE #
####################

def variant_annotate(vcf_path):
	format_annovar(vcf_path, 'anno_tmp.vcf')
	shell('table_annovar.pl anno_tmp.vcf '
		'/data/csb/tools/annovar-090513/humandb '
		'-buildver hg19 --remove --otherinfo --outfile annotated '
		'-operation g,f,f,f,f '
		'-protocol refGene,cosmic64,snp137,1000g2012feb_ALL,esp6500si_all')
	
	anno = open('annotated.hg19_multianno.txt')
	out = zopen('annotated.vcf.gz', 'w')
	anno.next()
	line = anno.next()
	headers = ['CHROM', 'POSITION', 'REF', 'ALT', 'FUNCTION', 'NEARBY_GENES',
		'EXONIC_FUNCTION', 'AA_CHANGE', 'COSMIC', 'dbSNP', '1000G', 'ESP6500']
	headers += line[:-1].split('\t')[13:]
	out.write('\t'.join(headers) + '\n')
	for line in anno:
		tokens = line[:-1].split('\t')
		out.write('\t'.join(tokens[0:2] + tokens[3:]))
		out.write('\n')
	out.close()
	
	
def format_annovar(vcf_path, out_path):
	out = open(out_path, 'w')
	
	for line in zopen(vcf_path):
		if line.startswith(('CHROM', '#')):
			headers = line[:-1].split('\t')
			headers[1] = 'START'
			headers.insert(2, 'END')
			out.write('\t'.join(headers) + '\n')
			continue
			
		tokens = line[:-1].split('\t')
		tokens.insert(2, tokens[1])
		
		if tokens[4][0] == '+':
			# Insertion in samtools mpileup format
			tokens[3] = '-'
			tokens[4] = tokens[4][1:]
			tokens[1] = str(int(tokens[1]) + 1)
			tokens[2] = str(int(tokens[2]) + 1)
		elif tokens[4][0] == '-':
			# Deletion in samtools mpileup format
			tokens[3] = tokens[4][1:]
			tokens[4] = '-'
			tokens[1] = str(int(tokens[1]) + 1)
			tokens[2] = str(int(tokens[2]) + len(tokens[3]))
			
		out.write('\t'.join(tokens))
		out.write('\n')
		
	out.close()






##################
# VARIANT FILTER #
##################

def variant_filter(vcf_path, nonsynonymous, no_1000g, min_refs):
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
		
		if min_refs > 0:
			gtypes = [gt.split(':')[0] for gt in cols[sample_col:]]
			if sum(gt_symbols.index(gt) == 1 for gt in gtypes) < min_refs:
				continue

		if nonsynonymous:
			if not cols[col_exonic_func].startswith('nonsy'): continue
			
		if no_1000g:
			if cols[col_1000g]: continue
				
		sys.stdout.write(line)



#######################################
# VARIANT FILTER WITH PAIRED CONTROLS #
#######################################

def variant_filter_with_paired_controls(vcf_path, sample_pairs):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
	samples = headers[sample_col:]
	
	# Convert sample pair names into index 2-tuples.
	sample_pairs = [pair.split(',') for pair in sample_pairs]
	if not all(len(pair) == 2 for pair in sample_pairs):
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
		cols = line[:-1].split('\t')
		gt_cols = cols[sample_col:]
		
		genotypes = [gt_symbols.index(g[:g.find(':')]) for g in gt_cols]
		
		somatic = [genotypes[pair[0]] >= 2 and genotypes[pair[1]] == 1
			for pair in sample_pairs]
		if not any(somatic): continue
			
		sys.stdout.write(line)




############################################
# VARIANT FILTER WITH INTERLEAVED CONTROLS #
############################################

def variant_filter_with_interleaved_controls(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1
		
	sys.stdout.write(line)
	for line in vcf_file:
		cols = line.rstrip().split('\t')
		gt = [gt_symbols.index(g[:g.find(':')]) for g in cols[sample_col:]]
		if not any(gt[s] >= 2 and gt[s+1] == 1 for s in range(0, len(gt), 2)):
			continue
		sys.stdout.write(line)




################################
# VARIANT FILTER WITH CONTROLS #
################################

def variant_filter_with_controls(vcf_path, control_samples):
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
		#gt_reads = [c.split(':')[1:] for c in cols]
		#alt_reads = np.array([float(gt[0]) for gt in gt_reads])
		#total_reads = np.array([float(gt[1]) for gt in gt_reads])
		if any(control & (gt > 1)) or not any(control & (gt == 1)):  
			continue
		sys.stdout.write(line)







###################################
# VARIANT FILTER BY BEST EVIDENCE #
###################################

def variant_filter_by_best_evidence(vcf_path, min_best):
	min_reads = int(min_best.split(':')[0])
	min_ratio = float(min_best.split(':')[1])

	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line.rstrip().split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	sys.stdout.write(line)
	for line in vcf_file:
		cols = line.rstrip().split('\t')[sample_col:]
		#gt = np.array([gt_symbols.index(c[:c.find(':')]) for c in cols])
		gt_reads = [c.split(':')[1:] for c in cols]
		alt_reads = np.array([float(gt[0]) for gt in gt_reads])
		total_reads = np.array([float(gt[1]) for gt in gt_reads])
		ratio = alt_reads / total_reads
		if not any((alt_reads >= min_reads) & (ratio >= min_ratio)):
			continue
		sys.stdout.write(line)








################
# VARIANT RANK #
################

def variant_rank(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	sys.stdout.write(line)

	lines = []
	for line in vcf_file:
		cols = line.rstrip().split('\t')[sample_col:]
		genotypes = [gt_symbols.index(gt[:gt.find(':')]) for gt in cols]
		lines.append((sum(gt >= 2 for gt in genotypes), line))

	lines.sort(key=lambda x: x[0], reverse=True)
	for line in lines:
		sys.stdout.write(line[1])





########################
# VARIANT CONSERVATION #
########################

def variant_conservation(vcf_path):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('##'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ALT') + 1
	samples = headers[sample_col:]
	
	chr = ''
	wsize = 500e3

	for line in vcf_file:
		cols = line[:-1].split('\t')
		genotypes = [gt[:gt.find(':')] for gt in cols[sample_col:]]
		genotypes = np.array([gt_symbols.index(gt) for gt in genotypes])
		
		if cols[0] != chr:
			chr = cols[0]
			window = (1, wsize)
		
		while int(cols[1]) > window[1]:
			window = (window[0] + wsize, window[1] + wsize)
			if window_genotypes:
				print('')
				window_genotypes = []
				
				
		
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
	
	col_gene = headers.index('NEARBY_GENES')
		
	gene_mutations = defaultdict(int)
	
	for line in vcf_file:
		cols = line[:-1].split('\t')
		gtypes = [gt.split(':')[0] for gt in cols[sample_col:]]
		
		gene_mutations[cols[col_gene]] += sum(
			gt_symbols.index(gt) > 1 for gt in gtypes)
	
	top_genes = sorted(gene_mutations.iteritems(), key=lambda x: x[1],
		reverse=True)
	for top in top_genes:
		print('%s\t%d' % (top[0], top[1]))
	



########################
# VARIANT KEEP SAMPLES #
########################

def variant_keep_samples(vcf_path, samples, discard=False):
	vcf_file = zopen(vcf_path)
	for line in vcf_file:
		if not line.startswith('#'): break

	headers = line[:-1].split('\t')
	sample_col = headers.index('ESP6500' if 'ESP6500' in headers else 'ALT')+1

	keep_col = [c < sample_col or (sample in samples) != discard
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
	elif args['filter'] and args['with'] and args['interleaved'] and args['controls']:
		variant_filter_with_interleaved_controls(args['<vcf_file>'])
	elif args['filter'] and args['with'] and args['paired'] and args['controls']:
		variant_filter_with_paired_controls(args['<vcf_file>'], args['<tumor,normal>'])
	elif args['filter'] and args['with'] and args['controls']:
		variant_filter_with_controls(args['<vcf_file>'],
			args['<control_samples>'])
	elif args['filter'] and args['by'] and args['best'] and args['evidence']:
		variant_filter_by_best_evidence(args['<vcf_file>'], args['<N:R>'])
	elif args['filter']:
		variant_filter(args['<vcf_file>'],
			nonsynonymous=args['--nonsynonymous'], no_1000g=args['--no-1000g'],
			min_refs=int(args['--min-refs']))
	elif args['rank']:
		variant_rank(args['<vcf_file>'])
	elif args['keep'] and args['samples']:
		variant_keep_samples(args['<vcf_file>'], args['<samples>'])
	elif args['discard'] and args['samples']:
		variant_keep_samples(args['<vcf_file>'], args['<samples>'], True)
	elif args['plot'] and args['evidence']:
		variant_plot_evidence(args['<vcf_file>'])
	elif args['statistics']:
		variant_statistics(args['<vcf_file>'])
	elif args['heterozygous'] and args['bases']:
		variant_heterozygous_bases(args['<vcf_file>'], args['<pos_file>'])
	elif args['allele'] and args['fractions']:
		variant_allele_fractions(args['<vcf_file>'], args['<pos_file>'])
		


