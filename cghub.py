#!/bin/env pypy

"""
Download sequencing data from the Cancer Genomics Hub.

Usage:
  cghub list <cancer> <library_type> [options]
  cghub download <cancer> <library_type> [options]

Options:
  -h --help                Show this screen.
  --type=TYPE              Only download specific file types (e.g. BAM).
  --filename=REGEXP        Regular expression that filename must match.
  --filename-in=PATH       List of filenames, everything else is excluded.
  --filename-not-in=PATH   List of filenames that must be excluded.
  --genome=VERSION         Only show data that matches the specified genome.
  --most-recent            Only download the most recent versions of BAM files.

Valid sequencing library types:
  WGS, WXS, RNA-Seq, ChIP-Seq, MeDIP-Seq, Bisulfite-Seq
"""

from __future__ import print_function
import sys, subprocess, docopt, re, os, urllib2, time
from pypette import shell, info, error, shell_stdout, Object

def identify_sample_type(sample):
	if any(f.endswith('.bam') for f in sample.files):
		return 'BAM'
	else:
		return 'Unknown'

def cghub_parse(output):
	samples = []
	sample = None
	for line in output:
		#sys.stdout.write(line)
		if re.search('^\s+Analysis \d+', line):
			sample = Object()
			sample.files = []
			sample.filesizes = []
			sample.ref_genome = ''
			continue
		
		if line.strip() == '' and sample != None:
			if sample.state == 'live':
				sample.type = identify_sample_type(sample)
				samples.append(sample)
			sample = None
			continue
		
		m = re.search('^\s+(\w+)\s+: (.*)', line)
		if not m: continue

		if m.group(1) == 'filename':
			sample.files.append(m.group(2))
		elif m.group(1) == 'filesize':
			sample.filesizes.append(int(m.group(2)))
		elif m.group(1) == 'upload_date':
			sample.upload_date = m.group(2)
		elif m.group(1) == 'aliquot_id':
			sample.aliquot = m.group(2)
		elif m.group(1) == 'center_name':
			sample.center = m.group(2)
		elif m.group(1) == 'legacy_sample_id':
			sample.legacy_sample_id = m.group(2)
		elif m.group(1) == 'analysis_data_uri':
			sample.analysis_data_uri = m.group(2)
		elif m.group(1) == 'state':
			sample.state = m.group(2)
		elif m.group(1) == 'refassem_short_name':
			sample.ref_genome = m.group(2)

	return samples



def cghub_list(samples):
	for s in samples:
		print('%s\t%s\t%s\t%s' % (s.files[0], s.legacy_sample_id,
			s.ref_genome, s.center))
		#print('%s\t%s\t%s' % (s.files[0], s.filesizes[0], s.center))

	info('Found a total of %d samples.' % len(samples))
	info('Total filesize: %.1f GB.' %
		(sum(s.filesizes[0] for s in samples) / 1e9))


	

def cghub_download(samples):
	for sample in samples:
		# Don't redownload files that are already present.
		existing = {}
		for root, dirnames, filenames in os.walk('.'):
			for f in filenames:
				path = os.path.join(root, f)
				existing[f] = os.stat(path).st_size
		
		filename = sample.files[0]
		filesize = sample.filesizes[0]

		if filename in existing and existing[filename] == filesize:
			info('%s has already been downloaded...' % filename)
			continue
		
		info('Downloading %s...' % filename)
		shell('gtdownload -v -d %s -c ~/tools/genetorrent*/cghub_2016.key' % sample.analysis_data_uri)
			
		

		
	
	
if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	predicates = [
		'disease_abbr=' + args['<cancer>'].upper(),
		'library_strategy=' + args['<library_type>']
	]
	if args['--genome']:
		predicates.append('refassem_short_name=' + args['--genome'])

	output = shell_stdout('cgquery "%s"' % '&'.join(predicates))
	samples = cghub_parse(output)

	# Filter the samples if the user has provided a whitelist
	if args['--filename']:
		rx = args['--filename']
		samples = [s for s in samples if re.search(rx, s.files[0])]

	# Filter the samples if the user has provided a filename whitelist
	if args['--filename-in']:
		whitelist = [line.strip() for line in open(args['--filename-in'])]
		samples = [s for s in samples if s.files[0] in whitelist]

	# Filter the samples if the user has provided a filename blacklist
	if args['--filename-not-in']:
		blacklist = [line.strip() for line in open(args['--filename-not-in'])]
		samples = [s for s in samples if not s.files[0] in blacklist]
	
	# Only download specific file types
	if args['--type']:
		samples = [s for s in samples if s.type == args['--type']]

	# Only keep the most recent versions of BAM files
	if args['--most-recent']:
		aliquots = {}
		for s in samples:
			if s.type != 'BAM': continue
			aliquots.setdefault(s.aliquot, []).append(s)
		for asamples in aliquots.itervalues():
			if len(asamples) == 1: continue
			samples_by_date = sorted(asamples, key=lambda s: s.upload_date)
			for s in samples_by_date[:-1]: samples.remove(s)

	if args['list']:
		cghub_list(samples)
	elif args['download']:
		cghub_download(samples)
	

