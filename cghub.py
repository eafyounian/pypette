#!/bin/env python

"""
Download sequencing data from the Cancer Genomics Hub.

Usage:
  cghub list <cancer> <library_type> [--samples=FILE]
  cghub download <cancer> <library_type> [--samples=FILE]

Options:
  -h --help     Show this screen.

Valid sequencing library types:
  WGS, WXS, RNA-Seq, ChIP-Seq, MeDIP-Seq, Bisulfite-Seq
"""

from __future__ import print_function
import sys, subprocess, docopt, re, os, urllib2, time
from pypette import shell, info, error

class Sample: pass

def cghub_query(query):
	output = subprocess.check_output(
		'python /data/csb/tools/cghub/bin/cgquery '
		'"disease_abbr=%s&library_strategy=%s"' %
		(query['cancer'], query['library_type']), shell=True)
	
	samples = []
	sample = None
	for line in output.split('\n'):
		if re.search('^\s+Analysis \d+', line):
			sample = Sample()
			sample.files = []
			continue
		
		if re.search('^\s*$', line) and sample != None:
			if sample.state == 'live': 
				samples.append(sample)
			sample = None
			continue
		
		m = re.search('^\s+filename\s+: (.*)', line)
		if m:
			sample.files.append((m.group(1), None))
			continue
		
		m = re.search('^\s+filesize\s+: (.*)', line)
		if m:
			prev = sample.files[-1]
			sample.files[-1] = (prev[0], int(m.group(1)))
			continue
		
		m = re.search('^\s+center_name\s+: (.*)', line)
		if m:
			sample.center = m.group(1)
			continue
		
		m = re.search('^\s+legacy_sample_id\s+: (.*)', line)
		if m:
			sample.legacy_sample_id = m.group(1)
			continue
		
		m = re.search('^\s+analysis_data_uri\s+: (.*)', line)
		if m:
			sample.analysis_data_uri = m.group(1)
			continue
		
		m = re.search('^\s+state\s+: (.*)', line)
		if m:
			sample.state = m.group(1)
			continue
	
	return samples



def cghub_list(samples):
	for s in samples:
		print('%s\t%s\t%s' % (s.files[0][0], s.legacy_sample_id, s.center))
		#print('%s\t%s\t%s' % (s.files[0][0], s.files[0][1], s.center))

	print('Found a total of %d samples.' % len(samples))


	

def cghub_download(samples):
	for sample in samples:
		# Don't redownload files that are already present.
		existing = {}
		for root, dirnames, filenames in os.walk('.'):
			for f in filenames:
				path = os.path.join(root, f)
				existing[f] = os.stat(path).st_size
		
		filename = sample.files[0][0]
		filesize = sample.files[0][1]

		if filename in existing and existing[filename] == filesize:
			info('%s has already been downloaded...' % filename)
			continue
		
		info('Downloading %s...' % filename)
		while not shell('/data/csb/tools/cghub/bin/gtdownload -v -d %s '
			'-c /data/csb/tools/cghub/share/cghub_20131029.key '
			'-C /data/csb/tools/cghub/share/GeneTorrent' % 
			sample.analysis_data_uri):
			time.sleep(60)
			
		

		
	
	
if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	query = {}
	query['cancer'] = args['<cancer>'].upper()
	query['library_type'] = args['<library_type>']
	
	samples = cghub_query(query)

	# Filter the samples if the user has provided a whitelist
	if args['--samples']:
		whitelist = [line.rstrip() for line in open(args['--samples'])]
		samples = [s for s in samples if
			any(re.search(rx, s.files[0][0]) for rx in whitelist if rx)]
	
	if args['list']:
		cghub_list(samples)
	elif args['download']:
		cghub_download(samples)
	

