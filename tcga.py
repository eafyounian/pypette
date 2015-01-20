#!/bin/env pypy

"""
Collection of software tools for the analysis of Cancer Genome Atlas datasets.

Usage:
  tcga tumor normal pairs <samples_file>

Options:
  -h --help      Show this screen.
"""

from __future__ import print_function
import docopt, sys, re
from pypette import error, info, zopen, Object





def tumor_normal_pairs(samples_path):
	samples = [line.strip() for line in zopen(samples_path)]
	patients = {}
	for s in samples:
		m = re.search('TCGA-..-....', s)
		if not m: continue
		patient = patients.setdefault(m.group(0), [])
		patient.append(s)

	for psamples in patients.values():
		tumors = [s for s in psamples if re.search('TCGA-..-....-0[12]', s)]
		normals = [s for s in psamples if re.search('TCGA-..-....-1[01]', s)]
		for tumor in tumors:
			for normal in normals:
				sys.stdout.write('%s,%s ' % (tumor, normal))
	print()

	




	
		
	
	
if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['tumor'] and args['normal']:
		tumor_normal_pairs(args['<samples_file>'])


