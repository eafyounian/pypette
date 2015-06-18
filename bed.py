#!/bin/env pypy

"""
Tools for analyzing BED files.

Usage:
  bed composite <bed_file>

Options:
  -h --help         Show this screen.
"""

from __future__ import print_function
import sys, docopt, re, os
from pypette import zopen, shell, info, error, shell_stdout


#################
# BED COMPOSITE #
#################

def bed_composite(bed_path):
	features = {}
	for line in zopen(bed_path):
		if line.startswith('#'): continue
		c = line.rstrip('\n').split('\t')
		chr, start, end, name = c[0], int(c[1]), int(c[2]), c[3]

		feature = features.setdefault(name, [chr, [(start, end)]])
	 	if chr != feature[0]: error('Chromosome mismatch.')
		
	 	segments = feature[1]
	 	overlapping = [seg for seg in segments
	 		if end >= seg[0] and start <= seg[1]]
	 	disjoint = [seg for seg in segments
	 		if not (end >= seg[0] and start <= seg[1])]
	 	disjoint.append((min([start] + [seg[0] for seg in overlapping]),
	 		max([end] + [seg[1] for seg in overlapping])))
	 	feature[1] = disjoint

	for name, feature in features.iteritems():
	 	segments = feature[1]
	 	for seg in segments:
	 		print('%s\t%d\t%d\t%s' % (feature[0], seg[0], seg[1], name))
	
	






#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['composite']:
		bed_composite(args['<bed_file>'])




