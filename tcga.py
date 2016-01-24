#!/bin/env pypy

"""
Collection of software tools for the analysis of Cancer Genome Atlas datasets.

Usage:
  tcga tumor normal pairs <samples_file>
  tcga partition <samples_file> <num_partitions>
  tcga samples by patient <samples_file>

Options:
  -h --help      Show this screen.
"""

from __future__ import print_function
import docopt, sys, re, math
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
		bloods = [s for s in psamples if re.search('TCGA-..-....-10', s)]
		adj_normals = [s for s in psamples if re.search('TCGA-..-....-11', s)]
		if bloods:
			best_normal = bloods[0]
		elif adj_normals:
			best_normal = adj_normals[0]
		for tumor in tumors:
			print('%s,%s' % (tumor, best_normal))

	


def partition(samples_path, num_partitions):
	samples = [line.strip() for line in zopen(samples_path)]
	part_size = float(len(samples)) / num_partitions
	partition_ends = [int((p+1) * part_size) for p in range(num_partitions)]
	print(partition_ends)

	patient_ids = []
	num_without_pid = 0
	for s in samples:
		m = re.search('TCGA-..-....', s)
		if not m: num_without_pid += 1
		patient_ids.append(m.group(0) if m else 'zzz' + s)

	if num_without_pid:
		info('WARNING: %d sample names did not contain a TCGA patient ID.' % num_without_pid)

	samples, patient_ids = zip(*sorted(zip(samples, patient_ids),
		key=lambda x: x[1]))

	partitions = []
	for p in range(num_partitions):
		first = sum(len(p) for p in partitions)
		last = partition_ends[p] - 1
		part = [s for s in samples[first:last+1]]
		while last + 1 < len(samples) and \
			patient_ids[last+1] == patient_ids[last]:
			part.append(samples[last+1])
			last += 1
		partitions.append(part)

	for idx, part in enumerate(partitions):
		out = open('batch_%d.txt' % (idx+1), 'w')
		for s in part: out.write('%s\n' % s)
		out.close()




def samples_by_patient(samples_path):
	samples = [line.strip() for line in zopen(samples_path)]
	patients = {}
	num_without_pid = 0
	for s in samples:
		m = re.search('TCGA-..-....', s)
		if not m:
			num_without_pid += 1
			continue
		psamples = patients.setdefault(m.group(0), [])
		psamples.append(s)

	if num_without_pid:
		info('WARNING: %d sample names did not contain a TCGA patient ID.' % num_without_pid)

	for patient, psamples in patients.iteritems():
		print('Patient %s (%d samples):' % (patient, len(psamples)))
		for sample in psamples: print('- %s' % sample)



		
	
	
if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['tumor'] and args['normal']:
		tumor_normal_pairs(args['<samples_file>'])
	elif args['partition']:
		partition(args['<samples_file>'], int(args['<num_partitions>']))
	elif args['samples'] and args['by'] and args['patient']:
		samples_by_patient(args['<samples_file>'])


