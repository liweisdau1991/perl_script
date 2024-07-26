#!/usr/bin/python

'''
**-- coding: utf-8 --**
**-- author: xinglongsheng@caas.cn --**
**-- usage: Processing kofamscan result to get the top hit --**
**-- change: Add e-value cutoff --**
'''

import os
import sys
import re
import gzip


### * evm.model.Chr10A.1        K13789  363.07  485.5  3.2e-145 geranylgeranyl diphosphate synthase, type II [EC:2.5.1.1 2.5.1.10 2.5.1.29]

def processKofamscanOutput(infile, outfile, cutoff):
	if infile.endswith('.gz'):
		ifh = gzip.open(infile, 'rt')
	else:
		ifh = open(infile, 'r')

	tid2result = {}
	while True:
		line = ifh.readline()
		if len(line) == 0:
			break
		line = line.rstrip()
		if line.startswith('#'):
			continue
		arr = re.split('\s+', line)
		query = arr[1]
		target = arr[2]
		evalue = float(arr[5])
		if query not in tid2result and evalue <= cutoff:
			tid2result[query] = target

	ifh.close()

	ofh = open(outfile, 'w')
	for tid in tid2result.keys():
		ofh.write('%s\t%s\n' % (tid, tid2result[tid])) 
	ofh.close()


if __name__ == '__main__':
	infile = sys.argv[1]
	outfile = sys.argv[2]
	if len(sys.argv) >= 4:
		cutoff = float(sys.argv[3])
	else:
		cutoff = 0.01
	processKofamscanOutput(infile, outfile, cutoff)
