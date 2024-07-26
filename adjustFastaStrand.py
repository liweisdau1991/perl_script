#!/usr/bin/env python
'''
**-- coding: utf-8 --**
**-- author: xinglongsheng2006@163.com
**-- usage: adjust the orientation of chromosomes in fasta sequence --**
'''

import os
import sys
import fastaReader

def adjustFastaStrand(chrids, infile, outfile):
	oriseqs = fastaReader.readFastaFile(infile)
	ofh = open(outfile, 'w')
	for header in oriseqs.keys():
		if header in chrids:
			newseq = fastaReader.revComplement(oriseqs[header])
			fastaReader.write2fh(header, newseq, ofh)
		else:
			fastaReader.write2fh(header, oriseqs[header], ofh)

	ofh.close()


if __name__ == '__main__':
	infile = sys.argv[1]
	chrlist = sys.argv[2]
	outfile = sys.argv[3]

	chrids = chrlist.split(',')
	adjustFastaStrand(chrids, infile, outfile)