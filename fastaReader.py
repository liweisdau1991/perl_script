#!/usr/bin/env python
'''
**-- coding: utf-8 --**
**-- author: xinglongsheng@caas.cn --**
**-- usage: read FASTA sequence --**
'''

import re
import sys
import gzip


###read in fasta file
def readFastaFile_alt(infile):
	'''read in multiple fasta sequences into a dictionary'''

	ifh = open(infile, 'r')

	seq = {}
	head = ''
	while True:
		line = ifh.readline()
		#line = line.strip()
		if len(line) == 0:
			break
		line = line.rstrip()
		matchObj = re.match(r'^>(\S+)',line,re.M|re.I)	
		matchObj2 = re.match(r'^(\w+)',line,re.M|re.I)
		if matchObj:
			head = matchObj.group(1)
			seq[head] = ''
		if matchObj2:
			seq[head] += matchObj2.group(1)

	ifh.close()
	return seq


def readFastaFile(infile):
	if infile.endswith('.gz'):
		ifh = gzip.open(infile, 'rb')
	else:
		ifh = open(infile, 'r')
	seq = {}
	head = ''
	temp = []
	while True:
		line = ifh.readline()
		if len(line) == 0:
			break
		line = line.rstrip()
		if line.startswith('>'):
			#head = line[1:].split(' ')[0]
			if head and temp:
				seq[head] = ''.join(temp)
			head = line[1:].split(' ')[0]
			temp = []
		elif len(line) > 0 and line[0].isalpha():
			temp.append(line)		

	ifh.close()
	seq[head] = ''.join(temp)

	return seq

### reverse complement ###
def revComplement(seq):
	'''obtain the reverse complement of the original sequence using map function'''
	resseq = ''
	for i in range(0,len(seq)):
		resseq += seq[len(seq) - i - 1]
	def complement(temp):
		return {'A':'T','T':'A','C':'G','G':'C','N':'N'}[temp]
	
	return ''.join(map(complement,resseq))


### write to file handler
def write2fh(head,seq,fh):
	fh.write('>%s\n' % head)
	i = 0
	while True:
		if i >= len(seq):
			break
		fh.write('%s\n' % seq[i:i+60])
		i += 60


if __name__ == '__main__':
	testSeq = 'ATCGGCTAATCGGCTA'
	revSeq = revComplement(testSeq)
	print('%s\t%s' % (testSeq,revSeq))
