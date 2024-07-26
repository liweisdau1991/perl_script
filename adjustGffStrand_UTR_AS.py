#!/usr/bin/env python
'''
**-- coding: utf-8 --**
**-- author: xinglongsheng2006@163.com
**-- usage: adjust the orientation of chromosomes in gff annotation 
     in support of UTR and splice isoforms --**
'''

import os
import sys
import re

def adjustChrOrientation(chrids, chr2size, origff, outgff):
	orifh = open(origff, 'r')
	outfh = open(outgff, 'w')
	
	while True:
		line = orifh.readline()
		if len(line) == 0:
			break
		line = line.rstrip()
		if len(line) == 0:
			continue
		if line.startswith('#'):
			outfh.write('{}\n'.format(line))
			continue
		arr = line.split('\t')
		startpos = int(arr[3])
		endpos = int(arr[4])
		strand = arr[6]
		if arr[0] in chrids:
			arr[3] = str(min(chr2size[arr[0]] - startpos + 1, chr2size[arr[0]] - endpos + 1))
			arr[4] = str(max(chr2size[arr[0]] - startpos + 1, chr2size[arr[0]] - endpos + 1))
			if strand == '+':
				arr[6] = '-'
			else:
				arr[6] = '+'
		outfh.write('{}\n'.format('\t'.join(arr)))
	
	orifh.close()	
	outfh.close()


def readGffFile(gfffile, chrids, chr2size):
	'''
	%prog: read gff3 annotation file 
	       add support for UTR and splice isoforms
	'''
	gene2tidDict = {}

	ifh = open(gfffile, 'r')
	gene_id = ''
	transcript_id = ''

	while True:
		line = ifh.readline()
		if len(line) == 0:
			break
		line = line.rstrip()
		if len(line) == 0:
			continue
		if line.startswith('#'):
			continue
		arr = line.split('\t')
		contig = arr[0]
		assembler = arr[1]
		strand = arr[6]
		if arr[2] == 'gene':
			gene_id = re.findall('ID=([^;]+)', arr[8])[0]
		elif arr[2] == 'mRNA':
			transcript_id = re.findall('ID=([^;]+)', arr[8])[0]
			if gene_id not in gene2tidDict:
				gene2tidDict[gene_id] = {}
				gene2tidDict[gene_id]['contig'] = contig
				gene2tidDict[gene_id]['assembler'] = assembler
				if contig in chrids:
					if strand == '+':
						gene2tidDict[gene_id]['strand'] = '-'
					else:
						gene2tidDict[gene_id]['strand'] = '+'
				else:
					gene2tidDict[gene_id]['strand'] = strand
				gene2tidDict[gene_id]['transcript'] = {}
			gene2tidDict[gene_id]['transcript'][transcript_id] = {}
			if contig in chrids:
				tempstart = int(arr[3])
				tempend = int(arr[4])
				gene2tidDict[gene_id]['transcript'][transcript_id]['start'] = min(chr2size[contig] - tempstart + 1, chr2size[contig] - tempend + 1)
				gene2tidDict[gene_id]['transcript'][transcript_id]['end'] = max(chr2size[contig] - tempstart + 1, chr2size[contig] - tempend + 1)
			else:	
				gene2tidDict[gene_id]['transcript'][transcript_id]['start'] = int(arr[3])
				gene2tidDict[gene_id]['transcript'][transcript_id]['end'] = int(arr[4])
			gene2tidDict[gene_id]['transcript'][transcript_id]['exon'] = []
			gene2tidDict[gene_id]['transcript'][transcript_id]['CDS'] = []
			gene2tidDict[gene_id]['transcript'][transcript_id]['5UTR'] = []
			gene2tidDict[gene_id]['transcript'][transcript_id]['3UTR'] = []
		elif arr[2] == 'exon':
			start = int(arr[3])
			end = int(arr[4])
			score = arr[5]
			temp = {}
			if contig in chrids:
				temp['start'] = min(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
				temp['end'] = max(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
			else:
				temp['start'] = start
				temp['end'] = end
			temp['score'] = score
			gene2tidDict[gene_id]['transcript'][transcript_id]['exon'].append(temp)
		elif arr[2] == 'CDS':
			start = int(arr[3])
			end = int(arr[4])
			score = arr[5]
			phase = arr[7]
			temp = {}
			if contig in chrids:
				temp['start'] = min(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
				temp['end'] = max(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
			else:
				temp['start'] = start
				temp['end'] = end
			temp['score'] = score
			temp['phase'] = phase
			gene2tidDict[gene_id]['transcript'][transcript_id]['CDS'].append(temp)
		elif arr[2] == 'five_prime_UTR':
			start = int(arr[3])
			end = int(arr[4])
			score = arr[5]
			temp = {}
			if contig in chrids:
				temp['start'] = min(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
				temp['end'] = max(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
			else:
				temp['start'] = start
				temp['end'] = end
			temp['score'] = score
			gene2tidDict[gene_id]['transcript'][transcript_id]['5UTR'].append(temp)
		elif arr[2] == 'three_prime_UTR':
			start = int(arr[3])
			end = int(arr[4])
			score = arr[5]
			temp = {}
			if contig in chrids:
				temp['start'] = min(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
				temp['end'] = max(chr2size[contig] - start + 1, chr2size[contig] - end + 1)
			else:
				temp['start'] = start
				temp['end'] = end
			temp['score'] = score
			gene2tidDict[gene_id]['transcript'][transcript_id]['3UTR'].append(temp)

	ifh.close()

	
	for gene in gene2tidDict.keys():
		gene2tidDict[gene]['start'] = 6000000000
		gene2tidDict[gene]['end'] = 0
		for tid in gene2tidDict[gene]['transcript'].keys():
			if gene2tidDict[gene]['transcript'][tid]['end'] > gene2tidDict[gene]['end']:
				gene2tidDict[gene]['end'] = gene2tidDict[gene]['transcript'][tid]['end']
			if gene2tidDict[gene]['transcript'][tid]['start'] < gene2tidDict[gene]['start']:
				gene2tidDict[gene]['start'] = gene2tidDict[gene]['transcript'][tid]['start']
	
	return gene2tidDict


def adjustGffStrand(gene2tst, outgff):
	ofh = open(outgff, 'w')
	for gene in sorted(gene2tst.keys(), key=lambda x: (re.findall('([A-Za-z]+)', gene2tst[x]['contig'])[0], int(re.findall('([0-9]+)', gene2tst[x]['contig'])[0]), int(gene2tst[x]['start']))):
		contig = gene2tst[gene]['contig']
		strand = gene2tst[gene]['strand']
		ofh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID={8};Name={9}\n'.format(contig,
                           gene2tst[gene]['assembler'], "gene", gene2tst[gene]['start'], gene2tst[gene]['end'],
                           ".", strand, ".", gene, gene))
		for tid in gene2tst[gene]['transcript'].keys():
			ofh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID={8};Name={9};Parent={10}\n'.format(contig,
                                  gene2tst[gene]['assembler'], "mRNA", gene2tst[gene]['transcript'][tid]['start'],
                                  gene2tst[gene]['transcript'][tid]['end'], ".", strand, ".", tid, tid, gene))
			if len(gene2tst[gene]['transcript'][tid]['5UTR']) != 0:
				sorted_arr_5utr = sorted(gene2tst[gene]['transcript'][tid]['5UTR'], key = lambda comp : comp['start'])
				for i in range(0, len(sorted_arr_5utr)):
					ofh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID={8}:utr5p{9};Name={8}:utr5p{9};Parent={8}\n'.format(contig,
                                          gene2tst[gene]['assembler'], "five_prime_UTR", sorted_arr_5utr[i]['start'], sorted_arr_5utr[i]['end'], sorted_arr_5utr[i]['score'],
                                          strand, ".", tid, (i+1)))

			sorted_arr = sorted(gene2tst[gene]['transcript'][tid]['exon'], key = lambda comp : comp['start'])
			for i in range(0, len(sorted_arr)):
				ofh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID={8}:exon{9};Name={8}:exon{9};Parent={8}\n'.format(contig,
                                          gene2tst[gene]['assembler'], "exon", sorted_arr[i]['start'], sorted_arr[i]['end'], sorted_arr[i]['score'],
                                          strand, ".", tid, (i+1)))
			sorted_arr2 = sorted(gene2tst[gene]['transcript'][tid]['CDS'], key = lambda comp : comp['start'])
			for i in range(0, len(sorted_arr2)):
				ofh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID=CDS{9}.{8};Name=CDS{9}.{8};Parent={8}\n'.format(contig,
                                          gene2tst[gene]['assembler'], "CDS", sorted_arr2[i]['start'], sorted_arr2[i]['end'], sorted_arr2[i]['score'],
	                                  strand, sorted_arr2[i]['phase'], tid, (i+1)))
			
			if len(gene2tst[gene]['transcript'][tid]['3UTR']) != 0:
				sorted_arr_3utr = sorted(gene2tst[gene]['transcript'][tid]['3UTR'], key = lambda comp : comp['start'])
				for i in range(0, len(sorted_arr_3utr)):
					ofh.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\tID={8}:utr3p{9};Name={8}:utr3p{9};Parent={8}\n'.format(contig,
                                          gene2tst[gene]['assembler'], "three_prime_UTR", sorted_arr_3utr[i]['start'], sorted_arr_3utr[i]['end'], sorted_arr_3utr[i]['score'],
                                          strand, ".", tid, (i+1)))

	ofh.close()


def readChrSize(infile):
	ifh = open(infile, 'r')

	chr2size = {}
	while True:
		line = ifh.readline()
		if len(line) == 0:
			break
		line = line.rstrip()
		arr = line.split('\t')
		chr2size[arr[0]] = int(arr[1])

	ifh.close()
	return chr2size


if __name__ == '__main__':
	origff = sys.argv[1]
	outgff = sys.argv[2]
	chrlist = sys.argv[3]
	chrsize = sys.argv[4]
	
	chrids = chrlist.split(',')
	chr2size = readChrSize(chrsize)
	gene2tst = readGffFile(origff, chrids, chr2size)
	adjustGffStrand(gene2tst, outgff)
