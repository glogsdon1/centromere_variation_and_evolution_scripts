#!/usr/bin/env python
import re
import argparse # for writing user-friendly command-line interfaces
import sys
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--input", help="input cigar.txt",  type=argparse.FileType('r'), default="/net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/ape_cens/tandemaligner/tandemAligner/chr5/cigar.txt")
parser.add_argument("-r", "--ref", help="chromosome you aligned", type=str, default=1)
parser.add_argument("-q", "--query", help="chromosome you aligned", type=str, default=1)
parser.add_argument("-o", "--out", help="output tsv", type=argparse.FileType('w'), default=sys.stdout)
args = parser.parse_args()

def get_change(cigar_change, ref_pos, bin_boundary, ref_pos_change):
	if (ref_pos - bin_boundary) < ref_pos_change:
		cigar_change = cigar_change+(ref_pos - bin_boundary)
	else:
		cigar_change = cigar_change+ref_pos_change
	return(cigar_change)

def len_ref(cigar):
	cigar = cigar
	q_start = 0
	r_start = 0
	pattern = re.compile('([MIDNSHPX=])')
	values = pattern.split(cigar)[:-1]
	paired = (values[n:n+2] for n in range(0, len(values), 2))
	i = q_start
	g = r_start
	#print(bin_start,bin_stop)
	for pair in paired:
		l = int(pair[0])
		t = pair[1]
		
		if t == 'M': ## if match, return consecutive coordinates
			i += l
			g += l
			#if g >= bin_start and g <= bin_stop:
			#	matches = get_change(matches, g, bin_start, l)
			#elif g >= bin_stop:
			#	matches = matches+bin_stop-(g-l)

		elif t == 'N': ## skipped region from the reference
			g += l
		elif t == 'D':
			g += l
			#if g >= bin_start and g <= bin_stop:
			#	deletions = deletions+1	
			#elif g >= bin_stop:
			#	deletions = deletions+1
		elif t == 'I':
			i += l
			#if g >= bin_start:
			#	insertions = insertions+1
		elif t == 'X':
			i += l
			g += l
			#if g >= bin_start and g <= bin_stop:
			#	mismatches = get_change(mismatches, g, bin_start, l)
			#elif g >= bin_stop:
			#	mismatches = matches+bin_stop-(g-l)
			#	cigar_keep = cigar_keep+str(bin_stop-(g-l))+t
		elif t == '=':
			i += l
			g += l
			#if g >= bin_start and g <= bin_stop:
			#	matches = get_change(matches, g, bin_start, l)
			#elif g >= bin_stop:
			#	matches = matches+bin_stop-(g-l)
		elif t == 'S': ## soft clipping (clipped sequences present in SEQ)
			i += l
			pass
		elif t == 'H': ## hard clipping (clipped sequences NOT present in SEQ)
			pass
		elif t == 'P': ## padding (silent deletion from padded reference)
			pass
		else:
			print("unknown symbol")
	return(g)

def split_cigar(cigar, bin_start, bin_stop):
	cigar = cigar
	q_start = 0
	r_start = 0
	bin_start = bin_start
	bin_stop = bin_stop
	cigar_keep=""
	g_keep=[]
	matches = 0
	mismatches = 0
	deletions = 0
	insertions = 0
	pattern = re.compile('([MIDNSHPX=])')
	values = pattern.split(cigar)[:-1]
	paired = (values[n:n+2] for n in range(0, len(values), 2))
	i = q_start
	g = r_start
	isSet = False
	query_start = None
	l = None
	#print(bin_start,bin_stop)
	for pair in paired:

		if g >= bin_start and not isSet:
			if l is None:
				query_start = 0
			elif l <= (bin_stop-bin_start):
				query_start = i-l
			else:
				query_start = i-(bin_stop-bin_start)
			query_start = i
			isSet=True
		if g >= bin_stop:
			perID_events = matches/(matches+mismatches+insertions+deletions)
			if query_start is None:
				if l <= (bin_stop-bin_start):
					query_start = i-l
				else:
					query_start = i-(bin_stop-bin_start)
			return({'perID_events': perID_events,'query_start':query_start,'query_end':i})


		l = int(pair[0])
		t = pair[1]
		if t == 'M': ## if match, return consecutive coordinates
			i += l
			g += l
			if g >= bin_start and g <= bin_stop:
				matches = get_change(matches, g, bin_start, l)
			elif g >= bin_stop:
				matches = matches+bin_stop-(g-l)

		elif t == 'N': ## skipped region from the reference
			g += l
		elif t == 'D':
			g += l
			if g >= bin_start and g <= bin_stop:
				deletions = deletions+1	
			elif g >= bin_stop:
				deletions = deletions+1
		elif t == 'I':
			i += l
			if g >= bin_start:
				insertions = insertions+1
		elif t == 'X':
			i += l
			g += l
			if g >= bin_start and g <= bin_stop:
				mismatches = get_change(mismatches, g, bin_start, l)
			elif g >= bin_stop:
				mismatches = matches+bin_stop-(g-l)
				cigar_keep = cigar_keep+str(bin_stop-(g-l))+t
		elif t == '=':
			i += l
			g += l
			if g >= bin_start and g <= bin_stop:
				matches = get_change(matches, g, bin_start, l)
			elif g >= bin_stop:
				matches = matches+bin_stop-(g-l)
		elif t == 'S': ## soft clipping (clipped sequences present in SEQ)
			i += l
			pass
		elif t == 'H': ## hard clipping (clipped sequences NOT present in SEQ)
			pass
		elif t == 'P': ## padding (silent deletion from padded reference)
			pass
		else:
			print("unknown symbol")

cigar_df = pd.read_csv(args.input, header=None, sep="\t", names=["cigar"])
cigar = cigar_df.at[0,"cigar"]
#print(cigar)
perID_df = pd.DataFrame(columns=['ref','ref_start','ref_end','ref_length','query','query_start','query_end','query_length','perID_events'])
window_size=(int(int(10000)/int(10000)))+1
print(len_ref(cigar))
seq = np.arange(int(0), int(len_ref(cigar))+int(10000), int(10000))

for i in range(len(seq) - window_size + 1):
	window = seq[i: i + window_size]
	keep_dict = split_cigar(cigar, window[0],window[-1])
	if keep_dict:
		keep_dict['ref'] = args.ref
		keep_dict['query'] = args.query
		keep_dict['ref_start'] = window[0]
		keep_dict['ref_end'] = window[-1]
		perID_df = perID_df.append(keep_dict, ignore_index=True)
perID_df['ref_length'] = perID_df['ref_end'] - perID_df['ref_start']
perID_df['query_length'] = perID_df['query_end'] - perID_df['query_start']

perID_df.to_csv(args.out, header=True, sep="\t", index=False)

