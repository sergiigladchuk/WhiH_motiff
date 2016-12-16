#!/usr/bin/python3

import argparse
import sys
import random

usage = '''This program extracts regions from bakterial genome (input as fasta) based on csv table of positions from ChIP-seq analysis'''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-v','--version',
		    action='version',
		    version='%(prog)s 1.1')

parser.add_argument('-l','--length_seq',
		    type=int,
		    metavar='LENGTH',
		    dest='lengthSeq',
		    help='lenght of extracted sequence',
		    required=True)

parser.add_argument('-g','--genome',
		    metavar='GENOME',
		    dest='genfile',
		    type=argparse.FileType('r'),
		    help='genome in fasta format',
		    required=True)

parser.add_argument('-i',
		    dest='infile',
		    metavar='INFILE',
		    type=argparse.FileType('r'),
		    help='input file with positions in csv format from ChIP-seq')

parser.add_argument('-r','--random',
			type=int,
			metavar='RANDOM_N',
			dest='randCount',
			help='if no input file, number of random sequence produced from genome')

parser.add_argument('-o',
		    dest='outfile',
		    metavar='OUTFILE',
		    type=argparse.FileType('w'),
		    default=sys.stdout)

args = parser.parse_args()
parser.parse_args()

#input check
if not (args.randCount or args.infile):
    parser.error('No action requested, add --random or -i as input peaks for region extraction')

#main routine

#get genome of bacteria
genome = ''
for genLine in args.genfile:
	if genLine[0] != '>':
		#acces sequence line
		genome += genLine.rstrip()


#loop through each line in csv and extract positions

csvLineNum =  0
positions = []
if (args.infile != None):
	for csvLine in args.infile:
		csvLineNum += 1
		
		#positions starts from raw 4 in first column
		if csvLineNum > 3:
			position = int(csvLine.split(',')[0])
			positions.append(position)
else:
	positions = random.sample(range(args.lengthSeq, len(genome) - args.lengthSeq), args.randCount)


#build new fastfile for MEME input
for pos in positions:
	#insert sequence with the specified length
	start = pos - args.lengthSeq // 2
	end = pos + args.lengthSeq // 2
	
	#in case genome is shorter
	if end <= len(genome):
	
		#build id line
		if (args.infile):
			print('>{} position from ChIP-seq'.format(pos), file=args.outfile)
		else:
			print('>{} random position from genome'.format(pos), file=args.outfile)
		
		print(genome[start : end], file=args.outfile)


