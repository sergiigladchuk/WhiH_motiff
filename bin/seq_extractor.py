#!/usr/bin/python3

import argparse
import sys

usage = '''This program extracts regions from bakterial genome (input as fasta) based on csv table from ChIP-seq analysis'''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-v','--version',
		    action='version',
		    version='%(prog)s 0.1')

parser.add_argument('-n',
		    type=int,
		    metavar='LINES',
		    dest='lines',
		    help='number of lines to analyze')

parser.add_argument('-i',
		    dest='infile',
		    metavar='INFILE',
		    type=argparse.FileType('r'),
		    required=True)

parser.add_argument('-o',
		    dest='outfile',
		    metavar='OUTFILE',
		    type=argparse.FileType('w'),
		    default=sys.stdout)

args = parser.parse_args()
parser.parse_args()

#main routine


