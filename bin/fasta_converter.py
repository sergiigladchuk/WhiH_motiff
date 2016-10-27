#!/usr/bin/python3

import sys

with open(sys.argv[1], 'r') as fIn, open(sys.argv[2], 'w') as fOut:
	for line in fIn:
		if line[0] == '>':
			print(line, end='', file=fOut)
		else:
			print(line.rstrip(), end='', file=fOut)