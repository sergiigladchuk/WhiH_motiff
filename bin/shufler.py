#!/usr/bin/python3
import random
import sys

with open(sys.argv[1], 'r') as inFile:
	
	for line in inFile:
		if line[0] == '>':
			print(line, end='')
		else:
			line = ''.join(random.sample(line[:-1],len(line[:-1])))
			print(line)