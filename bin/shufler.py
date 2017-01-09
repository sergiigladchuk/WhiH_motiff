#!/usr/bin/python3
import sys

with open(sys.argv[1], 'r') as inFile:
	
	genes = {};
	
	for line in inFile:
		#scip first line
		if line == 'Files':
			print(line, end='')
		else:
			line = ''.join(random.sample(line[:-1],len(line[:-1])))
			print(line)