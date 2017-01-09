#!/usr/bin/python3
import sys
import re

with open(sys.argv[1], 'r') as inFile:
	
	
	files = [];
	for line in inFile:
		#skip first line
		if line != 'Files\n':
			
			#add the filename
			files.append(line.rstrip())

genes = {};

for listFileName in files:
	#parse the file to get output column
	if re.search('predetec.*',listFileName) != None:
		sufix = 'pred'
	elif re.search('chip.*',listFileName) != None:
		sufix = 'chip'
	else:
		sufix = 'fimo'
	#get prefix
	split1 = listFileName.split('_')
	prefix = split1[len(split1) - 1].split('.')[0]

	
	#loop through file

	with open(listFileName, 'r') as listFile:
		for geneLine in listFile:
			
			if (geneLine[0:5] == '"SVEN'):
				geneVal = geneLine.rstrip().split('\t')
				geneName = geneVal[0]
				geneProduct = geneVal[5]
				
				if not (geneName in genes):
					genes[geneName] = {'product': geneProduct,
										'total_all' : 0,
										'8h_all' : 0,
										'10h_all' : 0,
										'12h_all' : 0,
										'14h_all' : 0,
										'16h_all' : 0,
										'18h_all' : 0,
										'20h_all' : 0,
										'total_fimo' : 0,
										'8h_fimo' : 0,
										'10h_fimo' : 0,
										'12h_fimo' : 0,
										'14h_fimo' : 0,
										'16h_fimo' : 0,
										'18h_fimo' : 0,
										'20h_fimo' : 0,
										'total_pred' : 0,
										'8h_pred' : 0,
										'10h_pred' : 0,
										'12h_pred' : 0,
										'14h_pred' : 0,
										'16h_pred' : 0,
										'18h_pred' : 0,
										'20h_pred' : 0,
										'total_chip' : 0,
										'8h_chip' : 0,
										'10h_chip' : 0,
										'12h_chip' : 0,
										'14h_chip' : 0,
										'16h_chip' : 0,
										'18h_chip' : 0,
										'20h_chip' : 0}
				
				genes[geneName]['total_all'] += 1
				genes[geneName][prefix + '_all'] += 1
				genes[geneName]['total_' + sufix] += 1
				genes[geneName][prefix + '_' + sufix] += 1

#sort genes base on total rank and output
#header
print('gene\tproduct\tTotal_all\t8h_all\t10h_all\t12h_all\t14h_all\t16h_all\t18h_all\t20h_all\tTotal_fimo\t8h_fimo\t10h_fimo\t12h_fimo\t14h_fimo\t16h_fimo\t18h_fimo\t20h_fimo\tTotal_pred\t8h_pred\t10h_pred\t12h_pred\t14h_pred\t16h_pred\t18h_pred\t20h_pred\tTotal_chip\t8h_chip\t10h_chip\t12h_chip\t14h_chip\t16h_chip\t18h_chip\t20h_chip')

sorted_genes = sorted(genes, key=lambda k: genes[k]['total_all'], reverse=True)

for gene in sorted_genes:
	print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(gene,genes[gene]['product'],genes[gene]['total_all'],genes[gene]['8h_all'],genes[gene]['10h_all'],genes[gene]['12h_all'],genes[gene]['14h_all'],genes[gene]['16h_all'],genes[gene]['18h_all'],genes[gene]['20h_all'],genes[gene]['total_fimo'],genes[gene]['8h_fimo'],genes[gene]['10h_fimo'],genes[gene]['12h_fimo'],genes[gene]['14h_fimo'],genes[gene]['16h_fimo'],genes[gene]['18h_fimo'],genes[gene]['20h_fimo'],genes[gene]['total_pred'],genes[gene]['8h_pred'],genes[gene]['10h_pred'],genes[gene]['12h_pred'],genes[gene]['14h_pred'],genes[gene]['16h_pred'],genes[gene]['18h_pred'],genes[gene]['20h_pred'],genes[gene]['total_chip'],genes[gene]['8h_chip'],genes[gene]['10h_chip'],genes[gene]['12h_chip'],genes[gene]['14h_chip'],genes[gene]['16h_chip'],genes[gene]['18h_chip'],genes[gene]['20h_chip']))





