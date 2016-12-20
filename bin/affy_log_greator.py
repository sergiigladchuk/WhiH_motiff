#!/usr/bin/python3

import argparse
import sys

usage = '''This program converts transcriptomic data, gff file with gene description and list of positions into gene table with affLogs and flags for further statistical analysis'''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-v','--version',
		    action='version',
		    version='%(prog)s 0.1')

parser.add_argument('-s','--start_cut',
		    type=int,
		    metavar='START_CUT',
		    dest='startCut',
		    help='cut-off start before gene start codon',
		    default=300)

parser.add_argument('-e','--end_cut',
		    type=int,
		    metavar='END_CUT',
		    dest='endCut',
		    help='cut-off end after gede start codon',
		    default=50)

parser.add_argument('-g','--gff',
		    metavar='GFF3',
		    dest='gffFile',
		    type=argparse.FileType('r'),
		    help='gene descr file in gff3 format',
		    required=True)

parser.add_argument('-i',
		    dest='infile',
		    metavar='INFILE',
		    type=argparse.FileType('r'),
		    help='input file with transcriptomics in csv format',
		    required=True)

parser.add_argument('-p',
		    dest='positions',
		    metavar='POSITIONS',
		    type=argparse.FileType('r'),
		    help='positions file for gene matching',
		    required=True)

parser.add_argument('-o',
		    dest='outfile',
		    metavar='OUTFILE',
		    type=argparse.FileType('w'),
		    default=sys.stdout)

args = parser.parse_args()
parser.parse_args()


#main routine

#create header for output
print('gene tag\tgene_old_tag\tgene_start\tgene_end\tstrand\tproduct\tclosest_match\tmatch_index\tmatch_value\tdistance\taffyLog_8h\taffyLog_10h\taffyLog_12h\taffyLog_14h\taffyLog_16h\taffyLog_18h\taffyLog_20h',file=args.outfile)

#convert FIMO to list
fimoList = []
fimoRank = 0
for fimoLine in args.positions:
	
	if fimoLine[0] == '1':
		fimoRank += 1 
		fimoValues = fimoLine.rstrip().split('\t')
		motiffMid = int((int(fimoValues[3])+int(fimoValues[2]))/2)
		fimoList.append({'rank' : fimoRank, 
						 'position' : motiffMid,
						 'score' : fimoValues[5],
						 'matched' : False})

#convert transcriptomic to dictionary
transDic = {}

indexColWT = [[94,95,96],[77,78,79],[80,81,82],[83,84,85],[86,87,88],[89,90,92],[93,94]]
indexColWhiH = [[158,159,160],[137,138,139],[140,141,142],[143,144,145],[146,147,148],[149,150,151],[152,153,154]]

def getAverage(rowList,indexes):
	outList = []
	for indGroup in indexes:
		outList.append(sum([float(oneVal) for oneVal in [rowList[x] for x in indGroup]])/len(indGroup))
	return outList

for transLine in args.infile:
	# check if data line
	if transLine[0:4] == 'SMD0' or transLine[0:4] == 'SMD1':
		transValues = transLine.rstrip().split(',')
		#record the average differances between log scores for times
		transDic[transValues[0]] = [WhiHlog - WTlog for WhiHlog, WTlog in zip(getAverage(transValues,indexColWhiH), getAverage(transValues,indexColWT))]
			
		
#loop through gff anotation file

for gffLine in args.gffFile:
	if gffLine[0] != '#' and len(gffLine) > 10:
		gffValues = gffLine.rstrip().split('\t')
		
		#acces gene line for old tag
		if gffValues[2] == 'gene':
			texts = dict(zip([fullEntry.split('=')[0] for fullEntry in gffValues[8].split(';')],[fullEntry.split('=')[1] for fullEntry in gffValues[8].split(';')]))
			tag = texts['locus_tag']
			if 'old_locus_tag' in texts:
				old_tag = texts['old_locus_tag']
			else:
				old_tag = '-'
			
		#access CDS
		if gffValues[2] == 'CDS':
			start = int(gffValues[3])
			end = int(gffValues[4])
			strand = gffValues[6]
			texts = dict(zip([fullEntry.split('=')[0] for fullEntry in gffValues[8].split(';')],[fullEntry.split('=')[1] for fullEntry in gffValues[8].split(';')]))
			
			product = texts['product']
			
			#get the closest match if any from fimo
			bestDist = args.startCut
			bestFimo = None
			
				
			for fimo in fimoList:
				
				if strand == '+':
					curDist = start - fimo['position']
				else:
					curDist = fimo['position'] - end
				#check if beter then previous
				
				if curDist < bestDist and curDist > 0:
					bestDist = curDist
					bestFimo = fimo
				elif abs(curDist) < args.endCut and abs(curDist) < bestDist:
					bestDist = abs(curDist)
					bestFimo = fimo
				
			if bestFimo != None:
				match = 1 
				matchIndex = bestFimo['rank']
				matchVal = bestFimo['score']
				matchDis = bestDist
				
			else:
				match = 0
				matchIndex = ''
				matchVal = ''
				matchDis = ''
			
			#add log values string
			if old_tag in transDic:
				logDifStr = '\t'.join(['{}'.format(x) for x in transDic[old_tag]])
			else:
				logDifStr = '\t\t\t\t\t\t'
			
			#output of the line
			print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(tag, old_tag, start, end, strand, product, match, matchIndex, matchVal, matchDis, logDifStr),file=args.outfile)
			
