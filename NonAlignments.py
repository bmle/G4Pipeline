# =============================================================================
# bmle
# GplexProject: NonAlignments.py
# Generates a 'non-alignment' file from a SAM file and genome annotation file
# =============================================================================

def main(alignPath, gffPath, outputPath):
	"""Generate a GFF file of non-aligned regions from a SAM alignment file.
	
	:param alignPath: the absolute path to the genome sam file (for alignments)
	:param gffPath: the absolute path to the genome annotation file (for sequence headers)
	:param outputPath: the absolute path to where the output file should be written
	:return: nothing
	"""
	import re
	from collections import OrderedDict
	from operator import itemgetter
	from GFF import loadSeqregs
	
	# Loads alignments from SAM file
	entries = []
	with open(alignPath, 'r') as alignFile:
		for line in alignFile:
			if not line.startswith('@'):
				entries.append(line.split('\t'))
	
	# Loads sequence headers from GFF file
	temp1 = [seq.split() for seq in sorted(loadSeqregs(gffPath))]
	seqs = OrderedDict([(seq[1], seq[3]) for seq in temp1])
	
	# Extracts coordinates of aligned sequences
	print('Extracting coordinates...')
	data = {k: [] for k in seqs}
	for line in entries:
		start = int(line[3])
		end = start + sum(map(int, re.findall('(\d+)(?=[MDNX=])', line[5]))) - 1
		data[line[2]].append((start, end))
	for k,v in data.items(): sorted(v, key=itemgetter(0))

	# Merges coordinates of overlapping aligned sequences
	print('Merging coordinates...')
	for k,v in data.items():
		i = 1
		while i < len(v):
			a = tuple([int(x) for x in v[i-1]])
			b = tuple([int(x) for x in v[i]])
			if not a[0] > b[1] and not a[1] < b[0]:
				data[k].remove(a)
				data[k].remove(b)
				data[k].insert(i-1, (min(a[0],b[0]), max(a[1],b[1])))
			else: i += 1
			
		# Adds filler annotations to beginning and end of sequence
		seqLen = str(int(seqs[k])+1)
		data[k].insert(0, (0, 0))
		data[k].append((seqLen, seqLen))
		
	# Inverts coordinates
	print('Inverting coordinates...')
	nalign = []
	for k,v in data.items():
		for i in range(1, len(v)):
			a = v[i-1]
			b = v[i]
			nalign.append((k, int(a[1])+1, int(b[0])-1))
	
	# Writes to GFF file
	print('Writing to file...')
	with open(outputPath, 'w') as outFile:
		outFile.write('##gff-version 3\n')
		for k,v in seqs.items(): outFile.write('##sequence-region ' + k + ' 1 ' + v + '\n')
		for i, item in enumerate(nalign):
			outFile.write(item[0] + '\tblastn\tnon-alignment\t' + str(item[1]) + '\t' + str(item[2]) + '\t.\t.\t.\tID=nal_' + str(i) + ';Name=nal_' + str(i) + ';Start=' + str(item[1]) + ';End=' + str(item[2]) + '\n')
	
	print('Finished!')
	
# =============================================================================

if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2], sys.argv[3])
	
	# For local testing purposes
	# from Paths import path
	# main(*path('NALs'))