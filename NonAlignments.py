# =============================================================================
# bmle
# GplexProject
# Generates a non-alignment file from a SAM file
# =============================================================================

def main(alignPath, outputPath):
	"""Generate a GFF file of non-aligned regions from a SAM alignment file.
	
	:param alignPath: the absolute path to the genome sam file
	:param outputPath: the absolute path to where the output file should be written
	:return: nothing
	"""
	import re
	from operator import itemgetter
	
	# Loads SAM file
	entries = []
	seqs = []
	with open(alignPath, 'r') as alignFile:
		for line in alignFile:
			temp = line.split('\t')
			if not temp[0].startswith('@'): entries.append(temp)
			elif temp[0].startswith('@SQ'): seqs.append((temp[1][3:], temp[2][3:].strip()))
	
	# Extracts coordinates of aligned sequences
	print('Extracting coordinates...')
	data = [[] for _ in range(len(seqs))]
	for line in entries:
		start = int(line[3])
		end = start + sum(map(int, re.findall('(\d+)(?=[MDNX=])', line[5]))) - 1
		data[[x[0] for x in seqs].index(line[2])].append((line[2], start, end))
		for row in data: sorted(row, key=itemgetter(0))
	
	# Merges coordinates of overlapping aligned sequences
	print('Merging coordinates...')
	for row in data:
		i = 1
		while i < len(row):
			a = row[i-1]
			b = row[i]
			if not a[1] > b[2] and not a[2] < b[1]:
				row.remove(a)
				row.remove(b)
				row.insert(i-1, (a[0], min(a[1],b[1]), max(a[2],b[2])))
			else: i += 1
			
		# Adds filler annotations to beginning and end of sequence
		e = seqs[data.index(row)]
		row.insert(0, (e[0], 0, 0))
		row.append((e[0], e[1], e[1]))
	
	# Inverts coordinates
	print('Inverting coordinates...')
	nalign = [[] for _ in range(len(seqs))]
	for i, row in enumerate(data):
		for j in range(1, len(row)):
			a = row[j-1]
			b = row[j]
			nalign[i].append((a[0], int(a[2])+1, int(b[1])-1))
	
	# Writes to GFF file
	print('Writing to file...')
	with open(outputPath, 'w') as outFile:
		outFile.write('##gff-version 3\n')
		for seq in seqs: outFile.write('##sequence-region ' + seq[0] + ' 1 ' + seq[1] + '\n')
		temp = [el for row in nalign for el in row]
		for i, item in enumerate(temp):
			outFile.write(item[0] + '\tblastn\tnon-alignment\t' + str(item[1]) + '\t' + str(item[2]) + '\t.\t.\t.\tID=nal_' + str(i) + ';Name=nal_' + str(i) + ';Start=' + str(item[1]) + ';End=' + str(item[2]) + '\n')
	
	print('Finished!')
	
# =============================================================================

if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2])
	
	# For local testing purposes
	# from Paths import path
	# main(*path('NALs'))