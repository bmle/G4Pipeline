# =============================================================================
# bmle
# GplexProject: BedToGff.py
# Reformats a QuadBase2-outputted BED file into a easier-parsable GFF file
# =============================================================================

def reformat(bedPath, gffPath, bedToGFFPath):
	"""Reformat a QuadBase2-outputted BED file to a more-parsable GFF file.
	
	:param bedPath: the absolute path to the QuadBase2 file
	:param gffPath: the absolute path to the genomic annotation gff file
	:param bedToGFFPath: the absolute path to where the new gff file should be written
	:return: nothing
	"""
	from GFF import loadSeqregs
	print('Reformatting...')

	# Prepares new GFF file
	bedToGFF = open(bedToGFFPath, 'w')
	bedToGFF.write('##gff-version 3\n')
	
	# Extracts sequence-regions from GFF file and writes to file
	for line in loadSeqregs(gffPath): bedToGFF.write(line)
	
	# Loads BED file into memory and sorts entries by sequence id and start position
	bed = []
	with open(bedPath) as bedFile:
		for line in bedFile: bed.append(line.split('\t'))
	bed.sort(key=lambda x: (x[0], int(x[1])))
	
	# Iterate over all G-plex entries in bed array
	for i, line in enumerate(bed):
		seqid = line[0]			# column 1: seqid
		source = 'QuadBase2'	# column 2: source
		typ = 'G_quartet'		# column 3: type
		start = int(line[1])+1	# column 4: start (converts from 0 to 1-indexing)
		end = int(line[2])		# column 5: end (end index stays the same)
		score = line[3]			# column 6: score
		def detStrand():		# column 7: strand
			if line[5].startswith('G'): return '+'
			elif line[5].startswith('C'): return '-'
			else: return '?'
		strand = detStrand()
		phase = '.'				# column 8: phase
		ID = i					# column 9: attributes
		nm = i
		motif = line[4]
		sequence = line[5].strip()
		
		# Write everything
		bedToGFF.write(seqid + '\t' + source + '\t' + typ + '\t' + str(start) + '\t' + str(
			end) + '\t' + score + '\t' + strand + '\t' + phase + '\tID=gplex_' + str(
			ID) + ';Name=gplex_' + str(
			nm) + ';motif=' + motif + ';sequence=' + sequence + ';start=' + str(
			start) + ';end=' + str(end) + '\n')
	
	print('Finished!')
	
# =============================================================================

if __name__ == '__main__':
	import sys
	reformat(sys.argv[1], sys.argv[2], sys.argv[3])
	
	# For local testing purposes
	# from Paths import path
	# reformat(*path('BedToGFF'))