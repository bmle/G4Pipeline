# =============================================================================
# bmle
# GplexProject: BedToGff.py
# Reformats a QuadBase2-outputted BED file into a easier-parsable GFF file
# =============================================================================

def reformat(bedPath, fastaPath, bedToGFFPath):
	"""Reformat a QuadBase2-outputted BED file to a more-parsable GFF file.
	
	:param bedPath: path to the BED-formatted QuadBase2 file
	:param fastaPath: path to the FASTA-formatted genomic sequence file
	:param bedToGFFPath: path where the reformatted file should be written to
	"""
	import os
	import errno
	from operator import itemgetter
	from natsort import natsorted
	from GFF import generateSeqRegs
	print('\nReformatting BED to GFF...')

	# Prepares new GFF file
	try:
		bedToGFF = open(bedToGFFPath, 'w')
	except IOError:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), bedToGFFPath)
	bedToGFF.write('##gff-version 3\n')
	
	# Extracts sequence-regions from GFF file and writes to file
	for line in generateSeqRegs(fastaPath): bedToGFF.write(line + '\n')
	
	# Loads BED file into memory and sorts entries by sequence id and start position
	bed = []
	with open(bedPath) as bedFile:
		for line in bedFile: bed.append(line.split('\t'))
	bed = natsorted(bed, key=itemgetter(0,1))
	
	# For determining the strand of a gplex
	def detStrand():
		if line[5].startswith('G'): return '+'
		elif line[5].startswith('C'): return '-'
		else: return '?'
	
	# Iterate over all G-plex entries in 'bed'
	for i, line in enumerate(bed):
		seqid = line[0]			# column 1: seqid
		source = 'QuadBase2'	# column 2: source
		typ = 'G_quartet'		# column 3: type
		start = int(line[1])+1	# column 4: start (converts from 0 to 1-indexing)
		end = int(line[2])		# column 5: end (end index stays the same)
		score = line[3]			# column 6: score
		strand = detStrand()	# column 7: strand
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
	
	bedToGFF.close()
	print('Finished writing output to ' + bedToGFFPath + '\nFinished reformatting!\n')
	
# =============================================================================

if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description='Reformat a QuadBase2-outputted BED file to a more-parsable GFF file.')
	parser.add_argument('bedPath',
						help='path to the BED-formatted QuadBase2 file')
	parser.add_argument('fastaPath',
						help='path to the FASTA-formatted genomic sequence file')
	parser.add_argument('bedToGFFPath',
						help='path where the reformatted file should be written to')
	args = parser.parse_args()
	
	reformat(args.bedPath, args.fastaPath, args.bedToGFFPath)