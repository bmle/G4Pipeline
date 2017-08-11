# =============================================================================
# bmle
# G4Pipeline: Splitter.py
# Splits FASTA files into smaller-sized chunks to feed into applications that have limited file sizes
# =============================================================================

def main(fastaPath, maxSize):
	"""Split a FASTA file into multiple files of a specified size.
	
	:param fastaPath: path to the FASTA-formatted genome file
	:param maxSize: max size (in MB) of each split file
	:return: writes a directory of files that represent pieces of fastaPath
	"""
	import os
	from Bio import SeqIO
	print('Splitting file...')
	
	actualMaxSize = float(maxSize)*1000000		# convert maxSize to bytes
	outPath = os.path.dirname(fastaPath) + '/splitFiles/'
	os.makedirs(outPath, exist_ok=True)
	
	# file path generator
	def fileGen(outStr):
		i = 0
		while True:
			yield outStr + str(i) + '.fasta'
			i += 1
	
	lenAcc = 0
	recordAcc = []
	fGen = fileGen(outPath)
	
	for rec in SeqIO.parse(fastaPath, 'fasta'):
		seqLen = len(rec.seq)
		seqLen += round(seqLen/60)	# Includes line breaks in length calculation
		
		if seqLen > actualMaxSize:
			raise ValueError("Sequence " + rec.id + " is larger than the specified max file size (" + str(maxSize) + "); please specify a larger size!")
		elif (lenAcc + seqLen) > actualMaxSize:
			SeqIO.write(recordAcc, next(fGen), 'fasta')
			lenAcc = seqLen
			recordAcc = [rec]
		else:
			lenAcc += seqLen
			recordAcc.append(rec)
	SeqIO.write(recordAcc, next(fGen), 'fasta')
	
	print('Finished!')

# =============================================================================

if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(
		description='Splits a FASTA file into multiple files of a specified size.')
	parser.add_argument('fastaPath',
						help='path to the FASTA-formatted genome file')
	parser.add_argument('size', type=float,
						help='max size (in MB) of each split file')
	args = parser.parse_args()

	main(args.fastaPath, args.size)
	