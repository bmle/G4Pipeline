# =============================================================================
# bmle
# GplexProject: Splitter.py
# Splits FASTA files into smaller-sized chunks to feed into applications that have limited file sizes=============================================================================

def main(fastaPath, maxSize):
	"""Split a FASTA file into multiple files of a specified size.
	
	:param fastaPath: path to the FASTA-formatted genome file
	:param maxSize: max size (in MB) of each split file
	:return: nothing
	"""
	import os
	from Bio import SeqIO
	print('Splitting file...')
	
	maxSize *= 1000000		# convert maxSize to bytes
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
		if seqLen > maxSize:
			raise ValueError("Sequence " + rec.id + " is larger than the specified max file size (" + str(maxSize/1000000) + "); please specify a larger size!")
		elif (lenAcc + seqLen) > maxSize:
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
	