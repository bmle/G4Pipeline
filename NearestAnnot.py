# =============================================================================
# bmle
# GplexProject: NearestAnnot.py
# Finds the nearest annotation to each Gplex
# =============================================================================

def generate(gplexPath, annotPath, dataPath):
	"""Generate a data file listing nearest annotations for each gplex.

		:param gplexPath: path to the GFF-formatted gplex file
		:param annotPath: path to the GFF-formatted gene annotation file
		:param dataPath: path to where the output data file should be written
		:return: nothing
	"""
	import os
	import math
	from operator import itemgetter
	from Utils import load
	print('\nGenerating data file...')
	
	# Loads data
	print('Loading files...')
	gplex = load(gplexPath)[2]
	annot = load(annotPath)[2]
	
	# Prepares list of closest ORFs
	orfListData = []
	orfListHeaders = [['gplex-id', 'seq-id', 'start', 'end', 'gplex strand', 'closest annot', 'location', 'gplex start', 'annot end', 'distance (bp)', 'annot strand']]

	# Iterate over all gplex entries
	print('Calculating stats for each G-quadruplex...')
	l = len(gplex)
	for i, line in enumerate(gplex):
		print('\tCalculating ' + str(i+1) + ' of ' + str(l) + '...')
		seqid = line[0]
		start = int(line[3])
		end = int(line[4])
		strand = line[6]
		ID = line[8][0][3:]
		
		# Calculates distance to the nearest annotation
		minDist = math.inf
		pos = 'n/a'
		annotid = 'n/a'
		annotStrand = 'n/a'
		
		# Iterate over every annotation
		gen = (line for line in annot if (line[0]==seqid and (line[2]=='gene' or line[
			2]=='non-alignment')))
		for annotLine in gen:
			dis = {'5\'-5\'': int(annotLine[3]) - start,
				'3\'-5\'': int(annotLine[3]) - end,
				'5\'-3\'': int(annotLine[4]) - start,
				'3\'-3\'': int(annotLine[4]) - end}
			key = min({k: abs(v) for k,v in dis.items()}.items(), key=itemgetter(1))[0]
			val = dis[key]
			
			# If the distance to this annotation is smaller than the currently-recorded annotation, replace the old annot with this one
			if abs(val) < abs(minDist):
				minDist = val
				pos = key
				annotid = annotLine[8][0][3:]
				annotStrand = annotLine[6]
			# This works because $gen (and $annot) are sorted by start + end
			# positions, so any successive entries will only get further away
			elif abs(val) > abs(minDist):
				break
		
		# Calculates location relative to its nearest annotation
		temp = pos.split('-')
		if temp[0] == 'n/a':
			location = 'n/a (no annotations on this sequence)'
			temp.append('n/a')
		elif temp[0] == temp[1]: location = 'Overlap'
		elif temp[0] == '3\'' and temp[1] == '5\'' and minDist > 0: location = 'Upstream'
		elif temp[0] == '5\'' and temp[1] == '3\'' and minDist < 0: location = 'Downstream'
		else: location = 'Overlap'

		# Appends nearest annot for this g-plex to the list of annot information
		orfListData.append([ID, seqid, start, end, strand, annotid, location, temp[0], temp[1], minDist, annotStrand])
	
	# Writes orfList to output file
	print('Writing to file...')
	orfList = orfListHeaders + orfListData
	os.makedirs(os.path.dirname(dataPath), exist_ok=True)
	
	with open(dataPath, 'w') as stats:
		col_width = [max(len(str(x)) + 2 for x in line) for line in zip(*orfList)]
		for row in orfList:	stats.write(''.join(str(word).ljust(col_width[i]) for i, word in enumerate(row)).rstrip() + '\n')
	print('Finished writing to ' + dataPath)
	print('Finished generating data file!\n')
	

def summarize(dataPath, fastaPath, summaryPath):
	"""Generate summary statistics for the data file previously written.
	
	:param dataPath: path to where the output data file is written
	:param fastaPath: path to the FASTA-formatted genomic sequence file
	:param summaryPath: path to where the output summary file should be written
	:return: nothing
	"""
	from Utils import generateSeqRegs
	print('\nGenerating summary file...')
	
	print('Loading data...')
	data = []
	with open(dataPath, 'r') as dataFile:
		next(dataFile)
		for line in dataFile:
			temp = line.split()
			temp[9] = int(temp[9]) if temp[9]!='on' else temp[9]
			data.append(temp)

	# Separates contents based on sequence regions
	dictlol = {}
	seqregs = [line.split(' ')[1] for line in generateSeqRegs(fastaPath)]
	for seq in seqregs:	dictlol[seq] = [row for row in data if row[1] == seq]
	
	# Writes to file
	with open(summaryPath, 'w') as sumFile:
		
		# ---------------------------------------------------------------------
		# Calculates total number of gplexes
		# ---------------------------------------------------------------------
		print('Calculating total number of G-quadruplexes...')
		countORFs = {key: len(value) for key, value in dictlol.items()}
		totalORFs = str(sum([value for key, value in countORFs.items()]))
		sumFile.write('Total number of G-quadruplexes:\n\tAll sequences: ' + totalORFs + '\n' + '\n'.join(['\t' + key + ': ' + str(value) for key, value in countORFs.items()]) + '\n\n')
		
		# ---------------------------------------------------------------------
		# Calculates strandedness of gplexes
		# ---------------------------------------------------------------------
		print('Calculating strandedness of G-quadruplexes...')
		toWrite = []
		signs = ['+', '-']
		totals = [0, 0]
		
		# Iterates over all sequences
		for key, value in dictlol.items():
			strands = []
			for s in signs: strands.append(sum(1 for x in value if x[4]==s))
			sums = sum(strands)
			totals[0] += strands[0]
			totals[1] += strands[1]
			
			# Check for lack of g-plexes in current sequence
			if sums == 0:
				toWrite.append([key+':', 'n/a (no G-quadruplexes in this sequence)', '', '', ''])
			else:
				pc = [str(round(st / sums * 100, 2)) for st in strands]
				toWrite.append([key+':', 'Sense:    ', str(strands[0]), str(sums), pc[0]])
				toWrite.append([' '*len(key), 'Antisense:', str(strands[1]), str(sums), pc[1]])
		
		# Calculates strandedness for all sequences
		totalSum = sum(totals)
		percents = [str(round(count / totalSum * 100, 2)) for count in totals]
		toWrite.insert(0, ['              ', 'Antisense:', str(totals[1]), str(totalSum), percents[1]])
		toWrite.insert(0, ['All sequences:', 'Sense:    ', str(totals[0]), str(totalSum), percents[0]])
		
		# Writes everything
		header = 'Strandedness of G-quadruplexes:'
		writer(sumFile, header, toWrite)
		sumFile.write('\n')
		
		# ---------------------------------------------------------------------
		# Calculates average distance to nearest annotation
		# ---------------------------------------------------------------------
		print('Calculating average distance to nearest annotation...')
		distances = []
		sumAll = 0
		count = 0
		
		# Calculates averages for each sequence
		def mean(a): return round(sum(map(abs, a)) / len(a), 2)
		for key, value in dictlol.items():
			if not value:
				distances.append([key + ':', 'n/a (no g-quadruplexes on this sequence)'])
			else:
				naList = list(zip(*[item for item in value if item[6] == 'n/a']))
				if naList:
					distances.append([key+':', 'n/a (no annotations for this sequence)'])
				else:
					tempList = list(zip(*[item for item in value if item[6]!='Overlap']))
					if not tempList:
						distances.append([key+':', 'n/a (no non-overlapping g-quadruplexes)'])
					else:
						distances.append([key+':', str(mean(tempList[9])) + ' bp'])
						sumAll += sum(map(abs, tempList[9]))
						count += len(tempList[9])
		
		# Calculates average across all sequences
		distances.insert(0, ['All sequences:', str(round(sumAll/count,2)) + ' bp'])
		
		# Writes everything
		width = [max(len(str(x))+1 for x in col) for col in zip(*distances)]
		sumFile.write('Average distance to nearest annotation (excluding overlapping g-plexes):\n')
		for row in distances: sumFile.write('\t' + row[0].ljust(width[0]) + row[1] + '\n')
		sumFile.write('\n')
		
		# ---------------------------------------------------------------------
		# Calculates relative locations of G-quadruplexes
		# ---------------------------------------------------------------------
		print('Calculating relative locations of G-quadruplexes...')
		locs = ['Upstream', 'Overlap', 'Downstream', 'n/a']
		toWrite = []
		sumAll = [0 for _ in locs]
		
		# Calculates locations for each sequence
		for key, value in dictlol.items():
			counts = [sum(1 for x in value if x[6] == loc) for loc in locs]
			tempSum = sum(counts)
			for i in range(len(locs)): sumAll[i] += counts[i]
			
			if tempSum == 0:
				toWrite.append([key+':', 'n/a (no G-quadruplexes in this sequence)', '', '', ''])
			elif counts[3] > 0:
				toWrite.append([key+':', 'n/a (no annotations in this sequence)', '', '', ''])
			else:
				percents = [str(round(count / tempSum * 100, 2)) for count in counts]
				toWrite.append([key+':', 'Upstream:  ', str(counts[0]), str(tempSum), percents[0]])
				toWrite.append([' '*(len(key)+1), 'Overlap:   ', str(counts[1]), str(tempSum), percents[1]])
				toWrite.append([' '*(len(key)+1), 'Downstream:', str(counts[2]), str(tempSum), percents[2]])
			
		# Calculates locations for all sequences
		sumTotal = sum(sumAll)
		percents = [str(round(count / sumTotal * 100, 2)) for count in sumAll]
		toWrite.insert(0, ['               ', 'n/a:       ', str(sumAll[3]), str(sumTotal), percents[3]])
		toWrite.insert(0, ['               ', 'Downstream:', str(sumAll[2]), str(sumTotal), percents[2]])
		toWrite.insert(0, ['               ', 'Overlap:   ', str(sumAll[1]), str(sumTotal), percents[1]])
		toWrite.insert(0, ['All sequences: ', 'Upstream:  ', str(sumAll[0]), str(sumTotal), percents[0]])
		
		# ---------------------------------------------------------------------
		# Writes everything
		# ---------------------------------------------------------------------
		print('Writing to file...')
		header = 'Locations of G-quadruplexes relative to their nearest annotations:'
		writer(sumFile, header, toWrite)
	
	print('Finished writing to ' + summaryPath)
	print('Finished generating summary file!\n')


def writer(fileObj, header, lol):
	"""Write the contents of a list of lists to fileObj
	
	:param fileObj: the fileObj to write to
	:param header: a header that should be written at the top of fileObj
	:param lol: the data to write to fileObj
		Takes the form of [sequence, category, amount, total, percentage]
	:return: nothing
	"""
	
	width = [max(len(str(x)) + 1 for x in line) for line in zip(*lol)]
	fileObj.write(header + '\n')
	for row in lol:
		fileObj.write('\t' + row[0].ljust(width[0]) + row[1])
		if not str(row[1]).startswith('n/a ('): fileObj.write(row[2].rjust(width[2]+1) + ' /' + row[3].rjust(width[3]) + ' (' + row[4] + '%)')
		fileObj.write('\n')
	
# =============================================================================

if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(
		description='Finds the nearest annotation to each Gplex.')
	parser.add_argument('gplexPath',
						help='path to the GFF-formatted gplex file')
	parser.add_argument('annotPath',
						help='path to the GFF-formatted gene annotation file')
	parser.add_argument('dataPath',
						help='path to where the output data file should be written')
	parser.add_argument('fastaPath',
						help='path to the FASTA-formatted genomic sequence file')
	parser.add_argument('summaryPath',
						help='path to where the output summary file should be written')
	args = parser.parse_args()
	
	generate(args.gplexPath, args.annotPath, args.dataPath)
	summarize(args.dataPath, args.fastaPath, args.summaryPath)