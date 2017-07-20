# =============================================================================
# bmle
# GplexProject
# Given a GFF file of G-quadruplexes and a GFF file of annotations, finds the
# nearest annotation to each Gplex
# =============================================================================

def generate(gplexPath, annotPath, dataPath):
	"""Generate a data file listing nearest annotations for each gplex.

		:param gplexPath: the absolute path to the GFF-formatted gplex file
		:param annotPath: the absolute path to the GFF-formatted genome annotation file
		:param dataPath: the absolute path to where the data file should be written
		:return: nothing
	"""
	import math
	from operator import itemgetter
	from GFFLoader import load
	print('Generating data file...')
	
	# Loads files into memory
	print('Loading files...')
	gplex = load(gplexPath)
	annot = load(annotPath)
	
	# Prepares list of closest ORFs
	orfListData = []
	orfListHeaders = [['gplex-id', 'seq-id', 'start', 'end', 'gplex strand', 'closest annot', 'location', 'gplex start', 'annot end', 'distance (bp)', 'annot strand']]
	
	# Extracts sequence-regions
	seqregs = []
	for line in annot:
		if line[0].startswith('##sequence-region'):	seqregs.append(line[0].split(' ')[1])
	
	# Iterate over all gplex entries
	print('Calculating stats for each G-quadruplex...')
	gen = (line for line in gplex if len(line) == 9)
	for line in gen:
		seqid = line[0]
		start = int(line[3])
		end = int(line[4])
		strand = line[6]
		ID = line[8][0][3:]
		
		# Calculates distance to the nearest annotation
		minDist = math.inf
		pos = '?'
		orfid = '?'
		orfStrand = '?'

		for annotLine in annot:
			if annotLine[0] == seqid:
				dis = {'5\'-5\'': int(annotLine[3]) - start,
					'3\'-5\'': int(annotLine[3]) - end,
					'5\'-3\'': int(annotLine[4]) - start,
					'3\'-3\'': int(annotLine[4]) - end}
				key = min({k: abs(v) for k,v in dis.items()}.items(), key=itemgetter(1))[0]
				val = dis[key]
				
				if abs(val) < abs(minDist):
					minDist = val
					pos = key
					orfid = annotLine[8][0][3:]
					orfStrand = annotLine[6]
		temp = pos.split('-')
	
		# Calculates location relative to its nearest annotation
		if temp[0] == temp[1]: location = 'Overlap'
		elif temp[0] == '3\'' and temp[1] == '5\'' and minDist > 0: location = 'Upstream'
		elif temp[0] == '5\'' and temp[1] == '3\'' and minDist < 0: location = 'Downstream'
		else: location = 'Overlap'

		# Appends nearest annot for this g-plex to the list of annot information
		orfListData.append([ID, seqid, start, end, strand, orfid, location, temp[0], temp[1], minDist, orfStrand])
	
	# Writes orfList to output file
	print('Writing to file...')
	orfList = orfListHeaders + orfListData
	with open(dataPath, 'w') as stats:
		col_width = [max(len(str(x)) + 2 for x in line) for line in zip(*orfList)]
		for row in orfList:	stats.write(''.join(str(word).ljust(col_width[i]) for i, word in enumerate(row)).rstrip() + '\n')
	print('Finished generating data file!\n')
	

def summarize(dataPath, summaryPath):
	"""Generate summary statistics for the data file previously written.
	
	:param dataPath: the absolute path to the data file generated beforehand
	:param summaryPath: the absolute path to where the summary file should be written
	:return: nothing
	"""
	print('Generating summary file...')
	
	print('Loading data...')
	data = []
	with open(dataPath, 'r') as dataFile:
		next(dataFile)
		for line in dataFile:
			temp = line.split()
			temp[9] = int(temp[9])
			data.append(temp)
	
	# Separates contents based on seqid
	dictlol = {}
	seqregs = sorted(set(list(zip(*data))[1]))
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
		
		# ---------------------------------------------------------------------
		# Calculates average distance to nearest annotation
		# ---------------------------------------------------------------------
		print('Calculating average distance to nearest annotation...')
		distances = []
		sumAll = 0
		count = 0
		
		# Calculates averages for each sequence
		for key, value in dictlol.items():
			def mean(a): return round(sum(map(abs, a)) / len(a), 2)
			tempList = list(zip(*[item for item in value if item[6]!='Overlap']))
			if not tempList: distances.append([key+':', 'n/a'])
			else:
				distances.append([key+':', str(mean(tempList[9]))])
				sumAll += sum(map(abs, tempList[9]))
				count += len(tempList[9])
		
		# Calculates average across all sequences
		distances.insert(0, ['All sequences:', str(round(sumAll/count,2))])
		
		# Prints everything
		width = [max(len(str(x))+1 for x in col) for col in zip(*distances)]
		sumFile.write('Average distance to nearest annotation (excluding overlapping g-plexes):\n')
		for row in distances: sumFile.write('\t' + row[0].ljust(width[0]) + row[1] + ' bp\n')
		sumFile.write('\n')
		
		# ---------------------------------------------------------------------
		# Calculates relative locations of G-quadruplexes
		# ---------------------------------------------------------------------
		print('Calculating relative locations of G-quadruplexes...')
		locs = ['Upstream', 'Overlap', 'Downstream']
		toWrite = []
		sumAll = [0, 0, 0]
		
		# Calculates locations for each sequence
		for key, value in dictlol.items():
			counts = [sum(1 for x in value if x[6] == loc) for loc in locs]
			for i, loc in enumerate(locs): sumAll[i] += counts[i]
			tempSum = sum(counts)
			
			if tempSum == 0:
				toWrite.append([key+':', 'n/a (no G-quadruplexes in this sequence)', '', '', ''])
			else:
				percents = [str(round(count / tempSum * 100, 2)) for count in counts]
				toWrite.append([key+':', 'Upstream:  ', str(counts[0]), str(tempSum), percents[0]])
				toWrite.append([' '*(len(key)+1), 'Overlap:   ', str(counts[1]), str(tempSum), percents[1]])
				toWrite.append([' '*(len(key)+1), 'Downstream:', str(counts[2]), str(tempSum), percents[2]])
			
		# Calculates locations for all sequences
		sumTotal = sum(sumAll)
		percents = [str(round(count / sumTotal * 100, 2)) for count in sumAll]
		toWrite.insert(0, ['               ', 'Downstream:', str(sumAll[2]), str(sumTotal), percents[2]])
		toWrite.insert(0, ['               ', 'Overlap:   ', str(sumAll[1]), str(sumTotal), percents[1]])
		toWrite.insert(0, ['All sequences: ', 'Upstream:  ', str(sumAll[0]), str(sumTotal), percents[0]])
		
		# ---------------------------------------------------------------------
		# Writes everything
		# ---------------------------------------------------------------------
		print('Writing to file...')
		header = 'Locations of G-quadruplexes relative to their nearest ORFs:'
		writer(sumFile, header, toWrite)
	print('Finished generating summary file!\n')


def writer(fileObj, header, lol):
	"""Write the contents of a list of lists to fileObj
	
	:param fileObj: the fileObj to write to
	:param header: a header that should be written at the top of fileObj
	:param lol: the data to write to fileObj
	:return: nothing
	"""
	
	width = [max(len(str(x)) + 1 for x in line) for line in zip(*lol)]
	fileObj.write(header + '\n')
	for row in lol:
		fileObj.write('\t' + row[0].ljust(width[0]) + row[1])
		if not str(row[1]).startswith('n/a'): fileObj.write(row[2].rjust(width[2]+1) + ' /' + row[3].rjust(width[3]) + ' (' + row[4] + '%)')
		fileObj.write('\n')
	fileObj.write('\n')
	
# =============================================================================

if __name__ == '__main__':
	import sys
	generate(sys.argv[1], sys.argv[2], sys.argv[3])
	summarize(sys.argv[3], sys.argv[4])
	
	# For local testing purposes
	# from Paths import path
	# strs = path('NearestAnnot-ORFs')
	# strs = path('NearestAnnot-NALs')
	# generate(strs[0], strs[1], strs[2])
	# summarize(strs[2], strs[3])
	