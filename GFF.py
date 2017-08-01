# =============================================================================
# bmle
# GplexProject: GFF.py
# Utilities for manipulating GFF files
# =============================================================================

def load(filePath):
	"""Load the contents of a GFF file.

	:param filePath: the absolute path to the GFF file
	:return: a list of lists representing the contents of the GFF file
	"""
	from operator import itemgetter
	from natsort import natsorted
	
	headerList = []
	seqregList = []
	dataList = []
	with open(filePath) as file:
		for line in file:
			temp = line.split('\t')
			if len(temp) == 9:
				temp[8] = temp[8].split(';')
				dataList.append(temp)
			elif temp[0].startswith('##sequence-region'):
				seqregList.append(line)
			elif temp[0].startswith('#'):
				headerList.append(line)
			else:
				raise AssertionError('Unknown line in GFF file: ' + line)
	
	# Sorts by: seqid -> start position -> end position
	dataList = natsorted(dataList, key=itemgetter(0,3,4))
	seqregList = natsorted(seqregList)
	return [headerList, seqregList, dataList]

def writeEntry(line):
	"""Convert a GFF-formatted entry into a string.
	GFF-formatted entry: [seqid, source, ..., [list of attributes]]
	
	:param line: A GFF-formatted entry
	:return: A string representation of the GFF-formatted entry
	"""
	toReturn = ''
	for item in line:
		if type(item) is list:
			for item2 in item: toReturn += str(item2) + ';'
		else:
			toReturn += str(item) + '\t'
	return toReturn[:-1]

def writeFile(filePath, header, data):
	"""Write data to a GFF file.

	:param filePath: the absolute path to the file to write to
	:param header: the headers of the GFF file
	:param data: the data for the file
	:return: nothing
	"""
	from natsort import natsorted
	from itertools import groupby
	
	with open(filePath, 'w') as file:
		for line in header: file.write(line)
		data = natsorted(data)
		dataFiltered = list(l for l,_ in groupby(data))	# removes duplicates from data
		for line in dataFiltered: file.write(writeEntry(line))

def generateSeqRegs(fastaPath):
	"""Generate sequence headers from a FASTA file.
	
	:param fastaPath: the absolute path to the FASTA file
	:return: a list of sequence-regions (formatted as strings)
	"""
	import re
	from natsort import natsorted
	
	tempList = []
	pattern = '[>\|,\s]+'
	f = open(fastaPath, 'r')
	
	# Prompts the user to specify the location of the sequence region name
	s = f.readline().strip()
	fline = re.split(pattern, s)
	print(fline)
	index = int(input('Index of position that contains sequence label: '))
	tempList.append([fline[index], 0])
	
	# Iterates over the rest of the strings
	for line in f:
		if not line.startswith('>'):
			tempList[-1][1] += len(line.strip())
		else:
			tempList.append([re.split(pattern, line.strip())[index], 0])
	f.close()
	
	# Builds the strings
	toReturn = []
	for pair in tempList:
		toReturn.append('##sequence-region ' + pair[0] + ' 1 ' + str(pair[1]))
	return natsorted(toReturn)
	