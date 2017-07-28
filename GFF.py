# =============================================================================
# bmle
# GplexProject: GFF.py
# Utilities for manipulating GFF files
# =============================================================================

def load(filePath):
	"""Load the contents of a GFF file into a list of lists.

	:param filePath: the absolute path to the GFF file
	:return: a list of lists representing the contents of the GFF file
	"""
	from operator import itemgetter
	
	headers = []
	data = []
	with open(filePath) as file:
		for line in file:
			temp = line.split('\t')
			if len(temp) == 9:
				temp[8] = temp[8].split(';')
				data.append(temp)
			else:
				headers.append(temp)
				
	data = sorted(data, key=itemgetter(0,3))
	headers.extend(data)
	return headers

def loadSeqregs(filePath):
	"""Load the sequence-regions from the GFF file into a list.
	
	:param filePath: the absolute path to the GFF file
	:return: a list of sequence-regions
	"""
	
	toReturn = []
	with open(filePath) as file:
		for line in file:
			if line.startswith('##sequence-region'): toReturn.append(line)
	return toReturn

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

	:param filePath: the path to the file to write to
	:param header: the headers of the GFF file
	:param data: the data for the file
	:return: nothing
	"""
	from itertools import groupby
	
	with open(filePath, 'w') as file:
		for line in header:
			if type(line) is list: file.write(line[0])
			else: file.write(line)
		data.sort()
		dataFiltered = list(l for l, _ in groupby(data))	# removes duplicates from data
		for line in dataFiltered: file.write(writeEntry(line))
