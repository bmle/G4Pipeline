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
	
	toReturn = []
	with open(filePath) as file:
		for line in file:
			temp = line.split('\t')
			if len(temp) == 9: temp[8] = temp[8].split(';')
			toReturn.append(temp)
	return toReturn

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