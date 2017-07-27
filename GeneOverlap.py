# =============================================================================
# bmle
# GplexProject: GeneOverlap.py
# Given gene annotations, gplex annotations, and non-alignment annotations,
# generate a GFF file of genes that overlap at least one gplex and at least one
# non-alignment at the same position.
# =============================================================================

def main(gffPath, gplexPath, nalPath, geneOPath, nalOPath, gplexOPath):
	"""Generate a GFF file of genes that overlap at least one gplex and at least one non-alignment.
	
	:param gffPath: the absolute path to the input gene annotation GFF file
	:param gplexPath: the absolute path to the input gplex GFF file
	:param nalPath: the absolute path to the input non-alignment GFF file
	:param geneOPath: the absolute path to the output gene GFF file
	:param nalOPath: the absolute path to the output non-alignment GFF file
	:param gplexOPath: the absolute path to the output gplex GFF file
	:return: nothing
	"""
	from GFF import load, writeFile
	
	print('Loading files...')
	nalsData = load(nalPath)
	gplexData = load(gplexPath)
	gffData = load(gffPath)
	
	# Separates annotation headers from data
	gffHeaders = []
	gffs = []
	for line in gffData:
		if line[0].startswith(('##', '#!')): gffHeaders.append(line)
		else: gffs.append(line)
	
	# Finds overlaps for each ORF
	print('Calculating overlaps...')
	genes = []
	nals = []
	gplexes = []
	minCov = 0.40	# minimum overlap required of non-aligned region
	maxDist = 25	# max number of base pairs away a gplex can be from a gene
	
	for gene in gffs:
		tempNals = []
		tempGplexes = []
		sumCov = 0
		gStart = int(gene[3])
		gEnd = int(gene[4])
		
		# Iterates over each non-alignment
		for nline in nalsData:
			if (nline not in tempNals) and (nline[0] == gene[0]):
				start = max(gStart, int(nline[3]))
				end = min(gEnd, int(nline[4]))
				if end > start:
					tempNals.append(nline)
					sumCov += (end-start)
		
		# Iterates over each gplex
		for gplex in gplexData:
			if gplex[0] == gene[0]:
				start = max(gStart, int(gplex[3]))
				end = min(gEnd, int(gplex[4]))
				if start-end <= maxDist:	tempGplexes.append(gplex)

		# If coverage is at least 'minCov' and there exists at least one gplex, add to data
		if (sumCov/(gEnd-gStart) > minCov) and (len(tempGplexes) > 0):
			genes.append(gene)
			nals.extend(tempNals)
			gplexes.extend(tempGplexes)
			
	# Write everything
	print('Writing to output files...')
	writeFile(geneOPath, gffHeaders, genes)
	writeFile(nalOPath, gffHeaders, nals)
	writeFile(gplexOPath, gffHeaders, gplexes)
	print('Finished!')

# =============================================================================

if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
	
	# For local testing purposes
	# from Paths import path
	# main(*path('GeneOverlap'))