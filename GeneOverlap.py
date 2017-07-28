# =============================================================================
# bmle
# GplexProject: GeneOverlap.py
# Given gene annotations, gplex annotations, and non-alignment annotations,
# generate a GFF file of genes that overlap at least one gplex and at least one
# non-alignment at the same position.
# =============================================================================

def main(gffPath, gplexPath, nalPath, minCov=0.4, maxDist=25):
	"""Generate a GFF file of genes that overlap at least one gplex and at least one non-alignment.
	
	:param gffPath: the absolute path to the input gene annotation GFF file
	:param gplexPath: the absolute path to the input gplex GFF file
	:param nalPath: the absolute path to the input non-alignment GFF file
	:param minCov: minimum overlap required of non-aligned region
	:param maxDist: max number of base pairs away a gplex can be from a gene
	:return: nothing
	"""
	import os
	from GFF import load, writeFile
	
	print('Loading files...')
	gplexData = load(gplexPath)[2]
	nalsData = load(nalPath)[2]
	tempGeneData = load(gffPath)
	headers = tempGeneData[0] + tempGeneData[1]
	geneData = tempGeneData[2]
	
	# Finds overlaps for each ORF
	print('Calculating overlaps...')
	genes = []
	nals = []
	gplexes = []
	
	for gene in geneData:
		tempNals = []
		tempGplexes = []
		sumCov = 0
		gStart = int(gene[3])
		gEnd = int(gene[4])
		
		# Iterates over each non-alignment
		for nline in nalsData:
			if nline[0] == gene[0]:
				start = max(gStart, int(nline[3]))
				end = min(gEnd, int(nline[4]))
				if end > start:
					tempNals.append(nline)
					sumCov += (end-start)
		
		# Iterates over each gplex; gplex is included if its distance to the orf doesn't exceed 'maxDist'
		for gplex in gplexData:
			if gplex[0] == gene[0]:
				start = max(gStart, int(gplex[3]))
				end = min(gEnd, int(gplex[4]))
				if (start-end) <= int(maxDist):	tempGplexes.append(gplex)

		# If coverage is at least 'minCov' and there exists at least one gplex, add to data
		if (sumCov/(gEnd-gStart) > int(minCov)) and (len(tempGplexes) > 0):
			genes.append(gene)
			nals.extend(tempNals)
			gplexes.extend(tempGplexes)
			
	# Write everything
	print('Writing to output files...')
	output = os.path.dirname(gffPath) + '/overlaps/'
	if not os.path.exists(output): os.makedirs(output)
	
	writeFile(output + 'genes.gff', headers, genes)
	writeFile(output + 'nals.gff', headers, nals)
	writeFile(output + 'gplexes.gff', headers, gplexes)
	print('Finished writing output to ' + output + '\nFinished!')

# =============================================================================

if __name__ == '__main__':
	# import sys
	# main(*sys.argv[1:])
	
	# For local testing purposes
	from Paths import path
	main(*path('GeneOverlap'))