# =============================================================================
# bmle
# GplexProject: GeneOverlap.py
# Given gene annotations, gplex annotations, and non-alignment annotations,
# generate a GFF file of genes that overlap at least one gplex and at least one
# non-alignment at the same position.
# =============================================================================

def main(gffPath, gplexPath, nalPath, minCov=0.5, maxDist=0):
	"""Generate a GFF file of genes that overlap at least one gplex and at least one non-alignment.
	
	:param gffPath: path to the GFF3-formatted gene annotation file
	:param gplexPath: path to the GFF3-formatted gplex file
	:param nalPath: path to the GFF3-formatted non-alignment file
	:param minCov: minimum overlap required of non-aligned region (default=0.5)
	:param maxDist: max number of base pairs separating a gplex and gene (default=0)
	:return: filters the three inputted files for entries that overlap each other into separate files
	"""
	import os
	from Utils import load, writeFile
	print('\nGenerating gene overlaps...')
	
	print('Loading files...')
	gplexData = load(gplexPath)[2]
	nalsData = load(nalPath)[2]
	tempGeneData = load(gffPath)
	headers = tempGeneData[0] + tempGeneData[1]
	geneData = []
	for line in tempGeneData[2]:
		if line[2]=='gene': geneData.append(line)
	
	# Finds overlaps for each ORF
	print('Calculating overlaps...')
	genes = []
	nals = []
	gplexes = []
	
	l = len(geneData)
	for i, gene in enumerate(geneData):
		print('\tCalculating ' + str(i+1) + ' of ' + str(l) + '...')
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
		if (sumCov/(gEnd-gStart) > float(minCov)) and (len(tempGplexes) > 0):
			genes.append(gene)
			nals.extend(tempNals)
			gplexes.extend(tempGplexes)
			
	# Write everything
	print('Writing to output files...')
	output = os.path.dirname(gffPath) + '/overlaps/'
	writeFile(output + 'genes.gff', headers, genes)
	writeFile(output + 'nals.gff', headers, nals)
	writeFile(output + 'gplexes.gff', headers, gplexes)
	print('Finished writing output to ' + output + '\nFinished!')

# =============================================================================

if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description='Generate a GFF file of genes that overlap at least one gplex and at least one non-alignment.')
	parser.add_argument('gffPath',
						help='path to the GFF-formatted gene annotation file')
	parser.add_argument('gplexPath',
						help='path to the GFF-formatted gplex file')
	parser.add_argument('nalPath',
						help='path to the GFF-formatted non-alignment file')
	parser.add_argument('--minCov', type=float, action='store', default=0.5,
						help='minimum overlap required of non-aligned region (default=0.5)')
	parser.add_argument('--maxDist', type=int, action='store', default=0,
						help='max number of base pairs separating a gplex and gene (default=0)')
	args = parser.parse_args()
	
	main(args.gffPath, args.gplexPath, args.nalPath, args.minCov, args.maxDist)
