# =============================================================================
# bmle
# G4Pipeline: Pipeline.py
# A wrapper module to run through all modules in the G4 annotation pipeline
# =============================================================================

def main():
	import Utils
	import BedToGFF
	import NearestAnnot
	import NonAlignments
	import GeneOverlap
	
	# =========================================================================
	# File input and formatting modules
	# =========================================================================
	
	prefix = input('Path to directory containing all files: ')
	if not prefix.endswith('/'): prefix += '/'
	
	fasta = prefix + input('Filename of genomic FASTA file: ')
	
	annot = prefix + input('Filename of genomic annotation GFF file: ')
	if annot.endswith('.gff'):
		with open(prefix + annot, 'r') as f:
			if not f.readline().startswith('##'):
				Utils.reformatGFF(annot, fasta)
				annot += '3'
	
	qb = prefix + input('Filename of QuadBase2 Tetraplex Finder BED file: ')
	BedToGFF.reformatBED(qb, fasta)
	gplex = qb[:-3] + 'gff3'
	
	sam = prefix + input('Filename of blastn SAM file: ')
	Utils.reformatSAM(sam, fasta)
	
	# =========================================================================
	# Main pipeline modules
	# =========================================================================
	
	NonAlignments.main(sam)
	nal = sam[:-3] + 'gff3'
	
	GeneOverlap.main(annot, gplex, nal)
	
	# =========================================================================
	# Summary data modules
	# =========================================================================
	
	output1 = prefix + 'analyses/gplex.txt'
	NearestAnnot.generate(gplex, annot, output1)
	NearestAnnot.summarize(output1, fasta)

	output2 = prefix + 'analyses/nal.txt'
	NearestAnnot.generate(gplex, nal, output2)
	NearestAnnot.summarize(output2, fasta)
	
# =============================================================================

if __name__=='__main__':
	main()