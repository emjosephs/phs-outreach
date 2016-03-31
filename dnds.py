
#generate codonTable
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codonTable = dict(zip(codons, amino_acids))
codonTable[0] = 0

#syn table
synList = [0.333,0.333,0.666,0.666,1,1,1,1,0.333,0.333,0.666,0.333,0.333,0.333,0.333,0,1,1,1.333,1.333,1,1,1,1,0.333,0.333,0.333,0.333,1,1,1,1,0.333,0.333,0.333,1,1,1,1,1,0.333,0.333,0.333,0.333,0.333,0.333,0.333,0.333,1,1,1,1,1,1,1,1,0.333,0.333,0.333,0.333,1,1,1,1]
synTable = dict(zip(codons,synList))

fasta = open("all_dumptruck_aln.fasta", 'r') #open fasta file

sequenceDictionary = {}

for line in fasta: #read through fasta file and stick the sequences into a dictionary
	if line[0] == ">": #is it a name line? 
		sequenceName = line.rstrip()
		sequenceDictionary[sequenceName] = []
	else: #this is the gene sequence
		for nucleotide in line.rstrip(): #read through the nucleotides
			sequenceDictionary[sequenceName].append(nucleotide)

seq1, seq2 = sequenceDictionary.values() #take the sequences out of a dictionary
trimmed1, trimmed2 = [],[]

#trim out the dashes
for i in range(0, len(seq1)):
	if seq1[i] == "-" or seq2[i] == "-": #this is a gap
		continue
	else:
		trimmed1.append(seq1[i]) #add the aligned sequences to the new lists
		trimmed2.append(seq2[i])


#calculate how many synonymous and nonsynonymous sites there are total

def countSites(seq):
	codons = ["".join(seq[i:i+3]) for i in range(0,len(seq)/3-1)]
	ns, s = 0,0
	for codon in codons:
		synCount = synTable[codon]	
		nonsynCount = 3-synCount
		ns += nonsynCount
		s += synCount

	return(ns,s)

nsCount1, sCount1 = countSites(trimmed1)
nsCount2, sCount2 = countSites(trimmed2)

codons1 = ["".join(trimmed1[i:i+3]) for i in range(0,len(trimmed1)/3-1)]
codons2 = ["".join(trimmed2[i:i+3]) for i in range(0,len(trimmed2)/3-1)]



