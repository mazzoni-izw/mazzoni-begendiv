import argparse
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser(description= "dereplicates reads in a fastq file, takes the highest quality value for each base")
parser.add_argument("fastqIn", help="fastq input file", type=str)
parser.add_argument("out",help="fastq output file(dereplicated)", type=str)
parser.add_argument("--g", "-gzip", action="store_true")


args = parser.parse_args()
if args.g:
	fqStream = gzip.open(args.fastqIn)
else:
	fqStream = open(args.fastqIn)

myWriter=open(args.out, "w")

readDICT={}
countDICT={}

print "Starting to dereplicate"
for rec in SeqIO.parse(fqStream, "fastq"):

	if not readDICT.has_key(rec.seq):
		readDICT[rec.seq] = rec
		countDICT[rec.seq] = 1

	else:
		countDICT[rec.seq] += 1
		oldrec = readDICT[rec.seq]

		bestqual =[]
		for i in range(0, len(rec.seq)):
			oldQual = oldrec.letter_annotations["phred_quality"][i]
			newQual = rec.letter_annotations["phred_quality"][i]

			if newQual > oldQual:
				bestqual.append(newQual)
				readDICT[rec.seq].letter_annotations["phred_quality"][i] = newQual
			else:
				bestqual.append(oldQual)

print "Dereplication done"


print "Starting to sort the reads by frequency"

for key, value in sorted(countDICT.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    #print "%s: %s" % (key, value)
    newID = readDICT[key].id + ";size=" + str(countDICT[key])
    readDICT[key].description = ""
    readDICT[key].id = newID
    SeqIO.write(readDICT[key], myWriter, "fastq")

print "Sorting finished"