import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="filters a fastq file to start and end with defined sequences (designed for merged ddRAD data)")
parser.add_argument("fqIn", help="fastq input file", type=str)
parser.add_argument("outCorrect", help="output file for sequences with correct restriction sites", type=str)
parser.add_argument("info", help="output file for numebr of reads with correct or incorrect restriction sites", type=str)
parser.add_argument("res1", help="overhang of 1. enzmye --> starting basepairs for all reads", type=str)
parser.add_argument("res2", help="overhang of 2. enzmye --> ending basepairs for all reads", type=str)

args = parser.parse_args()

myFQ = open(args.fqIn)

correctOut = open(args.outCorrect, "w")
errorOut = open(args.info, "w")

errorDICT_1 = {}
errorDICT_2 = {}

seqCount = 0
seqCorrect = 0
errorsRes1 = 0
errorsRes2 = 0

firstLineStorage = ""

for record in SeqIO.parse(myFQ, "fastq"):
	seqCount += 1

	curSeq = record.seq
	curRes1 = curSeq[0:len(args.res1)]
	curRes2 = curSeq[-len(args.res2):]



	if args.res1 != curRes1:
		if errorDICT_1.get(str(curRes1), "empty") == "empty":
			errorDICT_1[str(curRes1)] = 1
		else:
			errorDICT_1[str(curRes1)] += 1
		#print args.res1, curRes1
		#print record.id
		errorsRes1 += 1

	if args.res2 != curRes2:
		if errorDICT_2.get(str(curRes2), "empty") == "empty":
			errorDICT_2[str(curRes2)] = 1
		else:
			errorDICT_2[str(curRes2)] += 1

		#print args.res2, curRes2
		#print record.id
		errorsRes2 += 1

	elif args.res1 == curRes1 and args.res2 == curRes2:
		SeqIO.write(record, correctOut, "fastq")
		seqCorrect += 1
		#print "GOOD site"

errorOut.write("Total sequences: %i\n" %(seqCount))
errorOut.write("Sequences with correct sites written: %i\n" %(seqCorrect))
errorOut.write("Error sequences: %i\n\n" %(seqCount-seqCorrect))

errorOut.write("Total errors in cutting sites(combined): %i\n"%(errorsRes1+errorsRes2))
errorOut.write("Total errors in forward cutting site: %i\n"%(errorsRes1))
errorOut.write("Total errors in reverse cutting site: %i\n"%(errorsRes2))

errorOut.write("Correct forward cutting sites: %i\n"%(seqCorrect-errorsRes1))
errorOut.write("Correct reverse cutting sites: %i\n"%(seqCorrect-errorsRes2))
errorOut.write("\n\n")

errorOut.write("forwardErrors\tCounts\n")


sortr1 = [(i,errorDICT_1[i]) for i in errorDICT_1]
sortr1.sort(key=lambda x: x[1], reverse=True)
print sortr1

for val in sortr1:
	errorOut.write(val[0] +"\t" + str(val[1]) + "\n")
	#print key, errorDICT_1[key]

errorOut.write("reverseErrors\tCounts\n")

sortr2 = [(i,errorDICT_2[i]) for i in errorDICT_2]
sortr2.sort(key=lambda x: x[1], reverse=True)

for val in sortr2:
	errorOut.write(val[0] +"\t" + str(val[1]) + "\n")
	#print key, errorDICT_2[key]