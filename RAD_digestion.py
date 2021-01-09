#!/usr/bin/python

##################################################################################
#### property of BeGenDiv(https://begendiv.de/)								 #####
#### @authors: Harald Detering, Marie Jeschek, Maximilian Driller            #####
#### contact: drillermax@gmail.com                                           #####
##################################################################################


import argparse, os, sys
from Bio.Restriction import RestrictionBatch, Analysis
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from sets import Set
from collections import deque

from dig_functions import *






# parse arguments
parser = argparse.ArgumentParser(description="In-silico RAD digestion.", formatter_class=argparse.RawTextHelpFormatter)
#parser.add_argument("-g", required=True, metavar="genomeFile", type=argparse.FileType('r'),
#	help="Sequences to digest in fasta or fastq format (f.e. a genome file).")
parser.add_argument("-g", required=True, metavar="genomeFile", type=str,
	help="Sequences to digest in fasta or fastq format (f.e. a genome file).")
parser.add_argument("-e", required=True, metavar="enzymeFile", type=argparse.FileType('r'),
	help="List with enzymes/enzyme combinations including restriction pattern. Tab-separated.\nFormat:\nenzyme1\tpattern1\tenzyme2\tpattern2 -- the first enzyme is the forward enzyme, the second one in the list is reverse (format is as shown in the file downloaded by the get_rad_enzymes.py script)")
parser.add_argument("--dd", action='store_true',
	help="If argument is set, only fragments cut by two different enzymes will be returned.")
parser.add_argument("-o", default=sys.stdout, metavar="outputFile", type=mkOutDir,
	help="Output file prefix. [Default: stdout]")
parser.add_argument("--q", "-fastq", action='store_true',
	help="Genome file is in fastq format. Output will also be in fastq-format (qual. Phred33)")
parser.add_argument("--min", type=int,
	help="Minimum size for fragments.")
parser.add_argument("--max", type=int,
	help="Maximum size for fragments.")
parser.add_argument("--gz", "-gzip", action='store_true',
	help="Enable if input file is gzipped.")
parser.add_argument("--v", "-verbose", action='store_true',
	help="Enables progress information in console.")
parser.add_argument("--rad", "-radseq", action='store_true',
	help="Allow fragments to have Start(=first/forward enzyme in list) and End(=second/reverse enzyme in list) of a scaffold as a restriction side. Enable if digesting already digested sequences.")
args = parser.parse_args()


enzymeCombinations = readEnzymeListWithPattern(args.e)



# iterate over contigs (from genome or reads)
outFiles={}
#dictionary for number of fragments per enzyme combintaion
fragDict={}
#dictionary for number of fragments per length 
sizeDICT={}
#variable to keep track of the longest fragment --> needed for output later
maxFrag=0

if args.gz:
	import gzip
	contigParser=gzip.open(args.g)
else:
	contigParser=open(args.g)

if args.q:
	contigs=SeqIO.parse(contigParser, "fastq")
else:
	contigs=SeqIO.parse(contigParser, "fasta")

for contig in contigs:
	if args.v:
		print("Digesting sequence %s ..." % contig.id)
	
	# for each enzyme combination:
	for enzymeComb in enzymeCombinations:
		#initialze sizeDICT
		
		if args.v:
			print("\t... with enzyme(s) %s" % ", ".join([x.name for x in enzymeComb]))
		
		# find cutting sites of enzymes ...
		cuttingSites = {}
		for enzyme in enzymeComb:
			#print enzyme, enzyme.pPattern, enzyme.offset, enzyme.addBack
			fwdCuttingPositions = [m.start() for m in re.finditer(enzyme.pPattern, str(contig.seq[:-1]).upper())]
			#print fwdCuttingPositions		
			# position of the cut is the starting position of the match + the offset (determined by the cutting site)
			fwdCuttingSites = [myCuttingSite(x+enzyme.offset, enzyme, True) for x in fwdCuttingPositions]
			#print fwdCuttingSites

			#print fwdCuttingPositions
			#if not args.dd:
			# AAACCGCTT is reverse: AAGCGGTTT; length 9
			# pattern GCGG found at rev.pos 2
			# is equivalent to fwd.pos 9-2-4=3
			
			revCuttingPositions = [m.start() for m in re.finditer(enzyme.pPattern, str(contig.seq[:-1].reverse_complement()).upper())]
			revCuttingSites = [myCuttingSite(len(contig.seq)-x-len(enzyme.site)+enzyme.offset, enzyme, False) for x in revCuttingPositions]
			#revCuttingSites = [myCuttingSite(len(contig.seq)-x-len(enzyme.site), enzyme, False) for x in revCuttingPositions]#added offset of -1 caused 1 base missing
			#print revCuttingPositions


			# to return also fragments from the beginning/end of the contig to/from the first/last restriction site:
			if args.rad:
				# given Start and order of 0 = forward enzyme and End an order of 2 = reverse Enzyme --> so no fragments Start - forward enzyme and vice versa
				firstSite = myCuttingSite(0, myRestrictionEnzyme("Start","",0), True)
				#print(firstSite)
				if firstSite.position not in cuttingSites: cuttingSites[firstSite.position]=[firstSite]

				lastSite = myCuttingSite(len(contig), myRestrictionEnzyme("End","",2), True)
				if lastSite.position not in cuttingSites: cuttingSites[lastSite.position]=[lastSite]

			# rev cutting sited not used anymore
			for cuttingSite in fwdCuttingSites+revCuttingSites:
				#print(cuttingSite)
				cuttingSites.setdefault(cuttingSite.position, []).append(cuttingSite)
				#print cuttingSite, cuttingSite.position
			
			
			#for cuttingSite in fwdCuttingSites:
			#	cuttingSites.setdefault(cuttingSite.position, []).append(cuttingSite)
				#print cuttingSite, cuttingSite.position
				

		if args.v:
			# ... and return fragments
			print("\tWriting fragments.")


		outFileKey = "_".join([x.name for x in enzymeComb])
		if outFileKey not in outFiles:
			try:
				if args.q:
					if args.min and args.max:
						outFiles[outFileKey] = open(args.o+"%s_%i-%i.fastq" % ("_".join([x.name for x in enzymeComb]), args.min, args.max),"w")
					elif args.min:
						outFiles[outFileKey] = open(args.o+"%s_%i-.fastq" % ("_".join([x.name for x in enzymeComb]), args.min),"w")
					elif args.max:
						outFiles[outFileKey] = open(args.o+"%s_-%i.fastq" % ("_".join([x.name for x in enzymeComb]), args.max),"w")
					else:
						outFiles[outFileKey] = open(args.o+"%s.fastq" % "_".join([x.name for x in enzymeComb]),"w")

				else:
					if args.min and args.max:
						outFiles[outFileKey] = open(args.o+"%s_%i-%i.fasta" % ("_".join([x.name for x in enzymeComb]), args.min, args.max),"w")
					elif args.min:
						outFiles[outFileKey] = open(args.o+"%s_%i-.fasta" % ("_".join([x.name for x in enzymeComb]), args.min),"w")
					elif args.max:
						outFiles[outFileKey] = open(args.o+"%s_-%i.fasta" % ("_".join([x.name for x in enzymeComb]), args.max),"w")
					else:
						outFiles[outFileKey] = open(args.o+"%s.fasta" % "_".join([x.name for x in enzymeComb]),"w")

			except AttributeError:
				outFiles[outFileKey] = sys.stdout

		lastSites = None
		i=1
		for site in sorted(cuttingSites):
			unifiedSites=unifySites(cuttingSites[site])
			#print unifiedSites
			if lastSites:
				for start in lastSites:
					for end in unifiedSites:
						#print start, end


						#if not args.dd or start.enzyme!=end.enzyme:
						#comparing the order of the fragments so enzymes so Start - forward enz gets ignored and End - reverse also for option --rad
						if not args.dd or start.enzyme.order!=end.enzyme.order:
						#if not args.dd:
							
							#fragment = "%s%s" % (contig.seq[start.position:end.position], oHang)
							if end.enzyme.front==0:
								oHang=end.completeOverhang()
								fragment = contig.seq[start.position:end.position] + oHang
							else:
								fragment = contig.seq[start.position:end.position]

							if start.enzyme.front==1:
								oHang=start.completeOverhang()
								fragment= oHang+fragment

							#print(start.enzyme.name, end.enzyme.name)
							#print(start.completeOverhang(), end.completeOverhang())
							#print(start.enzyme.pattern, end.enzyme.pattern)
							#print(start.enzyme.offset, end.enzyme.offset)

							#sequence needs to be the reverse complement if the reverse enzyme cuts first
							if start.enzyme.order > end.enzyme.order:
									#print fragment
									mySeq = Seq.Seq(str(fragment), generic_dna)
							#		
									fragment = mySeq.reverse_complement()[0:]
									#fragment=fragment.reverse_complement()
									#print fragment

									# add string to second enzyme so it is in the name if sequence is reversed
									#end.enzyme.name=end.enzyme.name+"(reversed)"


							if args.q:
								#get quality string of sequence and trasnform the in vaues to phred33
								quals = contig.letter_annotations["phred_quality"][start.position:end.position+len(oHang)]
								#quals = contig.letter_annotations["phred_quality"][start.position:end.position]

								qual=""
								for bqual in quals:
									qual += chr(bqual+33)
							
							
							if not args.min or len(fragment)>=args.min:
								#print fragment
								if not args.max or len(fragment)<=args.max:
									
									if args.q:
										#add/count fragments in dictionary
										enzyme_combi = start.enzyme.name  + "_" + end.enzyme.name + ("(" + enzymeComb[0].name +"-"+ enzymeComb[1].name +")")
										enzNames=enzymeComb[0].name +"-"+ enzymeComb[1].name
										#if start.enzyme.name == "Start" or end.enzyme.name =="End":
										#enzyme_combi = enzyme_combi + ("(" + enzymeComb[0].name +"-"+ enzymeComb[1].name +")")
										if fragDict.get(enzyme_combi, "missing") == "missing":
											fragDict[enzyme_combi] = 1
										else:
											fragDict[enzyme_combi] += 1

										#outFiles[outFileKey].write("@%s|frag_%i|%i-%i|%s-%s\n%s\n%s\n%s\n" % (contig.id,i,start.position,end.position+len(oHang),start.enzyme.name,end.enzyme.name,fragment,"+",qual))
										outFiles[outFileKey].write("@%s|frag_%i|%i-%i|%s-%s\n%s\n%s\n%s\n" % (contig.id,i,start.position,end.position,start.enzyme.name,end.enzyme.name,fragment,"+",qual))
									
										if sizeDICT.get(enzNames, "no")=="no":
											sizeDICT[enzNames]={}

										fragLen=len(fragment)
										if sizeDICT[enzNames].get(fragLen, "no")=="no":
											sizeDICT[enzNames][fragLen]=1
										else:
											sizeDICT[enzNames][fragLen]+=1
										
										if fragLen>maxFrag:
											maxFrag=fragLen

										i+=1

									else:
																		
										#add/count fragments in dictionary
										enzyme_combi = start.enzyme.name  + "_" + end.enzyme.name + ("(" + enzymeComb[0].name +"-"+ enzymeComb[1].name +")")
										enzNames=enzymeComb[0].name +"-"+ enzymeComb[1].name
										#if start.enzyme.name == "Start" or end.enzyme.name =="End":
										#enzyme_combi = enzyme_combi + ("(" + enzymeComb[0].name +"-"+ enzymeComb[1].name +")")
										if fragDict.get(enzyme_combi, "missing") == "missing":
											fragDict[enzyme_combi] = 1
										else:
											fragDict[enzyme_combi] += 1
										

										#outFiles[outFileKey].write(">%s|frag_%i|%i-%i|%s-%s\n%s\n" % (contig.id,i,start.position,end.position+len(oHang),start.enzyme.name,end.enzyme.name,fragment))
										outFiles[outFileKey].write(">%s|frag_%i|%i-%i|%s-%s\n%s\n" % (contig.id,i,start.position,end.position,start.enzyme.name,end.enzyme.name,fragment))
										
										if sizeDICT.get(enzNames, "no")=="no":
											sizeDICT[enzNames]={}

										fragLen=len(fragment)
										if sizeDICT[enzNames].get(fragLen, "no")=="no":
											sizeDICT[enzNames][fragLen]=1
										else:
											sizeDICT[enzNames][fragLen]+=1

										if fragLen>maxFrag:
											maxFrag=fragLen

										i+=1
#										print fragment	


			lastSites=unifiedSites

try:
	#output a file containing the number of fragments for each enzyme combination
	overviewWriter = open(args.o + "digestion_statistics.tsv", "w")
	

	for key in fragDict:
		overviewWriter.write("%s\t%i\n" % (key, fragDict[key]))


	for eComb in sizeDICT:
		oWriter = open(args.o + "%s_fragmentSizes.tsv"%(eComb), "w")
		for i in range(0, maxFrag+1):
		
			if sizeDICT[eComb].get(i, "nope") == "nope":
				oWriter.write("%i\t%i\n"%(i, 0))
			else:
				oWriter.write("%i\t%i\n"%(i, sizeDICT[eComb][i]))
	
except AttributeError:
	overviewWriter = sys.stdout

finally:
	overviewWriter.close()
	if args.v:
		print("All Done.") 
