# reads list of restriction enzymes from file
# expected format:
# enzyme1	pattern1[	enzyme2	pattern2]
# returns list of list with Biopython restriction enzymes
def readEnzymeListWithPattern(eFile):
	enzymes = []
	
	for line in eFile:
		data = line.strip().split("\t")
		thisEnzymes = []
		
		for i in range(0,len(data),2):
			enzymeName=data[i]
			try:
				enzymePattern=data[i+1]
				thisEnzymes.append(myRestrictionEnzyme(enzymeName,enzymePattern,i))
			except IndexError:
				print("ERROR: Wrong format for enzyme file. No restriction pattern found for enzyme %s." % enzymeName)
				print("       The current line will be ignored: %s" % line)
				thisEnzymes=None
				break
		
		if thisEnzymes: enzymes.append(thisEnzymes)
	
	return(enzymes)


# custom class for restriction enzymes
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Data.IUPACData import ambiguous_dna_values
import re
class myRestrictionEnzyme:
	# self.pattern: the original pattern the user provided (including marks for cutting sites; string)
	# self.site: the DNA sequence to search for; might contain ambiguities (Seq)
	# self.pPattern: parsed site with ambiguities transformed to regex character classes (string)
	#self.order --> given the order of the enzyme in the parsed file --> first enzyme(forward)=1, second(reverse)=2
	def __init__(self, name, pattern, order):
		self.name=name
		self.pattern=pattern
		self.parsePattern()
		self.order=order
	def __repr__(self):
		return("%s: %s" % (self.name, self.pattern))
	def translAmbs(self):
		self.pPattern="".join(["%s" % ambiguous_dna_values[x] for x in str(self.site)])
	def parsePattern(self):
		# Biopython:
		# ^ refers to the position of the cut in the sense strand of the sequence
		# _ to the cut on the antisense or complementary strand
		tempPattern=self.pattern.upper().replace("*","^_")
		#print tempPattern
		fwdSite=tempPattern.find("^")
		revSite=tempPattern.find("_")

		self.site=Seq.Seq(re.sub("[^A-Z]","",tempPattern), IUPAC.ambiguous_dna)
		#print self.site


		self.addBack=Seq.Seq(tempPattern[min(fwdSite+1,revSite+1):max(fwdSite,revSite)], IUPAC.ambiguous_dna)
		#variable to define if a rev site is before the fwd site --> overlap needs to be added in the front of a seq
		if revSite>=fwdSite:
			self.front=0
			self.offset=fwdSite
		else:
			self.front=1
			#-1 for offset because 
			self.offset=fwdSite-1
		#taken out reverse seqs
		#print(tempPattern, revSite, self.addBack)
		#self.offset=min(fwdSite,revSite)
		#self.offset=fwdSite
		self.translAmbs()


# custom class for cutting site
from Bio import Seq
class myCuttingSite():
	def __init__(self, position, enzyme, forward):
		self.position=position
		self.enzyme=enzyme
		self.forward=forward
	def completeOverhang(self):
		if self.forward: return(self.enzyme.addBack)
		else: return(self.enzyme.addBack)
	def addLeftOfCut(self):
		if self.forward: return(self.enzyme.addBack)
		else: return(Seq.Seq(""))
	def addRightOfCut(self):
		if self.forward: return(Seq.Seq(""))
		else: return(self.enzyme.addBack)
	def __repr__(self):
		revFlag=["*",""]
		return("%i-%s%s" % (self.position, self.enzyme.name, revFlag[self.forward]))

# if cutting site is symmetric, it will be found in fwd and rev sequence and
# thus produce two fragments that are identical due to completed overhang
# to prevent this, sites at the same position from the same enzyme will be unified
def unifySites(sites):
	unified = []
	for site in sites:
		if not siteInListOfSites(site, unified):
			unified.append(site)
	return(unified)
def siteInListOfSites(site, siteList):
	for siteFromList in siteList:
		if siteFromList.enzyme==site.enzyme and siteFromList.position==site.position:
			return(True)
	return(False)

import os
def mkOutDir(path):
	try:
		os.makedirs(os.path.dirname(path))
	except OSError:
		print("Output directory already exists. Already existing files will be overwritten.")
	if os.path.basename(path):
		return(path+"_")
	else: return(path)




def digest(genome, freq_set, rare_set):
    # perform in-silico restriction
    rb = RestrictionBatch([list(rare_set)[0], list(freq_set)[0]])
    #rb = RestrictionBatch(['EcoRI', 'Bsp143I'])
    f = open("%s-%s.frag.len.csv" % ("_".join([str(x) for x in freq_set]), "_".join([str(x) for x in rare_set])), "w")
    for rec in genome:
        # find cutting sites in contig
        ana = Analysis(rb, rec.seq)
        rare_sites = ana.__dict__['mapping'].values()[0]
        freq_sites = ana.__dict__['mapping'].values()[1]
        if len(rare_sites) == 0 or len(freq_sites) == 0:
            continue
        # transform cuttings sites into queues)
        sites = [deque(freq_sites), deque(rare_sites)]
        lead = sites[0][0] > sites[1][0]
        # find neighbouring cutting sites for both enzymes
        lengths = []
        while len(sites[0])>0 and len(sites[1])>0:
            last_popped = sites[lead].popleft()
            if len(sites[lead]) == 0:
                break
            while sites[lead][0]<sites[not lead][0]:
                last_popped = sites[lead].popleft()
                if len(sites[lead])==0:
                    break
            if len(sites[lead])>0:
                lengths.append(sites[not lead][0]-last_popped)
                lead = not lead
        #sites = sorted(ana.__dict__['mapping'].values()[0] + ana.__dict__['mapping'].values()[1])
        #f.write("\n".join([str(b-a) for a,b in zip([0]+sites[:-1],sites)]) + "\n")
        f.write("\n".join([str(x) for x in lengths]) + "\n")





