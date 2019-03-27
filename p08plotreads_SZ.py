import pysam
import sys,os
import numpy as np
from collections import namedtuple
from collections import	Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.collections as collections


usage = """\
#============================================================================#
python realization of Murim's p08_plot_reads

currently only support hg19 and GATK's bam files.
quality score filter is 6 in GATK's bam file scale

Usage: python %s in.bam chr position sampleID 
e.g. python %s in.bam 17 8065350 sampleID

TO-DO: need to work on CIGAR strings!!! (siming---4.1.2015)
#============================================================================#
"""
def arrayinfo(pileupcolumn,reference,qsfilter=6):
	Array= namedtuple('Array', ['mutqpos','refseq','qualarray','ifrefarray','basearray'])
	pileups=pileupcolumn.pileups
	leftmostqpos=0
	Array.mutqpos=pileups[0].query_position
	rightmostqpos=max([Array.mutqpos-pileups[i].query_position+len(pileups[i].alignment.seq) for i in range(pileupcolumn.n)])
	Array.refseq=refnt(reference=reference,start=pileupcolumn.pos-Array.mutqpos,end=pileupcolumn.pos-Array.mutqpos+rightmostqpos,reffastafile=REFFASTA)
	Array.qualarray=np.empty([pileupcolumn.n,rightmostqpos]) # Array.qualarray records the quality for each pos vertically from the current pilupreads
	Array.qualarray[:]=np.NAN
	Array.ifrefarray=np.empty([pileupcolumn.n,rightmostqpos]) # Array.ifrefarray records if the called base is refbase for each pos vertically from the current pilupreads
	Array.ifrefarray[:]=np.NAN
	Array.basearray=np.empty([pileupcolumn.n,rightmostqpos],dtype="S1") # Array.basearray is a char array records the nucleotide (ACGT) type
	Array.basearray[:]='N'
	n=0
	for pileupread in pileups:
		if ord(pileupread.alignment.qual[pileupread.query_position])-33< qsfilter:
			continue
		for i in range(len(pileupread.alignment.seq)):
			#print(n,i)
			Array.qualarray[n,Array.mutqpos-pileupread.query_position+i]=ord(pileupread.alignment.qual[i])-33
			Array.ifrefarray[n,Array.mutqpos-pileupread.query_position+i]=1 if pileupread.alignment.seq[i]==Array.refseq[Array.mutqpos-pileupread.query_position+i] else 0
			Array.basearray[n,Array.mutqpos-pileupread.query_position+i]=pileupread.alignment.seq[i]
		n=n+1
	return Array
	
def refnt(reference=None,start=None,end=None, reffastafile=None):
	'''start and end, use python index'''
	reffile=pysam.Fastafile(reffastafile) 
	refs=reffile.fetch(reference,start,end).upper()
	return refs

def plot_pileupcolumn(Array,figurename): # need array object as obtained from arrayinfo function 
	highestqual=38
	colordict={'A':'red','C':'blue','T':'yellow','G':'green'}

	refqual=np.repeat(highestqual,len(Array.refseq))
	
	ntplotarray=np.vstack([np.array(list(Array.refseq)),Array.basearray])
	qualplotarray=np.vstack([refqual,Array.qualarray])
	np.savetxt('%s.txt'%figurename,ntplotarray,fmt='%s',delimiter='')
	
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xlim(-1,ntplotarray.shape[1])
	ax.set_ylim(-ntplotarray.shape[0],1)
    
	patches = {'A':[],'T':[],'C':[],'G':[]}

	for x in range(ntplotarray.shape[0]):
		for y in range(ntplotarray.shape[1]):
			if ntplotarray[x,y]!='N':
				rect = mpatches.Rectangle((y,-x),qualplotarray[x,y]/highestqual,1)
				patches[ntplotarray[x,y]].append(rect)
			
 	print("plotting %s..."%figurename)
 	for ntkey in patches.keys():	
	 	mycollection = collections.PatchCollection(patches[ntkey],facecolor=colordict[ntkey],edgecolor='none')
		ax.add_collection(mycollection)
		
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none') 
	ax.spines['left'].set_color('none')
	ax.spines['bottom'].set_color('none')  
	ax.xaxis.set_ticks_position('top')
	ax.yaxis.set_ticks_position('left')
	ax.tick_params(axis='x', direction='out',length=3) 
	ax.set_xticks([Array.mutqpos+0.5])
	ax.set_xticklabels(['Query'])
	ax.set_facecolor('black')
	fig.set_size_inches(8,ntplotarray.shape[0]*0.057)
	plt.savefig('%s.png'%figurename,dpi=300, bbox_inches='tight')
	
if __name__=='__main__':
	if len(sys.argv) !=5:
		sys.exit(usage % (os.path.basename(sys.argv[0]),os.path.basename(sys.argv[0])))
	else:
		bamfile=sys.argv[1]
		chr=sys.argv[2]
		start=int(sys.argv[3])-1
		end=start+1
		figurename=sys.argv[4]+'_chr'+chr+'_'+str(start+1)+'_plotreads'
		REFFASTA='/mnt/d/projects/Unix/GRCh38.d1.vd1.fa'
		
		samfile = pysam.Samfile(bamfile, "rb" )
		for pileupcolumn in samfile.pileup(chr, start, end):
			if pileupcolumn.pos==start:
				Array=arrayinfo(pileupcolumn,chr,qsfilter=6)
				print('Coverage for this position is %d...'%Array.basearray.shape[0])
				plot_pileupcolumn(Array,figurename)
	
