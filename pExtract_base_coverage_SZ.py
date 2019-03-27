
'''
Created on Mar 30, 2015

@author: siming

'''

import pysam
import sys,os
import numpy as np
from collections import namedtuple
from collections import    Counter
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.collections as collections


usage = """\
#============================================================================#

Extract base coverage

currently only support hg19 and GATK's bam files.
quality score filter is 6 in GATK's bam file scale

Usage: python %s inputfile in.bam sampleID 
inputfile has the following format:
chr\tpos\trefallele\tnonrefallele

TO-DO: need to work on CIGAR strings!!!

#============================================================================#
"""
def pileinfo(pileupcolumn,reference,qsfilter=6):
    pileups=pileupcolumn.pileups
    aqpos=pileupcolumn.pos
    myseq=[]
    for pileupread in pileups:
        readpos=aqpos-(pileupread.alignment.aend-pileupread.alignment.alen)
        if readpos<len(pileupread.alignment.seq):
            if ord(pileupread.alignment.qual[readpos])-33>= qsfilter:
                myseq.append(pileupread.alignment.seq[readpos])
    return myseq

    
def refnt(reference=None,start=None,end=None, reffastafile=None):
    '''start and end, use python index'''
    reffile=pysam.Fastafile(reffastafile) 
    refs=reffile.fetch(reference,start,end).upper()
    return refs

    
if __name__=='__main__':
    if len(sys.argv) !=4:
        sys.exit(usage % (os.path.basename(sys.argv[0]),os.path.basename(sys.argv[0])))
    else:
        bamfile=sys.argv[2]
        input=sys.argv[1]
        output=open(sys.argv[3]+'_cov.txt','w')
        REFFASTA='/mnt/d/projects/Unix/Homo_sapiens_assembly19.fasta'
        
        samfile = pysam.Samfile(bamfile, "rb" )
        inputfile=open(input,'r')
        for line in inputfile.readlines():
            item=line.strip('\n').split('\t')
            chr=item[0]
            start=int(item[1])-1
            end=start+1
            ref=item[2]
            nonref=item[3]
            for pileupcolumn in samfile.pileup(chr, start, end):
                if pileupcolumn.pos==start:
                    myseq=pileinfo(pileupcolumn,chr,qsfilter=6)
                    refcount=myseq.count(ref)
                    nonrefcount=myseq.count(nonref)
            #print('\t'.join([chr,str(start+1),ref,nonref,str(refcount),str(nonrefcount)])+'\n')
            output.write('\t'.join([chr,str(start+1),ref,nonref,str(refcount),str(nonrefcount)])+'\n')
        output.close()
                
            
    
