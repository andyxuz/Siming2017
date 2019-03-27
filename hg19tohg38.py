import pyliftover
import sys
if __name__=='__main__':
    chrom = str("chr"+sys.argv[1])  
    pos = int(sys.argv[2])
    lo = pyliftover.LiftOver('hg19ToHg38.over.chain')
    result = lo.convert_coordinate(chrom,pos)
    result = str(result[0]).replace('(','').replace(')','').replace("'","")
    result = [i.strip() for i in result.split(',')]
    print result[0]+','+result[1]
    
