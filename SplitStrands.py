import pysam
import numpy as np
import sys
#make an ipd matrix
def split_strands(ALIGNED_FILE,NEW_FILE):
    bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)
    forward_bamfile = pysam.AlignmentFile("{}.forward.bam".format(NEW_FILE),'wb',template=bamfile)
    reverse_bamfile = pysam.AlignmentFile("{}.reverse.bam".format(NEW_FILE),'wb',template=bamfile)
    for read in bamfile.fetch(until_eof=True):
        if read.is_reverse:
            reverse_bamfile.write(read)
        else:
            forward_bamfile.write(read)
    forward_bamfile.close()
    reverse_bamfile.close()
    bamfile.close()
        
if __name__ == "__main__":
    split_strands("{}.subreads.actc.bam".format(sys.argv[1]),"{}.subreads.actc".format(sys.argv[1]))
    
