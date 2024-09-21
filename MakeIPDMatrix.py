import pysam
import numpy as np
MAX_READ_LEN = 100000
MAX_PASS = 1000
#read a bam file and extract subreads with the same hole number
def align_by_read(ALIGNED_FILE,FOLDER_PATH):
    bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)
    prev_zm = 0
    min_pos = 0
    max_pos = 0
    ipd_matrix = np.zeros((MAX_PASS,MAX_READ_LEN))
    pass_count = 0
    count = 0
    for read in bamfile.fetch(until_eof=True):
        zm = int(read.get_tag('zm'))
        if (prev_zm != zm) & (prev_zm != 0):
            np.savetxt("{}/zm{}_ipdmatrix.txt".format(FOLDER_PATH,zm),ipd_matrix[:pass_count,:max_pos+1],fmt='%d')
            count += 1
            if count >= 1000:
                return
            ipd_matrix = np.zeros((MAX_PASS,MAX_READ_LEN))
            min_pos = 0
            max_pos = 0
            pass_count = 0
        prev_zm = zm
        aligned_pairs = np.array(read.get_aligned_pairs(matches_only=True),dtype=int)#list of tuples whose elements are (read_pos, ref_pos) for matched bases
        raw_ipd = np.array(read.get_tag('ip'),dtype=int)
        for i in range(len(aligned_pairs)):
            read_pos = aligned_pairs[i,0]
            ref_pos = aligned_pairs[i,1]
            ipd_matrix[pass_count,ref_pos] = raw_ipd[read_pos]
        pass_count += 1
        min_pos_subread = np.min(aligned_pairs[:,1])
        max_pos_subread = np.max(aligned_pairs[:,1])
        if min_pos > min_pos_subread:
            min_pos = min_pos_subread
        if max_pos < max_pos_subread:
            max_pos = max_pos_subread
        max_pos_subread = np.max(aligned_pairs[:,1])
    np.savetxt("{}/zm{}_ipdmatrix.txt".format(FOLDER_PATH,zm),ipd_matrix,fmt='%d')
    bamfile.close()
        
if __name__ == "__main__":
    ALIGNED_FILE = input()
    FOLDER_PATH = input()
    align_by_read(ALIGNED_FILE,FOLDER_PATH)
