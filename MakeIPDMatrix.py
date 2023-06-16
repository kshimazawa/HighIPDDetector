import pysam
import numpy as np
#make an ipd matrix
def make_matrix(zm,align_result,ipd,min_pos,max_pos,folder_pass):
    ipd_matrix = np.zeros((len(align_result),max_pos-min_pos+1))
    for p in range(len(align_result)):
        for pos in range(align_result[p].shape[0]):
            iter = align_result[p][pos,1]-min_pos
            iter_query = align_result[p][pos,0]
            ipd_matrix[p,iter] = ipd[p][iter_query]
    np.savetxt("{}/zm{}_ipdmatrix.txt".format(folder_pass,zm),ipd_matrix,fmt='%d')
#read a bam file and extract subreads with the same hole number
def align_by_read(ALIGNED_FILE,FOLDER_PASS):
    bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)
    prev_zm = 0
    align_result = []#shape:(pass,matched_bases,2)
    ipd = []#shape:(pass,subread_length)
    min_pos = 0
    max_pos = 0
    for read in bamfile.fetch(until_eof=True):
        zm = int(read.get_tag('zm'))
        if (prev_zm != zm) & (prev_zm != 0):
            make_matrix(prev_zm,align_result,ipd,min_pos,max_pos,FOLDER_PASS)
            align_result = []
            ipd = []
            min_pos = 0
            max_pos = 0
        prev_zm = zm
        aligned_pairs = np.array(read.get_aligned_pairs(matches_only=True))#list of tuples whose elements are (read_pos, ref_pos) for matched bases
        min_pos_subread = np.min(aligned_pairs[:,1])
        max_pos_subread = np.max(aligned_pairs[:,1])
        if min_pos > min_pos_subread:
            min_pos = min_pos_subread
        if max_pos < max_pos_subread:
            max_pos = max_pos_subread
        max_pos_subread = np.max(aligned_pairs[:,1])
        align_result.append(aligned_pairs)
        ipd.append(read.get_tag('ip'))
    make_matrix(prev_zm,align_result,ipd)
    bamfile.close()
        
