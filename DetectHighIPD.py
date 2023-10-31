import numpy as np
import pysam
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor
import os
from functools import partial
ALIGNED_FILE = "/work/gg57/g57016/data/m64121_220824_041525.subreads.actc.subsample.forward.bam"
RESULT_PATH = '/work/gg57/g57016/data/result_sub.txt'
MAX_READ_LEN = 100000
MAX_PASS = 1000
#detect consistently high IPD based on the binominal statistical testing
def combination(n,k):
    if k < 0:
        exit(1)
    if k == 0:
        return 1
    else:
        ans = 1
        for i in range(k):
            ans *= n-i
            ans /=(i + 1)
    return ans

#def detect(FILE_NAME,prob=0.01,alpha=0.01,min_passes=10):
def detect(ipd_mat,prob=0.01,alpha=0.01,min_passes=10):
    detected = []
    #ipd_mat = np.loadtxt(FILE_NAME,dtype='int')
    if (ipd_mat.shape[0] < min_passes) | (ipd_mat.ndim == 1):
        return [-1]
    #print(ipd_mat.shape,FILE_NAME)
    rank_array = np.zeros(ipd_mat.shape)
    for j in range(rank_array.shape[0]):
        rank_array[j] = np.argsort(np.argsort(ipd_mat[j]))

    passes = rank_array.shape[0]
    for pos in range(rank_array.shape[1]):
        cnt = np.sum(rank_array[:,pos] > (rank_array.shape[1]*(1-prob)))
        pval = 0
        comb = 1
        if cnt != 0:
            comb = combination(passes,cnt-1)
        for i in range(cnt,passes+1):
            if i == 0:
                pval += comb*(prob**i)*((1-prob)**(passes-i))
            else:
                pval += comb*(passes-i+1)/i*(prob**i)*((1-prob)**(passes-i))
        if pval < alpha/rank_array.shape[1]:
            detected.append(pos)
    return detected

def run_detect(ipd_mat_list,mat_count,zm_list):
    result_dict = {}
    result = []
    with ProcessPoolExecutor() as executor:
        result = list(executor.map(detect,ipd_mat_list))
    for i in range(mat_count):
        result_dict[zm_list[i]]=result[i]
    with open(RESULT_PATH,"a") as f:
        for k in result_dict:
            f.write("{}\t{}\n".format(k,result_dict[k]))
    print("until zm{} done".format(zm_list[-1]))

def align_by_read(ALIGNED_FILE,FOLDER_PATH=""):
    bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)
    prev_zm = 0
    min_pos = 0
    max_pos = 0
    ipd_mat_list = []
    zm_list = [] 
    mat_count = 0
    ipd_matrix = np.zeros((MAX_PASS,MAX_READ_LEN))
    pass_count = 0
    for read in bamfile.fetch(until_eof=True):
        zm = int(read.get_tag('zm'))
        if (prev_zm != zm) & (prev_zm != 0):
            if FOLDER_PATH != "":
                np.savetxt("{}/zm{}_ipdmatrix.txt".format(FOLDER_PATH,zm),ipd_matrix[:pass_count,:max_pos+1],fmt='%d')
            ipd_mat_list.append(ipd_matrix[:pass_count,:max_pos+1].copy())
            mat_count += 1
            zm_list.append(prev_zm)
            if mat_count >= 100:
                run_detect(ipd_mat_list,mat_count,zm_list)
                ipd_mat_list = []
                mat_count = 0
                zm_list = []
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
    if FOLDER_PATH != "":
        np.savetxt("{}/zm{}_ipdmatrix.txt".format(FOLDER_PATH,zm),ipd_matrix,fmt='%d')
    ipd_mat_list.append(ipd_matrix[:pass_count,:max_pos+1].copy())
    mat_count += 1
    zm_list.append(prev_zm)
    run_detect(ipd_mat_list,mat_count,zm_list)
    bamfile.close()
    return

if __name__ == "__main__":
    align_by_read(ALIGNED_FILE)
