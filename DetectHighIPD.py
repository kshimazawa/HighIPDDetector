import numpy as np
import pysam
from concurrent.futures import ProcessPoolExecutor
import os
ALIGNED_FILE = "#.bam"
RESULT_PATH = 'result.txt'
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
    with ProcessPoolExecutor() as executor:
        for i in range(mat_count):
            future = executor.submit(detect,ipd_mat = ipd_mat_list[i])
            try:
                result_dict[zm_list[i]]=future.result()
            except:
                print(result_dict)
                exit(1)
    with open(RESULT_PATH,"a") as f:
        for k in result_dict:
            f.write("{}\t{}\n".format(k,result_dict[k]))
    print("until zm{} done".format(zm_list[-1]))
#make an ipd matrix                                                                                                                                                                          
def make_matrix(zm,align_result,ipd,min_pos,max_pos):
    ipd_matrix = np.zeros((len(align_result),max_pos-min_pos+1))
    for p in range(len(align_result)):
        for pos in range(align_result[p].shape[0]):
            iter = align_result[p][pos,1]-min_pos
            iter_query = align_result[p][pos,0]
            ipd_matrix[p,iter] = ipd[p][iter_query]
    return ipd_matrix
#read a bam file and extract subreads with the same hole number                                                                                                                                             
def align_by_read(ALIGNED_FILE):
    bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)    
    prev_zm = 0
    align_result = []#shape:(pass,matched_bases,2)                                                                                                                                                         
    ipd = []#shape:(pass,subread_length)                                                                                                                                                                   
    min_pos = 0
    max_pos = 0
    ipd_mat_list = []
    zm_list = [] 
    mat_count = 0
    for read in bamfile.fetch(until_eof=True):
        zm = int(read.get_tag('zm'))
        if (prev_zm != zm) & (prev_zm != 0):
            ipd_mat_list.append(make_matrix(prev_zm,align_result,ipd,min_pos,max_pos))
            mat_count += 1
            zm_list.append(zm)
            if mat_count >= 100:
                run_detect(ipd_mat_list,mat_count,zm_list)
                ipd_mat_list = []
                mat_count = 0
                zm_list = []
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
    ipd_mat_list.append(make_matrix(prev_zm,align_result,ipd,min_pos,max_pos))
    mat_count += 1
    zm_list.append(prev_zm)
    run_detect(ipd_mat_list,mat_count,zm_list)
    bamfile.close()
    return

if __name__=='__main__':
    align_by_read(ALIGNED_FILE)
