import numpy as np
import pysam
from concurrent.futures import ProcessPoolExecutor

class DetectHighIPD:
    def __init__(self,ALIGNED_FILE,RESULT_PATH,MAX_READ_LEN,MAX_PASS,prob,alpha):
        self.ALIGNED_FILE = ALIGNED_FILE
        self.RESULT_PATH = RESULT_PATH
        self.MAX_READ_LEN = MAX_READ_LEN
        self.MAX_PASS = MAX_PASS
        self.prob = prob
        self.alpha = alpha
    #detect consistently high IPD based on the binominal statistical testing
    def combination(self,n,k):
        if (n <= 0) | (k < 0) | (n < k):
            exit(1)
        elif k == 0:
            return 1
        else:
            ans = 1
            for i in range(k):
                ans *= n-i
                ans /=(i + 1)
        return ans

    def detect(self,ipd_mat,min_passes=10):
        detected = []
        if (ipd_mat.shape[0] < min_passes) | (ipd_mat.ndim == 1):
            return [-1]
        rank_array = np.zeros(ipd_mat.shape)
        for j in range(rank_array.shape[0]):
            rank_array[j] = np.argsort(np.argsort(ipd_mat[j]))

        passes = rank_array.shape[0]
        for pos in range(rank_array.shape[1]):
            cnt = np.sum(rank_array[:,pos] > (rank_array.shape[1]*(1-self.prob)))
            pval = 0
            comb = 1
            if cnt != 0:
                comb = self.combination(passes,cnt-1)
            for i in range(cnt,passes+1):
                if i == 0:
                    pval += comb*(self.prob**i)*((1-self.prob)**(passes-i))
                else:
                    pval += comb*(passes-i+1)/i*(self.prob**i)*((1-self.prob)**(passes-i))
            if pval < self.alpha/rank_array.shape[1]:
                detected.append(pos)
        return detected

    def run_detect(self,ipd_mat_list,mat_count,zm_list):
        result_dict = {}
        result = []
        with ProcessPoolExecutor() as executor:
            result = list(executor.map(self.detect,ipd_mat_list[:mat_count]))
        for i in range(mat_count):
            result_dict[zm_list[i]]=result[i]
        with open(self.RESULT_PATH,"a") as f:
            for k in result_dict:
                f.write("{}\t{}\n".format(k,result_dict[k]))
                print("until zm{} done".format(zm_list[-1]))

    def align_by_read(self,FOLDER_PATH=""):
        bamfile = pysam.AlignmentFile(self.ALIGNED_FILE,'rb',check_sq=False)
        prev_zm = 0
        max_pos = 0
        ipd_mat_list = [None for i in range(100)]
        zm_list = [] 
        mat_count = 0
        ipd_matrix = np.zeros((self.MAX_PASS,self.MAX_READ_LEN))
        pass_count = 0
        for read in bamfile.fetch(until_eof=True):
            zm = int(read.get_tag('zm'))
            if (prev_zm != zm) & (prev_zm != 0):
                if FOLDER_PATH != "":
                    np.savetxt("{}/zm{}_ipdmatrix.txt".format(FOLDER_PATH,zm),ipd_matrix[:pass_count,:max_pos+1],fmt='%d')
                    ipd_mat_list[mat_count] = ipd_matrix[:pass_count,:max_pos+1]
                    mat_count += 1
                    zm_list.append(prev_zm)
                if mat_count >= 100:
                    self.run_detect(ipd_mat_list,mat_count,zm_list)
                    mat_count = 0
                    zm_list = []
                    ipd_matrix = np.zeros((self.MAX_PASS,self.MAX_READ_LEN))
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
            max_pos_subread = np.max(aligned_pairs[:,1])
            if max_pos < max_pos_subread:
                max_pos = max_pos_subread
        if FOLDER_PATH != "":
            np.savetxt("{}/zm{}_ipdmatrix.txt".format(FOLDER_PATH,zm),ipd_matrix,fmt='%d')
            ipd_mat_list.append(ipd_matrix[:pass_count,:max_pos+1].copy())
            mat_count += 1
            zm_list.append(prev_zm)
            self.run_detect(ipd_mat_list,mat_count,zm_list)
            bamfile.close()
        return



if __name__ == "__main__":
    ALIGNED_FILE,RESULT_PATH,MAX_READ_LEN,MAX_PASS,prob,alpha=input().split()
    dhipd = DetectHighIPD(ALIGNED_FILE,RESULT_PATH,MAX_READ_LEN,MAX_PASS,prob,alpha)
    dhipd.align_by_read()