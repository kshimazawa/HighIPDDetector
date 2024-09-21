import numpy as np
import pysam
from concurrent.futures import ProcessPoolExecutor

class NIPD:
    def __init__(self,ALIGNED_FILE,CCS_FILE,RESULT_PATH,MAX_READ_LEN,MAX_PASS,REF_FASTA):
        self.ALIGNED_FILE = ALIGNED_FILE
        self.CCS_FILE = CCS_FILE
        self.RESULT_PATH = RESULT_PATH
        self.MAX_READ_LEN = MAX_READ_LEN
        self.MAX_PASS = MAX_PASS
        self.REF_FASTA = REF_FASTA
        # self.stats = {}
        # ref_types = {}
        # with open(self.REF_FASTA,"r") as fr:
        #     ref_name = ""
        #     for line in fr:
        #         if ">" in line:
        #             ref_name = line.replace(">","").strip("\n")
        #             ref_types[ref_name] = 0
        #         else:
        #             ref_types[ref_name] += len(line.strip("\n"))
        # for k in ref_types:
        #     self.stats[k] = np.zeros((3,ref_types[k]))#smean,mean,num
    def calcNIPD(self,ipd_mat):
        # if (ipd_mat.shape[0] < min_passes) | (ipd_mat.ndim == 1):
        #    return
        n_ipd = np.zeros(ipd_mat.shape)
        if ipd_mat.ndim == 1:
            n_ipd = np.reshape(n_ipd,(1,ipd_mat.shape[0]))
        if ("reverse" in self.ALIGNED_FILE):#if the ipd values derive from forward strands
            for j in range(n_ipd.shape[0]):
                for pos in range(n_ipd.shape[1]):
                    if (ipd_mat[j,pos] == 0) | (ipd_mat[j,pos-1] == 0) | (pos == 0):
                        continue
                    n_ipd[j,pos] = ipd_mat[j,pos]*2/(ipd_mat[j,pos]+ipd_mat[j,pos-1])
        else:#if the ipd values derive from reverse strands
            for j in range(n_ipd.shape[0]):
                for pos in reversed(range(n_ipd.shape[1])):
                    if (ipd_mat[j,pos] == 0) | (ipd_mat[j,pos-1] == 0) | (pos == n_ipd.shape[1]-1):
                        continue
                    n_ipd[j,pos] = ipd_mat[j,pos]*2/(ipd_mat[j,pos]+ipd_mat[j,pos+1])
        return n_ipd
    # def map(self,zm_pass,findex,min_passes=5):
    #     zm_set = set([zm_pass[i][0] for i in range(len(zm_pass))])
    #     if (n_ipd.shape[0] < min_passes):
    #         return
    #     bamfile = pysam.AlignmentFile(self.CCS_FILE,'rb',check_sq=False)
    #     for read in bamfile.fetch(until_eof=True):
    #         zm = int(read.get_tag('zm'))
    #         if zm not in zm_set:
    #             continue
    #         if (("reverse" in self.ALIGNED_FILE) & (read.is_reverse)) | ("forward" in self.ALIGNED_FILE) & (read.is_reverse==False):
    #             continue
    #         line_to_skip = 0
    #         for i in range(len(zm_pass)):
    #             if zm == zm_pass[i][0]:
    #                 break
    #             line_to_skip += zm_pass[i][1]
    #         with open(self.RESULT_PATH+"tmp.{}.bam.nipd.csv".format(findex),"r") as fr:
    #             iter = 0
    #             length = 0
    #             for line in fr:
    #                 iter += 1
    #                 if iter < line_to_skip:
    #                     continue
    #                 row = line.split(',')
    #                 if length != len(row):
    #                     break
    #                 length = len(row)
                    

    #         if read.is_reverse:
    #             n_ipd = np.flip(n_ipd,axis=1)
    #         aligned_pairs = np.array(read.get_aligned_pairs(matches_only=True),dtype=int)
    #         ref_name = read.reference_name
    #         for i in range(aligned_pairs.shape[0]):
    #             read_pos = aligned_pairs[i,0]
    #             ref_pos = aligned_pairs[i,1]
    #             n_ipd_nonzero = n_ipd[:,read_pos].copy()
    #             n_ipd_nonzero = n_ipd_nonzero[n_ipd_nonzero>0]
    #             smean = np.mean(n_ipd_nonzero**n_ipd_nonzero)
    #             ave = np.mean(n_ipd_nonzero)
    #             cnt = n_ipd_nonzero.shape[0]
    #             if cnt == 0:
    #                 continue
    #             if self.stats[ref_name][2,ref_pos] == 0:
    #                 self.stats[ref_name][0,ref_pos] = smean
    #                 self.stats[ref_name][1,ref_pos] = ave
    #                 self.stats[ref_name][2,ref_pos] = cnt
    #             else:
    #                 self.stats[ref_name][0,ref_pos] = (self.stats[ref_name][0,ref_pos]*self.stats[ref_name][2,ref_pos]+smean*cnt)/(self.stats[ref_name][2,ref_pos]+cnt)
    #                 self.stats[ref_name][1,ref_pos] = (self.stats[ref_name][1,ref_pos]*self.stats[ref_name][2,ref_pos]+ave*cnt)/(self.stats[ref_name][2,ref_pos]+cnt)
    #                 self.stats[ref_name][2,ref_pos] = self.stats[ref_name][2,ref_pos] + cnt
    #         return

    def splitbam(self,by=10000):
        bamfile = pysam.AlignmentFile(self.ALIGNED_FILE,'rb',check_sq=False)
        read_cnt = 0
        num_files = 1
        split_file = pysam.AlignmentFile(self.RESULT_PATH+"tmp.{}.bam".format(read_cnt//by),'wb',template=bamfile)
        for read in bamfile.fetch(until_eof=True):
            if (read_cnt % by == 0) & (read_cnt != 0):
                split_file.close()
                split_file = pysam.AlignmentFile(self.RESULT_PATH+"tmp.{}.bam".format(read_cnt//by),'wb',template=bamfile)
                num_files += 1
            split_file.write(read)
            read_cnt += 1
        split_file.close()
        bamfile.close()
        return num_files
    
    def write_ipdmat(self,ipd_mat_list,file_name):
        with open(file_name,"w") as fw:
            for k in ipd_mat_list:
                np.savetxt(fw,ipd_mat_list[k],fmt='%.5f')
        return
            
    def align_by_read(self,FILE_PATH):
        bamfile = pysam.AlignmentFile(FILE_PATH,'rb',check_sq=False)
        prev_zm = 0
        max_pos = 0
        nipd_list = {}
        zm_list = []
        pass_list = []
        mat_count = 0
        ipd_matrix = np.zeros((self.MAX_PASS,self.MAX_READ_LEN))
        pass_count = 0
        for read in bamfile.fetch(until_eof=True):
            zm = int(read.get_tag('zm'))
            if (prev_zm != zm) & (prev_zm != 0):
                nipd = self.calcNIPD(ipd_matrix[:pass_count,:max_pos+1])
                nipd_list[mat_count] = nipd[:pass_count,:max_pos+1]
                mat_count += 1
                zm_list.append(prev_zm)
                pass_list.append(pass_count)
                ipd_matrix = np.zeros((self.MAX_PASS,self.MAX_READ_LEN))
                max_pos = 0
                pass_count = 0
            prev_zm = zm
            aligned_pairs = np.array(read.get_aligned_pairs(matches_only=True),dtype=int)#list of tuples whose elements are (read_pos, ref_pos) for matched bases
            raw_ipd = np.array(read.get_tag('ip'),dtype=int)
            if "reverse" in self.ALIGNED_FILE:
                raw_ipd = np.flip(raw_ipd)
            for i in range(len(aligned_pairs)):
                read_pos = aligned_pairs[i,0]
                ref_pos = aligned_pairs[i,1]
                ipd_matrix[pass_count,ref_pos] = raw_ipd[read_pos]
            pass_count += 1
            max_pos_subread = np.max(aligned_pairs[:,1])
            if max_pos < max_pos_subread:
                max_pos = max_pos_subread
        bamfile.close()
        nipd_list.append(ipd_matrix[:pass_count,:max_pos+1].copy())
        mat_count += 1
        zm_list.append(prev_zm)
        pass_list.append(pass_count)
        self.write_ipdmat(nipd_list,self.RESULT_PATH+"{}.nipd.csv".format(FILE_PATH)) 
        zm_pass = [(zm_list[i],pass_list[i]) for i in range(len(zm_list))]
        return zm_pass

    def __call__(self):
        num_files = self.splitbam()
        print("file split done")
        tmpfile_list = [self.RESULT_PATH+"tmp.{}.bam".format(i) for i in range(num_files)]
        with ProcessPoolExecutor() as executor:
            zm_pass = executor.map(self.align_by_read,tmpfile_list)
        # with ProcessPoolExecutor() as executor:
        #     executor.map(self.map,zm_pass,[i for i in range(num_files)])
        return


if __name__ == "__main__":
    ALIGNED_FILE=input()
    CCS_FILE = input()
    RESULT_PATH=input()
    MAX_READ_LEN=input()
    MAX_PASS=input()
    REF_FASTA=input()
    nipd = NIPD(ALIGNED_FILE,CCS_FILE,RESULT_PATH,int(MAX_READ_LEN),int(MAX_PASS),REF_FASTA)
    nipd()
