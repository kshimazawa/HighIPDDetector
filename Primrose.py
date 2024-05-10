import pysam
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
def primrose(ALIGNED_FILE,RESULT_PATH):
    kinetics = []
    ml_list = []
    ml_cnt = 0
    bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)
    for read in bamfile.fetch(until_eof=True):
        mean_ipd_f = np.array(read.get_tag('fi'),dtype=int)
        mean_ipd_r = np.flip(np.array(read.get_tag('ri'),dtype=int))
        mean_pw_f = np.array(read.get_tag('fp'),dtype=int)
        mean_pw_r = np.flip(np.array(read.get_tag('rp'),dtype=int))
        if (len(mean_ipd_f)==0) | (len(mean_ipd_r)==0):
            print("no values for either forward/reverse strand")
            continue
        try:
            ml = np.array(read.get_tag('ML'),dtype=int)
        except:
            print("no ml tag")
            continue
        mm = np.array(read.get_tag('MM').strip(";").split(",")[1:],dtype=int)
        seq = read.get_forward_sequence()
        iter  = 0
        mm_iter = 0
        for i in range(len(seq)):
            if mm_iter == len(mm):
                break
            if (seq[i] != "C"):
                continue
            if iter == np.sum(mm[:mm_iter+1])+mm_iter:
                mm_iter += 1
                f_start_iter = max(i-5,0)
                f_end_iter = min(i+11,len(seq))
                r_start_iter = max(i-8,0)
                r_end_iter = min(i + 8,len(seq))
                kinetics_row = np.zeros((4,16))
                if f_start_iter == 0:
                    kinetics_row[0][5-i:] = mean_ipd_f[f_start_iter:f_end_iter]
                    kinetics_row[2][5-i:] = mean_pw_f[f_start_iter:f_end_iter]
                elif f_end_iter == len(seq):
                    kinetics_row[0][:(f_end_iter-f_start_iter)] = mean_ipd_f[f_start_iter:f_end_iter]
                    kinetics_row[2][:(f_end_iter-f_start_iter)] = mean_pw_f[f_start_iter:f_end_iter]
                else:
                    kinetics_row[0] = mean_ipd_f[f_start_iter:f_end_iter]
                    kinetics_row[2] = mean_pw_f[f_start_iter:f_end_iter]
                if r_start_iter == 0:
                    kinetics_row[1][8-i:] = mean_ipd_r[r_start_iter:r_end_iter]
                    kinetics_row[3][8-i:] = mean_pw_r[r_start_iter:r_end_iter]
                elif r_end_iter == len(seq):
                    kinetics_row[1][:(r_end_iter-r_start_iter)] = mean_ipd_r[r_start_iter:r_end_iter]
                    kinetics_row[3][:(r_end_iter-r_start_iter)] = mean_pw_r[r_start_iter:r_end_iter]
                else:
                    kinetics_row[1] = mean_ipd_r[r_start_iter:r_end_iter]
                    kinetics_row[3] = mean_pw_r[r_start_iter:r_end_iter]
                kinetics.append(kinetics_row)
            iter += 1
        ml_list.append(ml)
        ml_cnt += len(ml)
    bamfile.close()
    kinetics = np.array(kinetics)
    ml_flatten = [v for ml in ml_list for v in ml]
    category =["forward_ipd","reverse_ipd","forward_pw","reverse_pw"]
    for kin in range(4):
        for pos in range(16):
            plt.figure(figsize=(10,10))
            print(len(ml_flatten),kinetics.shape)
            plt.hist2d(kinetics[:,kin,pos],ml_flatten,bins=[np.linspace(0,256,257),np.linspace(0,256,257)],density=True,cmap=cm.jet)
            plt.savefig(RESULT_PATH+"{}_pos{}.png".format(category[kin],pos))
            plt.close()
    return

        
if __name__ == "__main__":
    primrose(sys.argv[1],sys.argv[2])
    
