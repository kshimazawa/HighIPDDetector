import re
import numpy as np
import pysam
import matplotlib.pyplot as plt
import sys

def find_motif(fname,motif):
    pos = {}
    is_first = True
    bamfile = pysam.AlignmentFile(fname,'rb',check_sq=False)
    for read in bamfile.fetch(until_eof=True):
        if "aligned" in fname:
            if read.is_reverse == False:
                continue
        zm = int(read.get_tag('zm'))
        seq = read.get_forward_sequence()
        pos[zm] = [m.span() for m in re.finditer(motif,seq)]
    bamfile.close()
    #with open(fname,"r") as fr:
    #    seq = ""
    #    key = ""
    #    for line in fr:
    #        if (">" in line) & is_first:
    #            key = line.replace(">","").strip("\n")
    #            pos[key] = []
    #            is_first = False
    #        elif ">" in line:
    #            pos[key] = [m.span() for m in re.finditer(MOTIF, seq)]
    #            key = line.replace(">","").strip("\n")
    #            pos[key] = []
    #            seq = ""
    #        else:
    #            seq = seq + line.strip("\n")
    #    pos[key] = [m.span() for m in re.finditer(MOTIF, seq)]
    return pos

def detected(result_file,index,bam,folder,motif):
    detected_pos = {}
    with open(result_file,"r") as fr:
        for line in fr:
            line_list = line.split()
            if ('[]' in line_list[1]):
                detected_pos[int(line_list[0])] = []
                continue
            detected = []
            for l in range(1,len(line_list)):
                if ("[" in line_list[l]) | ("]" in line_list[l]):
                    line_list[l] =re.sub("[\[\]]","",line_list[l])
                if "nan" in line_list[l]:
                    detected.append(0)
                else:
                    detected.append(float(line_list[l].strip("\n").strip(",")))
            detected_pos[int(line_list[0])] = detected
    motif_pos = find_motif(bam,motif)
    nipds = []
    med = 0
    mean = 0
    var = 0
    for k in motif_pos:
     #   m = re.search(r"/(\d+)/",k)
     #   zm = int(m.group(1))
        zm = k
        if zm not in detected_pos:
            continue
        for p in range(len(motif_pos[k])):
            pp = int(motif_pos[k][p][0]) + index
            try:
                if detected_pos[zm][pp] == 0:
                    continue
            except:
                continue
            nipds.append(detected_pos[zm][pp])
    med = np.median(nipds)
    mean = np.mean(nipds)
    var = np.var(nipds)
    plt.hist(nipds,bins=np.arange(0,max(nipds),0.1))
    plt.savefig("{}{}{}.png".format(folder,motif,index))
    plt.close()
    return med,mean,var
    
if __name__ == '__main__':
    RESULT_FILE = input()
    BAM_FILE = input()
    #FASTA_FILE = input()
    RESULT_FOLDER = input()
    MOTIF = input()
    index = int(input())
    med,mean,var = detected(RESULT_FILE,index,BAM_FILE,RESULT_FOLDER,MOTIF)
    print(med, mean, var)
