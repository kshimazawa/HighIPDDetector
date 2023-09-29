import numpy as np
import re
FILE = "/work/00/gg57/share/genomes/ce11rel606/ce11rel606.fa"
#chr_list = ["chrI","chrII","chrIII","chrIV","chrV","chrX","chrM"]
def extract_seq(fname = FILE):
    pos = {}
    is_first = True
    with open(fname,"r") as fr:
        seq = ""
        key = ""
        for line in fr:
            if (">" in line) & is_first:
                key = line.replace(">","").strip("\n")
                pos[key] = []
                is_first = False
            elif ">" in line:
                iter_list = np.loadtxt("/work/gg57/g57016/data/{}".format(key)).astype(int)
                for i in iter_list:
                    pos[key].append(seq[i-20:i+21])
                key = line.replace(">","").strip("\n")
                pos[key] = []
                seq = ""
            else:
                seq = seq + line.strip("\n")
    with open("all.fasta","w") as fw:
        for k in pos:
#        with open("{}.fasta".format(k),'w') as fw:
            for i in range(len(pos[k])):
                fw.write(">{}_{}\n".format(k,i))
                fw.write(pos[k][i]+"\n")
    return

if __name__ == '__main__':
    extract_seq(FILE)
