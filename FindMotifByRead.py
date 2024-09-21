import re
import pysam
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
FILE = "/nfs/data05/kismzw/data/SRR15720539/ER3413.subreads.actc.fasta"
MOTIF = ['GATC','GTATA', 'TATA', '[Tt][Gg][Tt][Tt][Cc][Aa][ACGTacgt][ACGTacgt][Cc]','TGGTACC[ACGT]C']
IPD_POSITION = [2,5,4,7,7]
def find_motif(motif,fname=FILE):
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
                pos[key] = [m.span() for m in re.finditer(motif, seq)]
                key = line.replace(">","").strip("\n")
                pos[key] = []
                seq = ""
            else:
                seq = seq + line.strip("\n")
        pos[key] = [m.span() for m in re.finditer(motif, seq)]
#    for k in pos:
#        print(k)
#        if len(pos[k])==0:
#            print(-1)
#            continue
#        for i in range(len(pos[k])):
#            print(pos[k][i][0])
    return pos

def detected_ratio(result_file):
    detected_pos = {}
    detection_hmap = np.zeros((len(MOTIF),20))
    with open(result_file,"r") as fr:
        for line in fr:
            line_list = line.split()
            if ('[]' in line_list[1]) | ('[-1]' in line_list[1]):
                detected_pos[int(line_list[0])] = []
                continue
            detected = []
            for i in range(1,len(line_list)):
                m = re.search(r"(\d+)",line_list[i])
                detected.append(int(m.group(1)))
            detected_pos[int(line_list[0])] = set(detected)
    for i in range(len(MOTIF)):
        motif_pos = find_motif(motif=MOTIF[i])
        total = 0
        num_detected = 0
        for k in motif_pos:
            if len(motif_pos[k]) == 0:
                continue
            m = re.search(r"/(\d+)/",k)
            zm = int(m.group(1))
            if zm in detected_pos:
                total += len(motif_pos[k])
                for p in motif_pos[k]:
                    for pp in range(detection_hmap.shape[1]):
                        if (p[0]+IPD_POSITION[i]-10 + pp) in detected_pos[zm]:
                            detection_hmap[i,-pp-1] += 1
        detection_hmap[i] = detection_hmap[i]/total
        print(detection_hmap[i])
#                    print("detected",zm,p[0]+1)
#                else:
#                    print("not detected",zm,p[0]+2,detected_pos[zm])
    return detection_hmap


if __name__ == '__main__':
    RESULT_FILE = sys.argv[1]
    detection_hmap = detected_ratio(RESULT_FILE)
    yticks = ['GATC','NTATAC', 'NTATA','GNNTGAACA','GNGGTACCA']
    xticks = np.arange(-9,11)
    plt.figure(figsize=(10,5))
    sns.heatmap(detection_hmap,yticklabels=yticks,xticklabels=xticks,cmap='Oranges')
    plt.title("Detection Rates of Methylation-free Sample")
    plt.savefig("detection_ratio2.png")
    #print(ratio,num_detected,total)
