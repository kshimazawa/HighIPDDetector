# %%
import pandas as pd
import numpy as np
import re
import pysam

# %%
RESULT_FILE = "/work/gg57/g57016/data/result.txt"
BAM_FILE = "/work/gg57/g57016/data/m64121_220824_041525.ccs.aligned.bam"
CHIPD_FILE = "/work/gg57/g57016/data/range_forward.txt"
# %%
df = pd.read_table(RESULT_FILE,header=None)
zm_list = []
for i in range(len(df)):
    cnt = 0
    m = re.search(r"(\d+)",df.iloc[i][0].astype(str))
    zm_list.append(int(m.group()))
zm_list = set(zm_list)
df = []
bamfile = pysam.AlignmentFile(BAM_FILE,'rb',check_sq=False)
range_pos = {}
for read in bamfile.fetch(until_eof=True):
    zm = int(read.get_tag('zm'))
    if zm not in zm_list:
        continue
    aligned_pairs = np.array(read.get_aligned_pairs(matches_only=True)).astype(int)#list of tuples whose elements are (read_pos, ref_pos) for matched bases
    min_pos = np.min(aligned_pairs[:,1])
    max_pos = np.max(aligned_pairs[:,1])
    key_name = "{}_{}".format(zm,read.reference_name)
    range_pos[key_name] = (min_pos,max_pos)
#print(chipd_ref_pos["{}_{}".format(zm,read.reference_name)])
bamfile.close()
with open(CHIPD_FILE,"w") as fw:
    for k in range_pos:
        fw.write("{}".format(k))
        for v in range_pos[k]:
            fw.write(",{}".format(v))
        fw.write("\n")
        

