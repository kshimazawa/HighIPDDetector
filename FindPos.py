# %%
import pandas as pd
import numpy as np
import re
import pysam

# %%
RESULT_FILE  = input()
BAM_FILE = input()
OUTPUT_FOLDER = input()
reader = pd.read_table(RESULT_FILE,header=None,chunksize=10000)
nipd_ref_pos =  {}#key:ref_name, value:array whose shape is (ref_genome,2)
for df in reader:
    # %%
    chipd = {}
    for i in range(len(df)):
        cnt = 0
        m = re.search(r"(\d+)",df.iloc[i][0].astype(str))
        ml = re.findall(r"(\d+(\.\d+)?)",df.iloc[i][1])
        chipd[int(m.group())] = ml
    
    # %%
    bamfile = pysam.AlignmentFile(BAM_FILE,'rb',check_sq=False)
    for read in bamfile.fetch(until_eof=True):
        zm = int(read.get_tag('zm'))
        if zm not in chipd:
            continue
        tmp_ref_pos = []
        aligned_pairs = np.array(read.get_aligned_pairs(matches_only=True)).astype(int)#list of tuples whose elements are (read_pos, ref_pos) for matched bases
        ref_name = read.reference_name
        nipd_ref_pos.setdefault(ref_name,np.zeros((50000000,3)))
        for it in range(len(chipd[zm])):
            try:
                ref_pos = aligned_pairs[aligned_pairs[:,0]==it,1][0]
            except:
                continue
            nipd_ref_pos[ref_name][ref_pos,0] += float(chipd[zm][it][0])
            nipd_ref_pos[ref_name][ref_pos,1] += 1
            nipd_ref_pos[ref_name][ref_pos,2] = (nipd_ref_pos[ref_name][ref_pos,2]*max(1,nipd_ref_pos[ref_name][ref_pos,1]-1) + float(chipd[zm][it][0])**2)/nipd_ref_pos[ref_name][ref_pos,1]
            #print(nipd_ref_pos[ref_name][ref_pos])
    #print(chipd_ref_pos["{}_{}".format(zm,read.reference_name)])
    bamfile.close()
for k in nipd_ref_pos:
    np.savetxt(OUTPUT_FOLDER+k,nipd_ref_pos[k],fmt="%0.5f",delimiter=",")

                                                                                                                                                            
