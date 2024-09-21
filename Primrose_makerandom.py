import pysam
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import array as arr
import random
import string

def generate_from_regex(pattern):
    result = []
    i = 0
    while i < len(pattern):
        char = pattern[i]
        if char == '\\':
            # エスケープ文字の処理
            if i + 1 < len(pattern):
                next_char = pattern[i + 1]
                result.append(next_char)
                i += 1
        elif char == '.':
            # 任意の一文字
            result.append(random.choice(string.ascii_letters + string.digits))
        elif char == '[':
            # 文字クラス
            char_class = []
            i += 1
            if pattern[i] == '^':
                i += 1
            while pattern[i] != ']':
                char_class.append(pattern[i])
                i += 1
            result.append(random.choice(char_class))
        elif char == '{':
            # 繰り返し
            num = ''
            i += 1
            while pattern[i] != '}':
                num += pattern[i]
                i += 1
            count = int(num)
            result.extend([result[-1]] * (count - 1))
        elif char == '(':
            # グループ (基本的なサポート)
            sub_pattern = ''
            i += 1
            while pattern[i] != ')':
                sub_pattern += pattern[i]
                i += 1
            result.append(generate_from_regex(sub_pattern))
        else:
            result.append(char)
        i += 1

    return ''.join(result)

def makebam(ALIGNED_FILE):
    #kinetics = []
    context_list = []
    # ml_list = []
    ml_cnt = 0
    MOTIF = "CA[GCT]C[AG][ACG]C"
    NUM_READS = 100
    bamfile = pysam.AlignmentFile(ALIGNED_FILE,'rb',check_sq=False)
    with pysam.AlignmentFile("manual.bam","wb",template=bamfile,reference_names=None) as manualbam:
        for i in range(NUM_READS):
            seq = generate_from_regex(MOTIF)*100
            a = pysam.AlignedSegment()
            a.query_name = f"read/{i}/ccs"
            a.query_sequence=seq
            a.flag = 4
            a.mapping_quality = 255
            # a.cigar = ((7,len(seq)))
            a.template_length=0
            a.query_qualities = pysam.qualitystring_to_array("~"*len(seq))
            fi_array = np.random.randint(1,100,size=len(seq))
            fi_str = np.array2string(fi_array,separator=",").replace("[","").replace("]","").replace(" ","")
            ri_array = np.random.randint(1,100,size=len(seq))
            ri_str = np.array2string(ri_array,separator=",").replace("[","").replace("]","").replace(" ","")
            fp_array = np.random.randint(1,100,size=len(seq))
            fp_str = np.array2string(fp_array,separator=",").replace("[","").replace("]","").replace(" ","")
            rp_array = np.random.randint(1,100,size=len(seq))
            rp_str = np.array2string(rp_array,separator=",").replace("[","").replace("]","").replace(" ","")
            a.set_tag("RG", "14ecb5f8",value_type="Z")
            a.set_tag("ac", arr.array("i",[60,0,60,0]))
            a.set_tag('bc',arr.array("H",[5,5]))
            a.set_tag("bq",60,value_type="i")
            a.set_tag("ec",60.0,value_type="f")
            a.set_tag("fi",arr.array("B",list(fi_array)))
            a.set_tag("fn",30,value_type="i")
            a.set_tag("fp",arr.array("B",list(fp_array)))
            a.set_tag("ma",0,value_type="i")
            a.set_tag("np",60,value_type="i")
            a.set_tag("ri",arr.array("B",list(ri_array)))
            a.set_tag("rn",30,value_type="i")
            a.set_tag("rp",arr.array("B",list(rp_array)))
            a.set_tag("rq",1,value_type="f")
            a.set_tag("sn",arr.array("f",[8.00000,12.00000,2.00000,4.00000]))
            a.set_tag("zm",i,value_type="i")
            manualbam.write(a)
    bamfile.close()
    return

        
if __name__ == "__main__":
    makebam(sys.argv[1])
    
