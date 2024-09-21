
import re
FILE ='/home/yoshihiko_s/work/ccs_ipd/static/ce11rel606/ce11rel606.fa'
#MOTIF = '[gG][gG][ACGTacgt][gG][gG][ACGTacgt][gG][gG][ACGTacgt][gG][gG][ACGTacgt]'
#MOTIF = '[AGCTagct][Cc][Cc][ACGTacgt][Cc][cC][ACGTacgt][Cc][Cc][ACTGacgt][Cc][Cc]'
#MOTIF = '[Aa][Tt][Cc][Aa][Gg][Cc][Tt][Gg]'
#MOTIF = '[Tt][Gg][Aa][Cc][Gg][Tt][Cc][Aa]'
#MOTIF = '[Gg][Aa][Tt][Cc]'
#MOTIF = '[Aa][Tt][Aa][Tt][Aa][Tt][Aa][Tt][Aa][Tt][Aa][Tt][Aa][Tt][Aa][Tt][Aa][Tt][Aa][Tt]'
#MOTIF = '[Aa][Gg]'
#MOTIF = '[AT][TA][AT][TA][AG][TA][AT][GTA]GTATA[AT][GCTA][TACG][TA][TA]'
#MOTIF = '[TA][AT][AGC][TAC]TG[TC]TCA[GT][ACT]C[AG][TA][TCA]'
MOTIF = '[Gg][Gg][Aa][Gg]'
def find_motif(fname = FILE,motif = MOTIF):
    pos = {}
    is_first = True
    with open(fname,"r") as fr:
        seq = ""
        key = ""
        for line in fr:
            if (">" in line) & is_first:
                key = line.replace(">","")
                pos[key] = []
                is_first = False
            elif ">" in line:
                pos[key] = [m.span() for m in re.finditer(MOTIF, seq)]
                key = line.replace(">","")
                pos[key] = []
                seq = ""
            else:
                seq = seq + line.strip("\n")
        pos[key] = [m.span() for m in re.finditer(MOTIF, seq)]
    for k in pos:
        print(k)
        for i in range(len(pos[k])):
            print(pos[k][i][0])
    return

if __name__ == '__main__':
    find_motif(FILE,MOTIF)
