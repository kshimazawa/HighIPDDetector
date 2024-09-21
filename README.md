# HighIPDDetector
HighIPDDetector detects consistent signals of inter pulse duration for reads output by PacBio's next generation sequencing. This code extracts values of IPD, aligns them to make matrices, filters the matrices using heuristic algorithms,and conduct statistical testing to find the consistent signals.
## Requirements
Python (>3.9)
## How to use
Make sure you have a subread file that has the alignment information to its CCS read. If you do not have the alignment, run the following command before running this code.
```
ccs filename.subreads.bam filename.ccs.bam --hifi-kinetics
actc filename.subreads.bam filename.ccs.bam filename.subreads.aligned.bam
```
In addition, the created aligned subread file has to be divided into two files: forward and reverse complement strands to their CCS read. (filename.subreads.aligned.<forward/reverse>.bam)

After splitting the file into two, simply run the following command.
```
python3 dhipd.py <input.txt
```
Where the input.txt has all the parameters/paths as the following:
```
path_to_subreads.aligned.<forward/reverse>.bam
path_to_result_file
100000 #input a maximum read length as the third argument
1000 #input a maximum pass count as the fourth argument
0.01 #input the probability of having a high IPD value as the fifth argument
0.01 #input the statistical significance as the sixth argument
```
