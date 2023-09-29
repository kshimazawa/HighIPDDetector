#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L elapse=08:00:00
#PJM -L rscgrp=short

/work/00/gg57/g57016/tools/miniconda3/bin/meme /work/gg57/g57016/tools/chrM.fasta -dna -nmotifs 5 -maxsize 1000000 -o meme_out_chrM
