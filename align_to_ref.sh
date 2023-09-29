#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM --omp thread=56
#PJM -L elapse=24:00:00
#PJM -L rscgrp=regular

/work/gg57/g57016/tools/miniconda3/bin/pbmm2 align -m 2G  /work/00/gg57/share/genomes/ce11rel606/ce11rel606.fa /work/gg57/g57016/data/m64121_220824_041525.ccs.bam /work/gg57/g57016/data/m64121_220824_041525.ccs.aligned.bam
