#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM --omp thread=56
#PJM -L elapse=24:00:00
#PJM -L rscgrp=regular

/work/gg57/g57016/tools/miniconda3/bin/ccs --min-passes 10 /work/00/gg57/share/dataset/ce_hifi/reads_220824/m64121_220824_041525.subreads.bam /work/gg57/g57016/data/m64121_220824_041525.ccs.bam
