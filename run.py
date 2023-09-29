import subsample as ss
ALIGNED_FILE = "/work/gg57/g57016/data/m64121_220824_041525.subreads.actc.forward.bam"
NEW_FILE = "/work/gg57/g57016/data/m64121_220824_041525.subreads.actc.subsample"

ss.subsample(ALIGNED_FILE,NEW_FILE)
