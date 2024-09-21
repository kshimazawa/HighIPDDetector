import numpy as np
import matplotlib.pyplot as plt

SAMPLE1 = input()
SAMPLE2 = input()
SAVE_TO = input()
s1 = np.loadtxt(SAMPLE1,delimiter=",")
s2 = np.loadtxt(SAMPLE2,delimiter=",")
edge_s1 = 0
edge_s2 = 0
for i in reversed(range(s1.shape[0])):
    if s1[i,1] != 0:
        edge_s1 = i
        break
for i in reversed(range(s2.shape[0])):
    if s2[i,1] != 0:
        edge_s2 = i
        break
edge = max(edge_s1,edge_s2)
s1 = s1[:edge+1]
s2 = s2[:edge+1]
max_coverage_1 = int(np.max(s1[:,1]))
max_coverage_2 = int(np.max(s2[:,1]))
#fig = plt.figure()
#ax = fig.add_subplot(221)
#ax.hist(s1[:,1],bins=np.arange(100))
#ax = fig.add_subplot(222)
#ax.hist(s2[:,1],bins=np.arange(100))
#var_1 = np.zeros(max_coverage_1 + 1)
#var_2 = np.zeros(max_coverage_2 + 1)
#for i in range(max_coverage_1+1):
#    var_1[i] =  np.var(s1[s1[:,1]==i,0]/s1[s1[:,1]==i,1])
    #var_1[i] = np.mean(s1[s1[:,1]==i,2] - (s1[s1[:,1]==i,0]/s1[s1[:,1]==i,1])**2)
#ax = fig.add_subplot(223)
#ax.plot(np.arange(max_coverage_1+1),var_1)
#for i in range(max_coverage_2+1):
#    var_2[i] =  np.var(s2[s2[:,1]==i,0]/s2[s2[:,1]==i,1])
#    var_2[i] = np.mean(s2[s2[:,1]==i,2] - (s2[s2[:,1]==i,0]/s2[s2[:,1]==i,1])**2)
#ax = fig.add_subplot(224)
#ax.plot(np.arange(max_coverage_2+1),var_2)
s1_mean = np.zeros(edge + 1)
s2_mean = np.zeros(edge + 1)
for i in range(edge+1):
    if s1[i,1] != 0:
        s1_mean[i] = s1[i,0]/s1[i,1]
for i in range(edge+1):
    if s2[i,1] != 0:
        s2_mean[i] = s2[i,0]/s2[i,1]
ratio = np.zeros(edge+1)
for i in range(ratio.shape[0]):
    if s2_mean[i] != 0:
        ratio[i] = s1_mean[i]/s2_mean[i]
np.savetxt("{}ratio.txt".format(SAVE_TO),ratio,fmt="%0.5f")
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(221)
ax.scatter(np.arange(edge+1),ratio,s=0.1)
low_coverage = np.argwhere(s1[:,1]<=5)
ax = fig.add_subplot(223)
ax.hist(low_coverage,bins=100)
low_coverage = np.argwhere(s2[:,1]<=5)
ax = fig.add_subplot(224)
ax.hist(low_coverage,bins=100)
#ax.plot(np.arange(edge_s1+1),s1_mean)
#ax.plot(np.arange(edge_s2+1),s2_mean)
fig.savefig("{}comparison.png".format(SAVE_TO))
