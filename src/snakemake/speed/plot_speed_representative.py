import sys

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

def read_file(results, files):
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                mean = round(float(line.split('\t')[2]),2)
                stdev = round(float(line.split('\t')[3]),2)
                results.append((mean,stdev))
    return results

minimiser = read_file([], ["0_minimiser_hash_20_"+str(w)+"_speed.out" for w in [i for i in range(24,44,4)]])
minimiser_setw = read_file([], ["0_minimiser_hash_"+str(k)+"_40_speed.out" for k in [i for i in range(16,36,4)]])
# modmer
modmer = read_file([], ["0_modmer_hash_20_"+str(w)+"_speed.out" for w in [3,5,7,9,11]])
modmer_setw = read_file([], ["0_modmer_hash_"+str(k)+"_7_speed.out" for k in [i for i in range(16,36,4)]])
# syncmer
opensyncmer = read_file([], ["syncmer_hash_20_"+str(w)+"_0_0_speed.out" for w in [18,16,14,12,10]])
opensyncmer_setw = read_file([], ["syncmer_hash_"+str(k)+"_10_0_0_speed.out" for k in [i for i in range(22,12,-2)]])
closedsyncmer = read_file([], ["syncmer_hash_20_"+str(w)+"_0_6_speed.out" for w in [15,11,7,3,1]])
closedsyncmer_setw = read_file([], ["syncmer_hash_"+str(k)+"_3_0_6_speed.out" for k in [i for i in range(28,8,-4)]])


# Plot comparison between k-mers
k_size = [i for i in range(5)]
pos = [x+0.25 for x in range(len(k_size))]

fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#890015","#5cffca","#a13ff0","#ff9ba0"]
colors_error = ["#01d63a","#00e7e0","#fefea1","#748beb"]
plt.xlabel("w,m or s")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minimiser], color = colors[0], label='(w,20)-minimizer', linewidth=3.0)
plt.plot(pos, [x[0] for x in modmer], color = colors[1], label='(20,m)-modmer', linewidth=3.0)
plt.plot(pos, [x[0] for x in opensyncmer], color = colors[2], label='(20,s,[0],1)-syncmer', linewidth=3.0)
plt.plot(pos, [x[0] for x in closedsyncmer], color = colors[3], label='(20,s,[0,6],1)-syncmer',linewidth=3.0)

#plt.fill_between(pos, [x[0]-x[1] for x in minimiser], [x[0]+x[1] for x in minimiser], color = colors_error[0], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in modmer], [x[0]+x[1] for x in modmer], color = colors_error[1], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in opensyncmer], [x[0]+x[1] for x in opensyncmer], color = colors_error[2], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in closedsyncmer], [x[0]+x[1] for x in closedsyncmer], color = colors_error[3], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75),title="Methods")
plt.savefig("../results/Speed_representative.png", bbox_inches ='tight')


# Plot comparison between k-mers
k_size = [i for i in range(16,36,4)]
pos = [x+0.25 for x in range(len(k_size))]

fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#890015","#5cffca","#a13ff0","#ff9ba0"]
colors_error = ["#01d63a","#00e7e0","#fefea1","#748beb"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minimiser_setw], color = colors[0], label='(40,k)-minimizer',linewidth=3.0)
plt.plot(pos, [x[0] for x in modmer_setw], color = colors[1], label='(k,7)-modmer',linewidth=3.0)
plt.plot(pos, [x[0] for x in opensyncmer_setw], color = colors[2], label='(k,10,[0],1)-syncmer',linewidth=3.0)
plt.plot(pos, [x[0] for x in closedsyncmer_setw], color = colors[3], label='(k,3,[0,6],1)-syncmer',linewidth=3.0)

#plt.fill_between(pos, [x[0]-x[1] for x in minimiser_setw], [x[0]+x[1] for x in minimiser_setw], color = colors_error[0], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in modmer_setw], [x[0]+x[1] for x in modmer_setw], color = colors_error[1], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in opensyncmer_setw], [x[0]+x[1] for x in opensyncmer_setw], color = colors_error[2], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in closedsyncmer_setw], [x[0]+x[1] for x in closedsyncmer_setw], color = colors_error[3], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75),title="Methods")
plt.savefig("../results/Speed_representative_setw.png", bbox_inches ='tight')
