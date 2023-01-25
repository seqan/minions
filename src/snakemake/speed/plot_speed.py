import sys

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

k_size = [16,20,24,28,32]
pos = [x+0.25 for x in range(len(k_size))]

def read_file(results, files):
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                mean = round(float(line.split('\t')[2]),2)
                stdev = round(float(line.split('\t')[3]),2)
                results.append((mean,stdev))
    return results

# Read all files
kmers = read_file([], ["0_kmer_hash_"+str(k)+"_speed.out" for k in k_size])
shapes4 = ["36607","933855","14548847","3758077695","3169577727"]
gapped4_kmers = read_file([], [shapes4[i] + "_kmer_hash_"+str(k_size[i])+"_speed.out" for i in range(len(k_size))])
shapes8 = ["51755","975475","13954519","241004285","241004285"]
gapped8_kmers = read_file([], [shapes8[i] + "_kmer_hash_"+str(k_size[i])+"_speed.out" for i in range(len(k_size))])
minstrobemers = read_file([], ["Original_minstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(4+k)+"_speed.out" for k in strobe_range])
hybridstrobemers = read_file([],["Original_hybridstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(5+k)+"_speed.out" for k in strobe_range])
randstrobemers2 = read_file([], ["Original_randstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(4+k)+"_speed.out" for k in strobe_range])
randstrobemers3 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(4+k)+"_speed.out" for k in [9,12,15]])
minstrobemers8 = read_file([], ["Original_minstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
hybridstrobemers8 = read_file([],["Original_hybridstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
randstrobemers28 = read_file([], ["Original_randstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
randstrobemers38 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in [9,12,15]])

# Plot comparison between k-mers
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#00ba32","#00d6e7","#fad100"]
colors_error = ["#01d63a","#00e7e0","#fefea1"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='0')
plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4')
plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8')


plt.fill_between(pos, [x[0]-x[1] for x in kmers], [x[0]+x[1] for x in kmers], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in gapped4_kmers], [x[0]+x[1] for x in gapped4_kmers], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in gapped8_kmers], [x[0]+x[1] for x in gapped8_kmers], color = colors_error[2], alpha=0.7)

plt.legend(title="Number of gaps")
plt.savefig("Speed_kmers.png")

# Plot comparison between strobemers 4 gaps
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1"]
colors_error = ["#748beb","#e47585","#b261c2"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minstrobemers], color = colors[0], label='minstrobemers')
plt.plot(pos, [x[0] for x in hybridstrobemers], color = colors[1], label='hybridstrobemers')
plt.plot(pos, [x[0] for x in randstrobemers2], color = colors[2], label='randstrobemers')

plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers], [x[0]+x[1] for x in minstrobemers], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers], [x[0]+x[1] for x in hybridstrobemers], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers2], [x[0]+x[1] for x in randstrobemers2], color = colors_error[2], alpha=0.7)

#plt.legend(bbox_to_anchor=(1.25, 0.75), title="Methods")
plt.savefig("Speed_strobemers4.png")

# Plot comparison between strobemers 8 gaps
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1"]
colors_error = ["#748beb","#e47585","#b261c2"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minstrobemers8], color = colors[0], label='minstrobemers')
plt.plot(pos, [x[0] for x in hybridstrobemers8], color = colors[1], label='hybridstrobemers')
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors[2], label='randstrobemers')

plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers8], [x[0]+x[1] for x in minstrobemers8], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers8], [x[0]+x[1] for x in hybridstrobemers8], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers28], [x[0]+x[1] for x in randstrobemers28], color = colors_error[2], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("Speed_strobemers8.png", bbox_inches='tight')

# Plot comparison between all
fig = plt.figure()
X = np.arange(len(k_size))

colors = ["#00ba32","#00d6e7","#fad100","#697ed5","#c76674","#9350a1"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in minstrobemers], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in hybridstrobemers], color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in randstrobemers2], color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("Speed_all.png",bbox_inches='tight')

# Plot comparison between all with 8
fig = plt.figure()
X = np.arange(len(k_size))

colors = ["#00ba32","#00d6e7","#fad100","#697ed5","#c76674","#9350a1"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in minstrobemers8], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in hybridstrobemers8], color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("Speed_all8.png",bbox_inches='tight')


# Plot comparison between strobemers all gaps
k_size = [16,18,20,22,24,26,28,30,32]
pos = [x+0.25 for x in range(len(k_size))]
pos_order3 = [1.25,4.25,7.25]
strobe_range = [k for k in range(8,17)]

minstrobemers = read_file([], ["Original_minstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(4+k)+"_speed.out" for k in strobe_range])
hybridstrobemers = read_file([],["Original_hybridstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(5+k)+"_speed.out" for k in strobe_range])
randstrobemers2 = read_file([], ["Original_randstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(4+k)+"_speed.out" for k in strobe_range])
randstrobemers3 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(4+k)+"_speed.out" for k in [9,12,15]])
minstrobemers8 = read_file([], ["Original_minstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
hybridstrobemers8 = read_file([],["Original_hybridstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
randstrobemers28 = read_file([], ["Original_randstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
randstrobemers38 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in [9,12,15]])

fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1","#00ba32","#00d6e7","#fad100"]
colors_error = ["#748beb","#e47585","#b261c2","#01d63a","#00e7e0","#fefea1"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minstrobemers], color = colors[0], label='4 minstrobemers')
plt.plot(pos, [x[0] for x in hybridstrobemers], color = colors[1], label='4 hybridstrobemers')
plt.plot(pos, [x[0] for x in randstrobemers2], color = colors[2], label='4 randstrobemers')
plt.plot(pos_order3, [x[0] for x in randstrobemers3], color = colors[2], label='4 randstrobemers3',linestyle="dashed")
plt.plot(pos, [x[0] for x in minstrobemers8], color = colors[3], label='8 minstrobemers')
plt.plot(pos, [x[0] for x in hybridstrobemers8], color = colors[4], label='8 hybridstrobemers')
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors[5], label='8 randstrobemers')
plt.plot(pos_order3, [x[0] for x in randstrobemers38], color = colors[5], label='8 randstrobemers3',linestyle="dashed")

plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers], [x[0]+x[1] for x in minstrobemers], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers], [x[0]+x[1] for x in hybridstrobemers], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers2], [x[0]+x[1] for x in randstrobemers2], color = colors_error[2], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers3], [x[0]+x[1] for x in randstrobemers3], color = colors_error[2], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers8], [x[0]+x[1] for x in minstrobemers8], color = colors_error[3], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers8], [x[0]+x[1] for x in hybridstrobemers8], color = colors_error[4], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers28], [x[0]+x[1] for x in randstrobemers28], color = colors_error[5], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers38], [x[0]+x[1] for x in randstrobemers38], color = colors_error[5], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("Speed_strobemers.png", bbox_inches='tight')
