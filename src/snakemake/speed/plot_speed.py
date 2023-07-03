import sys

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

#k_size = [16,20,24,28,32]
#pos = [x+0.25 for x in range(len(k_size))]
#strobe_range = [int(k/2) for k in k_size]
k_size = [16,18,20,22,24,26,28,30,32]
pos = [x+0.25 for x in range(len(k_size))]
pos_order3 = [1.25,4.25,7.25]
k_order3 = [6,8,10]
k_size_order3 = [i*3 for i in k_order3]
strobe_range = [k for k in range(8,17)]

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
#shapes4 = ["36607","233469","933855","4192891","14548847","62257151","234879855","805287931","3169577727"]
shapes4=['777695', '2621175', '16252901', '50196477', '251620351', '905838335', '4286578095', '13958643693', '66035113981']
gapped4_kmers = read_file([], [shapes4[i] + "_kmer_hash_"+str(k_size[i]+4)+"_speed.out" for i in range(len(k_size))])
#shapes8 = ["51755","246365","975475","3669089","13954519","66560815","241004285","1004529051","3856068575"]
shapes8 = ['14021527', '45607667', '180082591', '1068161519', '3522001919', '13957854679', '64423783901', '205814423455', '1094946651927']
gapped8_kmers = read_file([], [shapes8[i] + "_kmer_hash_"+str(k_size[i]+8)+"_speed.out" for i in range(len(k_size))])

kmers_order3 = read_file([], ["0_kmer_hash_"+str(k)+"_speed.out" for k in k_size_order3])
#shapes4_order3 = ["233469","14548847","805287931"]
shapes4_order3 = ['2621175', '251620351', '13958643693']
gapped4_order3 = read_file([], [shapes4_order3[i] + "_kmer_hash_"+str(k_size_order3[i])+"_speed.out" for i in range(len(k_order3))])
#shapes8_order3 = ["246365","13954519","1004529051"]
shapes8_order3 = ['45607667', '3522001919', '205814423455']
gapped8_order3 = read_file([], [shapes8_order3[i] + "_kmer_hash_"+str(k_size_order3[i])+"_speed.out" for i in range(len(k_order3))])

minstrobemers2 = read_file([], ["minstrobemers_"+str(k)+"_2_"+str(1)+"_"+str(4+k)+"_speed.out" for k in strobe_range])
minstrobemers3 = read_file([], ["minstrobemers_"+str(k)+"_3_"+str(1)+"_"+str(4+k)+"_speed.out" for k in k_order3])
hybridstrobemers2 = read_file([],["hybridstrobemers_"+str(k)+"_2_"+str(1)+"_"+str(4+k)+"_speed.out" for k in strobe_range])
hybridstrobemers3 = read_file([],["hybridstrobemers_"+str(k)+"_3_"+str(1)+"_"+str(4+k)+"_speed.out" for k in k_order3])
randstrobemers2 = read_file([], ["randstrobemers_"+str(k)+"_2_"+str(1)+"_"+str(4+k)+"_speed.out" for k in strobe_range])
randstrobemers3 = read_file([],["randstrobemers_"+str(k)+"_3_"+str(1)+"_"+str(4+k)+"_speed.out" for k in k_order3])
minstrobemers28 = read_file([], ["minstrobemers_"+str(k)+"_2_"+str(1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
minstrobemers38 = read_file([], ["minstrobemers_"+str(k)+"_3_"+str(1)+"_"+str(8+k)+"_speed.out" for k in k_order3])
hybridstrobemers28 = read_file([],["hybridstrobemers_"+str(k)+"_2_"+str(1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
hybridstrobemers38 = read_file([],["hybridstrobemers_"+str(k)+"_3_"+str(1)+"_"+str(8+k)+"_speed.out" for k in k_order3])
randstrobemers28 = read_file([], ["randstrobemers_"+str(k)+"_2_"+str(1)+"_"+str(8+k)+"_speed.out" for k in strobe_range])
randstrobemers38 = read_file([],["randstrobemers_"+str(k)+"_3_"+str(1)+"_"+str(8+k)+"_speed.out" for k in k_order3])
original_randstrobemers2 = read_file([], ["Original_randstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(5+k)+"_speed.out" for k in strobe_range])
original_randstrobemers38 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(9+k)+"_speed.out" for k in k_order3])

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
plt.savefig("../results/Speed_kmers.png")

# Plot comparison between strobemers 4 gaps
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1"]
colors_error = ["#748beb","#e47585","#b261c2"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minstrobemers2], color = colors[0], label='minstrobemers')
plt.plot(pos, [x[0] for x in hybridstrobemers2], color = colors[1], label='hybridstrobemers')
plt.plot(pos, [x[0] for x in randstrobemers2], color = colors[2], label='randstrobemers')

plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers2], [x[0]+x[1] for x in minstrobemers2], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers2], [x[0]+x[1] for x in hybridstrobemers2], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers2], [x[0]+x[1] for x in randstrobemers2], color = colors_error[2], alpha=0.7)

#plt.legend(bbox_to_anchor=(1.25, 0.75), title="Methods")
plt.savefig("../results/Speed_strobemers4.png")

# Plot comparison between strobemers 8 gaps
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1"]
colors_error = ["#748beb","#e47585","#b261c2"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minstrobemers28], color = colors[0], label='minstrobemers')
plt.plot(pos, [x[0] for x in hybridstrobemers28], color = colors[1], label='hybridstrobemers')
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors[2], label='randstrobemers')

plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers28], [x[0]+x[1] for x in minstrobemers28], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers28], [x[0]+x[1] for x in hybridstrobemers28], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers28], [x[0]+x[1] for x in randstrobemers28], color = colors_error[2], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_strobemers8.png", bbox_inches='tight')

# Plot comparison between strobemers 4 gaps order 3
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1"]
colors_error = ["#748beb","#e47585","#b261c2"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos_order3, [x[0] for x in minstrobemers3], color = colors[0], label='minstrobemers')
plt.plot(pos_order3, [x[0] for x in hybridstrobemers3], color = colors[1], label='hybridstrobemers')
plt.plot(pos_order3, [x[0] for x in randstrobemers3], color = colors[2], label='randstrobemers')

plt.fill_between(pos_order3, [x[0]-x[1] for x in minstrobemers3], [x[0]+x[1] for x in minstrobemers3], color = colors_error[0], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in hybridstrobemers3], [x[0]+x[1] for x in hybridstrobemers3], color = colors_error[1], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers3], [x[0]+x[1] for x in randstrobemers3], color = colors_error[2], alpha=0.7)

#plt.legend(bbox_to_anchor=(1.25, 0.75), title="Methods")
plt.savefig("../results/Speed_strobemers4_order3.png")

# Plot comparison between strobemers 8 gaps order 3
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1"]
colors_error = ["#748beb","#e47585","#b261c2"]
plt.xlabel("k")
plt.xticks(pos_order3, k_order3)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos_order3, [x[0] for x in minstrobemers38], color = colors[0], label='minstrobemers')
plt.plot(pos_order3, [x[0] for x in hybridstrobemers38], color = colors[1], label='hybridstrobemers')
plt.plot(pos_order3, [x[0] for x in randstrobemers38], color = colors[2], label='randstrobemers')

plt.fill_between(pos_order3, [x[0]-x[1] for x in minstrobemers38], [x[0]+x[1] for x in minstrobemers38], color = colors_error[0], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in hybridstrobemers38], [x[0]+x[1] for x in hybridstrobemers38], color = colors_error[1], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers38], [x[0]+x[1] for x in randstrobemers38], color = colors_error[2], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_strobemers8_order3.png", bbox_inches='tight')

# Plot comparison between all
fig = plt.figure()
X = np.arange(len(k_size))

colors = ["#004c6d","#009dbe","#00f6ff","#fdcc8a","#fc8d59","#d7301f"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in minstrobemers2], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in hybridstrobemers2], color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in randstrobemers2], color = colors[5], label='randstrobemers',linewidth=3.0)
#plt.plot(pos, [x[0] for x in original_randstrobemers2], color = colors[5], label='randstrobemers', linewidth=3.0)
#plt.plot(pos_order3, [x[0] for x in original_randstrobemers38], color = colors[1], label='8 ori')

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_all.png",bbox_inches='tight')

# Plot comparison between all with 8
fig = plt.figure()
X = np.arange(len(k_size))

plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in minstrobemers28], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in hybridstrobemers28], color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors[5], label='randstrobemers',linewidth=3.0)

#plt.plot(pos, [x[0] for x in original_randstrobemers28], color = colors[5], label='ranstrobemers', linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_all8.png",bbox_inches='tight')

# Plot comparison between all order 3
fig = plt.figure()
X = np.arange(len(k_size_order3))

plt.xlabel("k")
plt.xticks(pos_order3, k_size_order3)
plt.ylabel("Speed in microseconds") # in microseconds

original_randstrobemers3 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(4+k)+"_speed.out" for k in k_order3])

plt.plot(pos_order3, [x[0] for x in kmers_order3], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped4_order3], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped8_order3], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in minstrobemers3], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in hybridstrobemers3], color = colors[4], label='hybridstrobemers',linewidth=3.0)
#plt.plot(pos_order3, [x[0] for x in randstrobemers3], color = colors[5], label='randstrobemers',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in original_randstrobemers3], color = colors[5], label='randstrobemers',linewidth=3.0)
#plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_all_order3.png",bbox_inches='tight')

# Plot comparison between all with 8 order 3
fig = plt.figure()
X = np.arange(len(k_size_order3))

plt.xlabel("k")
plt.xticks(pos_order3, k_size_order3)
plt.ylabel("Speed in microseconds") # in microseconds

original_randstrobemers38 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(8+k)+"_speed.out" for k in k_order3])

plt.plot(pos_order3, [x[0] for x in kmers_order3], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped4_order3], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped8_order3], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in minstrobemers38], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in hybridstrobemers38], color = colors[4], label='hybridstrobemers',linewidth=3.0)
#plt.plot(pos_order3, [x[0] for x in randstrobemers38], color = colors[5], label='randstrobemers',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in original_randstrobemers38], color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_all8_order3.png",bbox_inches='tight')

# Plot comparison between strobemers all gaps
fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1","#00ba32","#00d6e7","#fad100"]
colors_error = ["#748beb","#e47585","#b261c2","#01d63a","#00e7e0","#fefea1"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in minstrobemers2], color = colors[0], label='4 minstrobemers')
plt.plot(pos_order3, [x[0] for x in minstrobemers3], color = colors[0], label='4 randstrobemers3',linestyle="dashed")
plt.plot(pos, [x[0] for x in hybridstrobemers2], color = colors[1], label='4 hybridstrobemers')
plt.plot(pos_order3, [x[0] for x in hybridstrobemers3], color = colors[1], label='4 randstrobemers3',linestyle="dashed")
plt.plot(pos, [x[0] for x in randstrobemers2], color = colors[2], label='4 randstrobemers')
plt.plot(pos_order3, [x[0] for x in randstrobemers3], color = colors[2], label='4 randstrobemers3',linestyle="dashed")
plt.plot(pos, [x[0] for x in minstrobemers28], color = colors[3], label='8 minstrobemers')
plt.plot(pos_order3, [x[0] for x in minstrobemers38], color = colors[3], label='4 randstrobemers3',linestyle="dashed")
plt.plot(pos, [x[0] for x in hybridstrobemers28], color = colors[4], label='8 hybridstrobemers')
plt.plot(pos_order3, [x[0] for x in hybridstrobemers38], color = colors[4], label='4 randstrobemers3',linestyle="dashed")
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors[5], label='8 randstrobemers')
plt.plot(pos_order3, [x[0] for x in randstrobemers38], color = colors[5], label='8 randstrobemers3',linestyle="dashed")

plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers2], [x[0]+x[1] for x in minstrobemers2], color = colors_error[0], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in minstrobemers3], [x[0]+x[1] for x in minstrobemers3], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers2], [x[0]+x[1] for x in hybridstrobemers2], color = colors_error[1], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in hybridstrobemers3], [x[0]+x[1] for x in hybridstrobemers3], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers2], [x[0]+x[1] for x in randstrobemers2], color = colors_error[2], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers3], [x[0]+x[1] for x in randstrobemers3], color = colors_error[2], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers28], [x[0]+x[1] for x in minstrobemers28], color = colors_error[3], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in minstrobemers38], [x[0]+x[1] for x in minstrobemers38], color = colors_error[3], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers28], [x[0]+x[1] for x in hybridstrobemers28], color = colors_error[4], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in hybridstrobemers38], [x[0]+x[1] for x in hybridstrobemers38], color = colors_error[4], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers28], [x[0]+x[1] for x in randstrobemers28], color = colors_error[5], alpha=0.7)
plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers38], [x[0]+x[1] for x in randstrobemers38], color = colors_error[5], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_strobemers.png", bbox_inches='tight')


# Plot comparison between strobemers all gaps
original_minstrobemers2 = read_file([], ["Original_minstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(5+k)+"_speed.out" for k in strobe_range])
original_hybridstrobemers2 = read_file([],["Original_hybridstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(5+k)+"_speed.out" for k in strobe_range])
original_randstrobemers2 = read_file([], ["Original_randstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(5+k)+"_speed.out" for k in strobe_range])
original_randstrobemers3 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(5+k)+"_speed.out" for k in k_order3])
original_minstrobemers28 = read_file([], ["Original_minstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(9+k)+"_speed.out" for k in strobe_range])
original_hybridstrobemers28 = read_file([],["Original_hybridstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(9+k)+"_speed.out" for k in strobe_range])
original_randstrobemers28 = read_file([], ["Original_randstrobemers_"+str(k)+"_2_"+str(k+1)+"_"+str(9+k)+"_speed.out" for k in strobe_range])
original_randstrobemers38 = read_file([],["Original_randstrobemers_"+str(k)+"_3_"+str(k+1)+"_"+str(9+k)+"_speed.out" for k in k_order3])

fig = plt.figure()
X = np.arange(len(k_size))
colors = ["#697ed5","#c76674","#9350a1","#00ba32","#00d6e7","#fad100"]
colors_error = ["#748beb","#e47585","#b261c2","#01d63a","#00e7e0","#fefea1"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in original_minstrobemers2], color = colors[0], label='4 minstrobemers')
plt.plot(pos, [x[0] for x in original_hybridstrobemers2], color = colors[1], label='4 hybridstrobemers')
plt.plot(pos, [x[0] for x in original_randstrobemers2], color = colors[2], label='4 randstrobemers')
plt.plot(pos, [x[0] for x in original_minstrobemers28], color = colors[3], label='8 minstrobemers')
plt.plot(pos, [x[0] for x in original_hybridstrobemers28], color = colors[4], label='8 hybridstrobemers')
plt.plot(pos, [x[0] for x in original_randstrobemers28], color = colors[5], label='8 randstrobemers')

plt.fill_between(pos, [x[0]-x[1] for x in original_minstrobemers2], [x[0]+x[1] for x in original_minstrobemers2], color = colors_error[0], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in original_hybridstrobemers2], [x[0]+x[1] for x in original_hybridstrobemers2], color = colors_error[1], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in original_randstrobemers2], [x[0]+x[1] for x in original_randstrobemers2], color = colors_error[2], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in original_minstrobemers28], [x[0]+x[1] for x in original_minstrobemers28], color = colors_error[3], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in original_hybridstrobemers28], [x[0]+x[1] for x in original_hybridstrobemers28], color = colors_error[4], alpha=0.7)
plt.fill_between(pos, [x[0]-x[1] for x in original_randstrobemers28], [x[0]+x[1] for x in original_randstrobemers28], color = colors_error[5], alpha=0.7)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Speed_strobemers_original_all.png", bbox_inches='tight')


fig = plt.figure()
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds
colors_ori = ["#bae4bc","#7bccc4","#43a2ca","#0868ac"]

plt.plot(pos, [x[0] for x in minstrobemers2], color = colors_ori[0], label='4')
plt.plot(pos, [x[0] for x in minstrobemers28], color = colors_ori[1], label='8')
plt.plot(pos, [x[0] for x in original_minstrobemers2], color = colors_ori[2], label='4 ori')
plt.plot(pos, [x[0] for x in original_minstrobemers28], color = colors_ori[3], label='8 ori')

#plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers2], [x[0]+x[1] for x in minstrobemers2], color = colors_error[2], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in minstrobemers28], [x[0]+x[1] for x in minstrobemers28], color = colors_error[5], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in original_minstrobemers2], [x[0]+x[1] for x in original_minstrobemers2], color = colors_error[0], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in original_minstrobemers28], [x[0]+x[1] for x in original_minstrobemers28], color = colors_error[1], alpha=0.7)

plt.legend(title="Methods")
plt.savefig("../results/Speed_minstrobemers_original.png", bbox_inches='tight')

fig = plt.figure()
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in hybridstrobemers2], color = colors_ori[0], label='4')
plt.plot(pos, [x[0] for x in hybridstrobemers28], color = colors_ori[1], label='8')
plt.plot(pos, [x[0] for x in original_hybridstrobemers2], color = colors_ori[2], label='4 ori')
plt.plot(pos, [x[0] for x in original_hybridstrobemers28], color = colors_ori[3], label='8 ori')

#plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers2], [x[0]+x[1] for x in hybridstrobemers2], color = colors_error[2], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in hybridstrobemers28], [x[0]+x[1] for x in hybridstrobemers28], color = colors_error[5], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in original_hybridstrobemers2], [x[0]+x[1] for x in original_hybridstrobemers2], color = colors_error[0], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in original_hybridstrobemers28], [x[0]+x[1] for x in original_hybridstrobemers28], color = colors_error[1], alpha=0.7)

plt.legend(title="Methods")
plt.savefig("../results/Speed_hybridstrobemers_original.png", bbox_inches='tight')

fig = plt.figure()
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos, [x[0] for x in randstrobemers2], color = colors_ori[0], label='4')
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors_ori[1], label='8')
plt.plot(pos, [x[0] for x in original_randstrobemers2], color = colors_ori[2], label='4 ori')
plt.plot(pos, [x[0] for x in original_randstrobemers28], color = colors_ori[3], label='8 ori')

#plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers2], [x[0]+x[1] for x in randstrobemers2], color = colors_error[2], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in randstrobemers28], [x[0]+x[1] for x in randstrobemers28], color = colors_error[5], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in original_randstrobemers2], [x[0]+x[1] for x in original_randstrobemers2], color = colors_error[0], alpha=0.7)
#plt.fill_between(pos, [x[0]-x[1] for x in original_randstrobemers28], [x[0]+x[1] for x in original_randstrobemers28], color = colors_error[1], alpha=0.7)

plt.legend(title="Methods")
plt.savefig("../results/Speed_randstrobemers_original.png", bbox_inches='tight')

fig = plt.figure()
plt.xlabel("k")
plt.xticks(pos_order3, k_order3)
plt.ylabel("Speed in microseconds") # in microseconds

plt.plot(pos_order3, [x[0] for x in randstrobemers3], color = colors_ori[0], label='4')
plt.plot(pos_order3, [x[0] for x in randstrobemers38], color = colors_ori[1], label='8')
plt.plot(pos_order3, [x[0] for x in original_randstrobemers3], color = colors_ori[2], label='4 ori')
plt.plot(pos_order3, [x[0] for x in original_randstrobemers38], color = colors_ori[3], label='8 ori')

#plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers3], [x[0]+x[1] for x in randstrobemers3], color = colors_error[2], alpha=0.7)
#plt.fill_between(pos_order3, [x[0]-x[1] for x in randstrobemers38], [x[0]+x[1] for x in randstrobemers38], color = colors_error[5], alpha=0.7)
#plt.fill_between(pos_order3, [x[0]-x[1] for x in original_randstrobemers3], [x[0]+x[1] for x in original_randstrobemers3], color = colors_error[0], alpha=0.7)
#plt.fill_between(pos_order3, [x[0]-x[1] for x in original_randstrobemers38], [x[0]+x[1] for x in original_randstrobemers38], color = colors_error[1], alpha=0.7)

plt.legend(title="Methods")
plt.savefig("../results/Speed_randstrobemers_original_order3.png", bbox_inches='tight')
