import sys

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

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
kmers = read_file([], ["0_minimiser_hash_"+str(k)+"_"+str(k)+"_counts.out" for k in k_size])
#shapes4 = ["36607","233469","933855","4192891","14548847","62257151","234879855","805287931","3169577727"]
shapes4=['777695', '2621175', '16252901', '50196477', '251620351', '905838335', '4286578095', '13958643693', '66035113981']
gapped4_kmers = read_file([], [shapes4[i] + "_minimiser_hash_"+str(k_size[i]+4)+"_"+str(k_size[i]+4)+"_counts.out" for i in range(len(k_size))])
#shapes8 = ["51755","246365","975475","3669089","13954519","66560815","241004285","1004529051","3856068575"]
shapes8 = ['14021527', '45607667', '180082591', '1068161519', '3522001919', '13957854679', '64423783901', '205814423455', '1094946651927']
gapped8_kmers = read_file([], [shapes8[i] + "_minimiser_hash_"+str(k_size[i]+8)+"_"+str(k_size[i]+8)+"_counts.out" for i in range(len(k_size))])

kmers_order3 = read_file([], ["0_minimiser_hash_"+str(k)+"_"+str(k)+"_counts.out" for k in k_size_order3])
#shapes4_order3 = ["233469","14548847","805287931"]
shapes4_order3 = ['2621175', '251620351', '13958643693']
gapped4_order3 = read_file([], [shapes4_order3[i] + "_minimiser_hash_"+str(k_size_order3[i]+4)+"_"+str(k_size_order3[i]+4)+"_counts.out" for i in range(len(k_order3))])
#shapes8_order3 = ["246365","13954519","1004529051"]
shapes8_order3 = ['45607667', '3522001919', '205814423455']
gapped8_order3 = read_file([], [shapes8_order3[i] + "_minimiser_hash_"+str(k_size_order3[i]+8)+"_"+str(k_size_order3[i]+8)+"_counts.out" for i in range(len(k_order3))])

minstrobemers2 = read_file([], ["minstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(3+k)+"_counts.out" for k in strobe_range])
minstrobemers3 = read_file([], ["minstrobemers_"+str(k)+"_3_"+str(0)+"_"+str(3+k)+"_counts.out" for k in k_order3])
hybridstrobemers2 = read_file([],["hybridstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(4+k)+"_counts.out" for k in strobe_range])
hybridstrobemers3 = read_file([],["hybridstrobemers_"+str(k)+"_3_"+str(0)+"_"+str(4+k)+"_counts.out" for k in k_order3])
randstrobemers2 = read_file([], ["randstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(3+k)+"_counts.out" for k in strobe_range])
randstrobemers3 = read_file([],["randstrobemers_"+str(k)+"_3_"+str(0)+"_"+str(3+k)+"_counts.out" for k in k_order3])
minstrobemers28 = read_file([], ["minstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_counts.out" for k in strobe_range])
minstrobemers38 = read_file([], ["minstrobemers_"+str(k)+"_3_"+str(0)+"_"+str(7+k)+"_counts.out" for k in k_order3])
hybridstrobemers28 = read_file([],["hybridstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_counts.out" for k in strobe_range])
hybridstrobemers38 = read_file([],["hybridstrobemers_"+str(k)+"_3_"+str(0)+"_"+str(7+k)+"_counts.out" for k in k_order3])
randstrobemers28 = read_file([], ["randstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_counts.out" for k in strobe_range])
randstrobemers38 = read_file([],["randstrobemers_"+str(k)+"_3_"+str(0)+"_"+str(7+k)+"_counts.out" for k in k_order3])

# Plot comparison between all
fig = plt.figure()
X = np.arange(len(k_size))

colors = ["#004c6d","#009dbe","#00f6ff","#fdcc8a","#fc8d59","#d7301f"]
plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("# of Submers")
plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in minstrobemers2], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in hybridstrobemers2], color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in randstrobemers2], color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Count_all.png",bbox_inches='tight')

# Plot comparison between all
fig = plt.figure()
X = np.arange(len(k_size))

plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("% of Submers")
divs = [0.5*(4**k) for k in k_size]
divs4 = [0.5*(4**(k-4)) for k in k_size]
divs8 = [0.5*(4**(k-8)) for k in k_size]
plt.plot(pos, [kmers[i][0]/divs[i] for i in range(len(k_size))], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, [gapped4_kmers[i][0]/divs4[i] for i in range(len(k_size))], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, [gapped8_kmers[i][0]/divs8[i] for i in range(len(k_size))], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, [minstrobemers2[i][0]/divs[i] for i in range(len(k_size))], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, [hybridstrobemers2[i][0]/divs[i] for i in range(len(k_size))], color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, [randstrobemers2[i][0]/divs[i] for i in range(len(k_size))], color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Count_all_percentage.png",bbox_inches='tight')

# Plot comparison between all with 8
fig = plt.figure()
X = np.arange(len(k_size))

plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("# of Submers")

plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, [x[0] for x in minstrobemers28], color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in hybridstrobemers28], color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, [x[0] for x in randstrobemers28], color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Count_all8.png",bbox_inches='tight')

# Plot comparison between all
fig = plt.figure()
X = np.arange(len(k_order3))

plt.xlabel("k")
plt.xticks(pos_order3, k_order3)
plt.ylabel("# of Submers")
plt.plot(pos_order3, [x[0] for x in kmers_order3], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped4_order3], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped8_order3], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in minstrobemers3], color = colors[3], label='4 minstrobemers3',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in hybridstrobemers3], color = colors[4], label='4 hybridstrobemers3',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in randstrobemers3], color = colors[5], label='4 randstrobemers3',linewidth=3.0)

#plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Count_all3.png",bbox_inches='tight')

# Plot comparison between all
fig = plt.figure()
X = np.arange(len(k_order3))

plt.xlabel("k")
plt.xticks(pos_order3, k_order3)
plt.ylabel("# of Submers")
plt.plot(pos_order3, [x[0] for x in kmers_order3], color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped4_order3], color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in gapped8_order3], color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in minstrobemers38], color = colors[3], label='4 minstrobemers3',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in hybridstrobemers38], color = colors[4], label='4 hybridstrobemers3',linewidth=3.0)
plt.plot(pos_order3, [x[0] for x in randstrobemers38], color = colors[5], label='4 randstrobemers3',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Count_all38.png",bbox_inches='tight')

# Plot Uniqueness
kmers = []
gapped4 = []
gapped8 = []
kmers_order3 = []
gapped4_order3 = []
gapped8_order3 = []
minstrobemers2 = []
minstrobemers3 = []
hybridstrobemers2 = []
hybridstrobemers3 = []
randstrobemers2 = []
randstrobemers3 = []
minstrobemers28 = []
minstrobemers38 = []
hybridstrobemers28 = []
hybridstrobemers38 = []
randstrobemers28 = []
randstrobemers38 = []

it = 0
with open("../results/Unique.out", 'r') as f:
    for line in f:
        number = float(line.split()[1])
        if (it < 27):
            mod = it % 3
        if (it < 9):
            kmers.append(number)
        elif (it < 18):
            gapped4.append(number)
        elif (it < 27):
            gapped8.append(number)
        elif (it < 36):
            minstrobemers2.append(number)
        elif (it < 39):
            minstrobemers3.append(number)
        elif (it < 48):
            hybridstrobemers2.append(number)
        elif (it < 51):
            hybridstrobemers3.append(number)
        elif (it < 60):
            randstrobemers2.append(number)
        elif (it < 63):
            randstrobemers3.append(number)
        elif (it < 72):
            minstrobemers28.append(number)
        elif (it < 75):
            minstrobemers38.append(number)
        elif (it < 84):
            hybridstrobemers28.append(number)
        elif (it < 87):
            hybridstrobemers38.append(number)
        elif (it < 96):
            randstrobemers28.append(number)
        elif (it < 99):
            randstrobemers38.append(number)
        it += 1

kmers_order3 = [kmers[1],kmers[4],kmers[7]]
gapped4_order3 = [gapped4[1],gapped4[4],gapped4[7]]
gapped8_order3 = [gapped8[1],gapped8[4],gapped8[7]]

# Plot comparison between all
fig = plt.figure()
X = np.arange(len(k_size))

plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("% of unique submers")
plt.plot(pos, kmers, color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, gapped4, color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, gapped8, color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, minstrobemers2, color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, hybridstrobemers2, color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, randstrobemers2, color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Unique_all.png",bbox_inches='tight')

# Plot comparison between all 8
fig = plt.figure()
X = np.arange(len(k_size))

plt.xlabel("k")
plt.xticks(pos, k_size)
plt.ylabel("% of unique submers")
plt.plot(pos, kmers, color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos, gapped4, color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos, gapped8, color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos, minstrobemers28, color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos, hybridstrobemers28, color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos, randstrobemers28, color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Unique_all8.png",bbox_inches='tight')

# Plot comparison between all order 3
fig = plt.figure()
X = np.arange(len(k_order3))

plt.xlabel("k")
plt.xticks(pos_order3, k_order3)
plt.ylabel("% of unique submers")
plt.plot(pos_order3, kmers_order3, color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos_order3, gapped4_order3, color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos_order3, gapped8_order3, color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos_order3, minstrobemers3, color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos_order3, hybridstrobemers3, color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos_order3, randstrobemers3, color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Unique_all_order3.png",bbox_inches='tight')

# Plot comparison between all order 3 8
fig = plt.figure()
X = np.arange(len(k_order3))

plt.xlabel("k")
plt.xticks(pos_order3, k_order3)
plt.ylabel("% of unique submers")
plt.plot(pos_order3, kmers_order3, color = colors[0], label='k-mer', linewidth=3.0)
plt.plot(pos_order3, gapped4_order3, color = colors[1], label='4 k-mer',linewidth=3.0)
plt.plot(pos_order3, gapped8_order3, color = colors[2], label='8 k-mer',linewidth=3.0)
plt.plot(pos_order3, minstrobemers38, color = colors[3], label='minstrobemers',linewidth=3.0)
plt.plot(pos_order3, hybridstrobemers38, color = colors[4], label='hybridstrobemers',linewidth=3.0)
plt.plot(pos_order3, randstrobemers38, color = colors[5], label='randstrobemers',linewidth=3.0)

plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
plt.savefig("../results/Unique_all_order38.png",bbox_inches='tight')
