import sys

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

#k_size = [16,20,24,28,32]
#pos = [x+0.25 for x in range(len(k_size))]
#strobe_range = [int(k/2) for k in k_size]
k_size = [16,20,24,28,32]
pos = [x+0.25 for x in range(len(k_size))]
pos_order3 = [1.25,4.25,7.25]
k_order3 = [9,12,15]
k_size_order3 = [i*2 for i in k_order3]
strobe_range = [k for k in range(8,17,2)]

def read_file(results, files):
    cov = 0.0
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                if (line[:7]=="Match C"):
                    cov = round(float(line.split()[2]),2)
                if (line[:7]=="Islands"):
                    mean = round(float(line.split('\t')[2]),2)
                    stdev = round(float(line.split('\t')[3]),2)
                    results.append((mean,stdev,cov))
    return results

# Read all files for an error
for error in [1,2,5,10]:
    kmers = read_file([], ["0_minimiser_hash_"+str(k)+"_"+str(k)+"_match_"+str(error)+".out" for k in range(16,36,4)])
    shapes4 = ["36607","933855","14548847","234879855","3169577727"]
    gapped4_kmers = read_file([], [shapes4[i] + "_minimiser_hash_"+str(k_size[i])+"_"+str(k_size[i])+"_match_"+str(error)+".out" for i in range(len(k_size))])
    shapes8 = ["51755","975475","13954519","241004285","3856068575"]
    gapped8_kmers = read_file([], [shapes8[i] + "_minimiser_hash_"+str(k_size[i])+"_"+str(k_size[i])+"_match_"+str(error)+".out" for i in range(len(k_size))])

    minstrobemers2 = read_file([], ["minstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(3+k)+"_match_"+str(error)+".out" for k in strobe_range])
    hybridstrobemers2 = read_file([],["hybridstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(4+k)+"_match_"+str(error)+".out" for k in strobe_range])
    randstrobemers2 = read_file([], ["randstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(3+k)+"_match_"+str(error)+".out" for k in strobe_range])
    minstrobemers28 = read_file([], ["minstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_match_"+str(error)+".out" for k in strobe_range])
    hybridstrobemers28 = read_file([],["hybridstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_match_"+str(error)+".out" for k in strobe_range])
    randstrobemers28 = read_file([], ["randstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_match_"+str(error)+".out" for k in strobe_range])


    # Plot comparison between all Island size
    fig = plt.figure()
    X = np.arange(len(k_size))

    colors = ["#00ba32","#00d6e7","#fad100","#697ed5","#c76674","#9350a1"]
    plt.xlabel("k")
    plt.xticks(pos, k_size)
    plt.ylabel("Average island size")

    plt.plot(pos, [x[0] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
    plt.plot(pos, [x[0] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
    plt.plot(pos, [x[0] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
    plt.plot(pos, [x[0] for x in minstrobemers2], color = colors[3], label='minstrobemers',linewidth=3.0)
    plt.plot(pos, [x[0] for x in hybridstrobemers2], color = colors[4], label='hybridstrobemers',linewidth=3.0)
    plt.plot(pos, [x[2] for x in randstrobemers2], color = colors[5], label='randstrobemers',linewidth=3.0)

    plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
    plt.savefig("../results/Match_island_"+str(error)+".png",bbox_inches='tight')

    # Plot comparison between all match coverage
    fig = plt.figure()
    X = np.arange(len(k_size))

    colors = ["#00ba32","#00d6e7","#fad100","#697ed5","#c76674","#9350a1"]
    plt.xlabel("k")
    plt.xticks(pos, k_size)
    plt.ylabel("Match coverage")

    plt.plot(pos, [x[2] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
    plt.plot(pos, [x[2] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
    plt.plot(pos, [x[2] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
    plt.plot(pos, [x[2] for x in minstrobemers2], color = colors[3], label='minstrobemers',linewidth=3.0)
    plt.plot(pos, [x[2] for x in hybridstrobemers2], color = colors[4], label='hybridstrobemers',linewidth=3.0)
    plt.plot(pos, [x[2] for x in randstrobemers2], color = colors[5], label='randstrobemers',linewidth=3.0)


    plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
    plt.savefig("../results/Match_cov_"+str(error)+".png",bbox_inches='tight')
