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
    minimiser = read_file([], ["0_minimiser_hash_20_"+str(w)+"_match_"+str(error)+".out" for w in range(24,44,4)])
    modmer = read_file([], ["0_modmer_hash_20_"+str(w)+"_match_"+str(error)+".out" for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["syncmer_hash_20_"+str(w)+"_0_0_match_"+str(error)+".out" for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["syncmer_hash_20_"+str(w)+"_0_6_match_"+str(error)+".out" for w in [15,11,7,3,1]])

    # Plot comparison between all Island size
    fig = plt.figure()
    X = np.arange(len(k_size))

    colors = ["#01d63a","#00e7e0","#fefea1","#748beb"]
    plt.xlabel("k")
    plt.xticks(pos, k_size)
    print(opensyncmer)#plt.ylabel("Average island size")

    plt.plot(pos, [x[0] for x in minimiser], color = colors[0], label='(w,20)-minimizer', linewidth=3.0)
    plt.plot(pos, [x[0] for x in modmer], color = colors[1], label='(20,m)-modmer',linewidth=3.0)
    plt.plot(pos, [x[0] for x in opensyncmer], color = colors[2], label='(20,s,[0],1)-syncmer',linewidth=3.0)
    plt.plot(pos, [x[0] for x in closedsyncmer], color = colors[3], label='(20,s,[0,6],1)-syncmer',linewidth=3.0)

    plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
    plt.savefig("../results/Match_island_representative"+str(error)+".png",bbox_inches='tight')

    # Plot comparison between all match coverage
    fig = plt.figure()
    X = np.arange(len(k_size))

    colors = ["#00ba32","#00d6e7","#fad100","#697ed5","#c76674","#9350a1"]
    plt.xlabel("k")
    plt.xticks(pos, k_size)
    plt.ylabel("Match coverage")

    plt.plot(pos, [x[2] for x in minimiser], color = colors[0], label='(w,20)-minimizer', linewidth=3.0)
    plt.plot(pos, [x[2] for x in modmer], color = colors[1], label='(20,m)-modmer',linewidth=3.0)
    plt.plot(pos, [x[2] for x in opensyncmer], color = colors[2], label='(20,s,[0],1)-syncmer',linewidth=3.0)
    plt.plot(pos, [x[2] for x in closedsyncmer], color = colors[3], label='(20,s,[0,6],1)-syncmer',linewidth=3.0)


    plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
    plt.savefig("../results/Match_cov_representative"+str(error)+".png",bbox_inches='tight')
