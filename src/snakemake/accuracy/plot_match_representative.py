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

def read_file(results, files,normalization = False):
    cov = 0.0
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                if (line[:7] == "Matches"):
                    matches = int(line.split(' ')[1].split('\t')[0])
                    missed = int(line.split('\t')[1].split(' ')[1].strip())
                if (line[:7]=="Match C"):
                    cov = round(float(line.split()[2]),2)
                if (line[:7]=="Islands"):
                    mean = round(float(line.split('\t')[1]),2)
                    stdev = round(float(line.split('\t')[2]),2)
                    if (normalization):
                        results.append((mean,stdev,cov/(matches)))
                    else:
                        results.append((mean,stdev,cov))
    return results

def plot_match(minimiser, modmer, opensyncmer, closedsyncmer, outfile1, outfile2):
    # Plot comparison between all Island size
    fig = plt.figure()
    X = np.arange(len(k_size))

    colors = ["#890015","#5cffca","#a13ff0","#ff9ba0"]
    plt.xlabel("k")
    plt.xticks(pos, k_size)
    plt.ylabel("Average island size")

    plt.plot(pos, [x[0] for x in minimiser], color = colors[0], label='(w,20)-minimizer', linewidth=3.0)
    plt.plot(pos, [x[0] for x in modmer], color = colors[1], label='(20,m)-modmer',linewidth=3.0)
    plt.plot(pos, [x[0] for x in opensyncmer], color = colors[2], label='(20,s,[0],1)-syncmer',linewidth=3.0)
    plt.plot(pos, [x[0] for x in closedsyncmer], color = colors[3], label='(20,s,[0,20-s],1)-syncmer',linewidth=3.0)

    plt.legend(loc="upper left", title="Methods")
    plt.savefig(outfile1,bbox_inches='tight')

    # Plot comparison between all match coverage
    fig = plt.figure()
    X = np.arange(len(k_size))
    plt.xlabel("k")
    plt.xticks(pos, k_size)
    plt.ylabel("Match coverage") # /# of matches

    plt.plot(pos, [x[2] for x in minimiser], color = colors[0], label='(w,20)-minimizer', linewidth=3.0)
    plt.plot(pos, [x[2] for x in modmer], color = colors[1], label='(20,m)-modmer',linewidth=3.0)
    plt.plot(pos, [x[2] for x in opensyncmer], color = colors[2], label='(20,s,[0],1)-syncmer',linewidth=3.0)
    plt.plot(pos, [x[2] for x in closedsyncmer], color = colors[3], label='(20,s,[0,20-s],1)-syncmer',linewidth=3.0)


    plt.legend(loc="lower left", title="Methods")
    plt.savefig(outfile2,bbox_inches='tight')


# Read all files for an error
for error in [1,2,5,10]:
    minimiser = read_file([], ["0_minimiser_hash_20_"+str(w)+"_match_"+str(error)+".out" for w in range(24,44,4)])
    modmer = read_file([], ["0_modmer_hash_20_"+str(w)+"_match_"+str(error)+".out" for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["syncmer_hash_20_"+str(w)+"_0_0_match_"+str(error)+".out" for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["syncmer_hash_20_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for w in [15,11,7,3,1]])

    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative"+str(error)+".png","../results/Match_cov_representative"+str(error)+".png")

    minimiser = read_file([], ["0_minimiser_hash_20_"+str(w)+"_match_"+str(error)+".out" for w in range(24,44,4)], True)
    modmer = read_file([], ["0_modmer_hash_20_"+str(w)+"_match_"+str(error)+".out" for w in [3,5,7,9,11]], True)
    opensyncmer = read_file([],["syncmer_hash_20_"+str(w)+"_0_0_match_"+str(error)+".out" for w in [18,16,14,12,10]], True)
    closedsyncmer = read_file([], ["syncmer_hash_20_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for w in [15,11,7,3,1]], True)
    print(minimiser)
    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative"+str(error)+".png","../results/Match_cov_representative_corrected"+str(error)+".png")

    minimiser = read_file([], ["hybridstrobemers_2_" +str(1)+"_"+str(4+k)+"_minimiser_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in range(24,44,4)])
    modmer = read_file([], ["hybridstrobemers_2_" +str(1)+"_"+str(4+k)+"_modmer_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["hybridstrobemers_2_" +str(1)+"_"+str(4+k)+"_syncmer_hash_10_"+str(w)+"_0_0_match_"+str(error)+".out" for k in [10] for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["hybridstrobemers_2_" +str(1)+"_"+str(4+k)+"_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for k in [10] for w in [15,11,7,3,1]])
    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative_hybrid1_"+str(error)+".png","../results/Match_cov_representative_hybrid1_"+str(error)+".png")

    minimiser = read_file([], ["minstrobemers_2_" +str(1)+"_"+str(4+k)+"_minimiser_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in range(24,44,4)])
    modmer = read_file([], ["minstrobemers_2_" +str(1)+"_"+str(4+k)+"_modmer_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["minstrobemers_2_" +str(1)+"_"+str(4+k)+"_syncmer_hash_10_"+str(w)+"_0_0_match_"+str(error)+".out" for k in [10] for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["minstrobemers_2_" +str(1)+"_"+str(4+k)+"_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for k in [10] for w in [15,11,7,3,1]])
    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative_min1_"+str(error)+".png","../results/Match_cov_representative_min1_"+str(error)+".png")

    minimiser = read_file([], ["randstrobemers_2_" +str(1)+"_"+str(4+k)+"_minimiser_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in range(24,44,4)])
    modmer = read_file([], ["randstrobemers_2_" +str(1)+"_"+str(4+k)+"_modmer_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["randstrobemers_2_" +str(1)+"_"+str(4+k)+"_syncmer_hash_10_"+str(w)+"_0_0_match_"+str(error)+".out" for k in [10] for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["randstrobemers_2_" +str(1)+"_"+str(4+k)+"_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for k in [10] for w in [15,11,7,3,1]])
    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative_rand1_"+str(error)+".png","../results/Match_cov_representative_rand1_"+str(error)+".png")

    minimiser = read_file([], ["hybridstrobemers_2_" +str(1)+"_"+str(8+k)+"_minimiser_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in range(24,44,4)])
    modmer = read_file([], ["hybridstrobemers_2_" +str(1)+"_"+str(8+k)+"_modmer_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["hybridstrobemers_2_" +str(1)+"_"+str(8+k)+"_syncmer_hash_10_"+str(w)+"_0_0_match_"+str(error)+".out" for k in [10] for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["hybridstrobemers_2_" +str(1)+"_"+str(8+k)+"_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for k in [10] for w in [15,11,7,3,1]])
    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative_hybrid_"+str(error)+".png","../results/Match_cov_representative_hybrid_"+str(error)+".png")

    minimiser = read_file([], ["minstrobemers_2_" +str(1)+"_"+str(8+k)+"_minimiser_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in range(24,44,4)])
    modmer = read_file([], ["minstrobemers_2_" +str(1)+"_"+str(8+k)+"_modmer_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["minstrobemers_2_" +str(1)+"_"+str(8+k)+"_syncmer_hash_10_"+str(w)+"_0_0_match_"+str(error)+".out" for k in [10] for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["minstrobemers_2_" +str(1)+"_"+str(8+k)+"_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for k in [10] for w in [15,11,7,3,1]])
    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative_min_"+str(error)+".png","../results/Match_cov_representative_min_"+str(error)+".png")

    minimiser = read_file([], ["randstrobemers_2_" +str(1)+"_"+str(8+k)+"_minimiser_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in range(24,44,4)])
    modmer = read_file([], ["randstrobemers_2_" +str(1)+"_"+str(8+k)+"_modmer_hash_10_"+str(w)+"_match_"+str(error)+".out" for k in [10] for w in [3,5,7,9,11]])
    opensyncmer = read_file([],["randstrobemers_2_" +str(1)+"_"+str(8+k)+"_syncmer_hash_10_"+str(w)+"_0_0_match_"+str(error)+".out" for k in [10] for w in [18,16,14,12,10]])
    closedsyncmer = read_file([], ["randstrobemers_2_" +str(1)+"_"+str(8+k)+"_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_match_"+str(error)+".out" for k in [10] for w in [15,11,7,3,1]])
    plot_match(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Match_island_representative_rand_"+str(error)+".png","../results/Match_cov_representative_rand_"+str(error)+".png")
