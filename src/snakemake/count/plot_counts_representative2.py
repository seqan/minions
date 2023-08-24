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


def create_plot(minimiser, modmer, opensyncmer, closedsyncmer, outfile,number_elem=5):
    # Plot comparison between k-mers
    k_size = [i for i in range(number_elem)]
    pos = [x+0.25 for x in range(len(k_size))]

    fig = plt.figure()
    X = np.arange(len(k_size))
    colors = ["#890015","#5cffca","#a13ff0","#ff9ba0"]
    colors_error = ["#01d63a","#00e7e0","#fefea1","#748beb"]
    plt.xlabel("w,m or s")
    plt.xticks(pos, k_size)
    plt.ylabel("# of submers")

    if (number_elem == 5):
        plt.plot(pos, [x[0] for x in minimiser], color = colors[0], label='(w,20)-minimizer',linewidth=3.0)
        plt.plot(pos, [x[0] for x in modmer], color = colors[1], label='(20,m)-modmer')
        plt.plot(pos, [x[0] for x in opensyncmer], color = colors[2], label='(20,s,[0],1)-syncmer',linewidth=3.0)
        plt.plot(pos, [x[0] for x in closedsyncmer], color = colors[3], label='(20,s,[0,20-s],1)-syncmer',linewidth=3.0)
    else:
        plt.plot(pos, [x[0] for x in minimiser], color = colors[0], label='(w,27)-minimizer',linewidth=3.0)
        plt.plot(pos, [x[0] for x in modmer], color = colors[1], label='(27,m)-modmer',linewidth=3.0)
        plt.plot(pos, [x[0] for x in opensyncmer], color = colors[2], label='(27,s,[0],1)-syncmer',linewidth=3.0)
        plt.plot(pos, [x[0] for x in closedsyncmer], color = colors[3], label='(27,s,[0,27-s],1)-syncmer',linewidth=3.0)

    plt.legend(title="Methods")
    plt.savefig(outfile, bbox_inches='tight')

minimiser = read_file([], ["Rep2_min_2_1_14_Strobemer_minimiser_hash_10_"+str(w)+"_counts.out" for w in [i for i in range(24,44,4)]])
modmer = read_file([], ["Rep2_min_2_1_14_Strobemer_modmer_hash_10_"+str(w)+"_counts.out" for w in [3,5,7,9,11]])
opensyncmer = read_file([], ["min_2_1_14_Strobemer_syncmer_hash_10_"+str(w)+"_0_0_counts.out" for w in [18,16,14,12,10]])
closedsyncmer = read_file([], ["min_2_1_14_Strobemer_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_counts.out" for w in [15,11,7,3,1]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_min1.png")

minimiser = read_file([], ["Rep2_rand_2_1_14_Strobemer_minimiser_hash_10_"+str(w)+"_counts.out" for w in [i for i in range(24,44,4)]])
modmer = read_file([], ["Rep2_rand_2_1_14_Strobemer_modmer_hash_10_"+str(w)+"_counts.out" for w in [3,5,7,9,11]])
opensyncmer = read_file([], ["rand_2_1_14_Strobemer_syncmer_hash_10_"+str(w)+"_0_0_counts.out" for w in [18,16,14,12,10]])
closedsyncmer = read_file([], ["rand_2_1_14_Strobemer_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_counts.out" for w in [15,11,7,3,1]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_rand1.png")

minimiser = read_file([], ["Rep2_hybrid_2_1_14_Strobemer_minimiser_hash_10_"+str(w)+"_counts.out" for w in [i for i in range(24,44,4)]])
modmer = read_file([], ["Rep2_hybrid_2_1_14_Strobemer_modmer_hash_10_"+str(w)+"_counts.out" for w in [3,5,7,9,11]])
opensyncmer = read_file([], ["hybrid_2_1_14_Strobemer_syncmer_hash_10_"+str(w)+"_0_0_counts.out" for w in [18,16,14,12,10]])
closedsyncmer = read_file([], ["hybrid_2_1_14_Strobemer_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_counts.out" for w in [15,11,7,3,1]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_hybrid1.png")

minimiser = read_file([], ["Rep2_min_2_1_18_Strobemer_minimiser_hash_10_"+str(w)+"_counts.out" for w in [i for i in range(24,44,4)]])
modmer = read_file([], ["Rep2_min_2_1_18_Strobemer_modmer_hash_10_"+str(w)+"_counts.out" for w in [3,5,7,9,11]])
opensyncmer = read_file([], ["min_2_1_18_Strobemer_syncmer_hash_10_"+str(w)+"_0_0_counts.out" for w in [18,16,14,12,10]])
closedsyncmer = read_file([], ["min_2_1_18_Strobemer_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_counts.out" for w in [15,11,7,3,1]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_min.png")


minimiser = read_file([], ["Rep2_rand_2_1_18_Strobemer_minimiser_hash_10_"+str(w)+"_counts.out" for w in [i for i in range(24,44,4)]])
modmer = read_file([], ["Rep2_rand_2_1_18_Strobemer_modmer_hash_10_"+str(w)+"_counts.out" for w in [3,5,7,9,11]])
opensyncmer = read_file([], ["rand_2_1_18_Strobemer_syncmer_hash_10_"+str(w)+"_0_0_counts.out" for w in [18,16,14,12,10]])
closedsyncmer = read_file([], ["rand_2_1_18_Strobemer_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_counts.out" for w in [15,11,7,3,1]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_rand.png")

minimiser = read_file([], ["Rep2_hybrid_2_1_18_Strobemer_minimiser_hash_10_"+str(w)+"_counts.out" for w in [i for i in range(24,44,4)]])
modmer = read_file([], ["Rep2_hybrid_2_1_18_Strobemer_modmer_hash_10_"+str(w)+"_counts.out" for w in [3,5,7,9,11]])
opensyncmer = read_file([], ["hybrid_2_1_18_Strobemer_syncmer_hash_10_"+str(w)+"_0_0_counts.out" for w in [18,16,14,12,10]])
closedsyncmer = read_file([], ["hybrid_2_1_18_Strobemer_syncmer_hash_10_"+str(w)+"_0_"+str(20-w)+"_counts.out" for w in [15,11,7,3,1]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_hybrid.png")

minimiser = read_file([], ["Rep2_min_3_1_13_Strobemer_minimiser_hash_9_"+str(w)+"_counts.out" for w in [i for i in range(29,44,4)]])
modmer = read_file([], ["Rep2_min_3_1_13_Strobemer_modmer_hash_9_"+str(w)+"_counts.out" for w in [2,4,6,8]])
opensyncmer = read_file([], ["min_3_1_13_Strobemer_syncmer_hash_9_"+str(w)+"_0_0_counts.out" for w in [26,24,22,20]])
closedsyncmer = read_file([], ["min_3_1_13_Strobemer_syncmer_hash_9_"+str(w)+"_0_"+str(27-w)+"_counts.out" for w in [24,20,16,12]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_min31.png",4)

minimiser = read_file([], ["Rep2_rand_3_1_13_Strobemer_minimiser_hash_9_"+str(w)+"_counts.out" for w in [i for i in range(29,44,4)]])
modmer = read_file([], ["Rep2_rand_3_1_13_Strobemer_modmer_hash_9_"+str(w)+"_counts.out" for w in [2,4,6,8]])
opensyncmer = read_file([], ["rand_3_1_13_Strobemer_syncmer_hash_9_"+str(w)+"_0_0_counts.out" for w in [26,24,22,20]])
closedsyncmer = read_file([], ["rand_3_1_13_Strobemer_syncmer_hash_9_"+str(w)+"_0_"+str(27-w)+"_counts.out" for w in [24,20,16,12]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representatie_rand31.png", 4)

minimiser = read_file([], ["Rep2_hybrid_3_1_13_Strobemer_minimiser_hash_9_"+str(w)+"_counts.out" for w in [i for i in range(29,44,4)]])
modmer = read_file([], ["Rep2_hybrid_3_1_13_Strobemer_modmer_hash_9_"+str(w)+"_counts.out" for w in [2,4,6,8]])
opensyncmer = read_file([], ["hybrid_3_1_13_Strobemer_syncmer_hash_9_"+str(w)+"_0_0_counts.out" for w in [26,24,22,20]])
closedsyncmer = read_file([], ["hybrid_3_1_13_Strobemer_syncmer_hash_9_"+str(w)+"_0_"+str(27-w)+"_counts.out" for w in [24,20,16,12]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_hybrid31.png",4)

minimiser = read_file([], ["Rep2_min_3_1_17_Strobemer_minimiser_hash_9_"+str(w)+"_counts.out" for w in [i for i in range(29,44,4)]])
modmer = read_file([], ["Rep2_min_3_1_17_Strobemer_modmer_hash_9_"+str(w)+"_counts.out" for w in [2,4,6,8]])
opensyncmer = read_file([], ["min_3_1_17_Strobemer_syncmer_hash_9_"+str(w)+"_0_0_counts.out" for w in [26,24,22,20]])
closedsyncmer = read_file([], ["min_3_1_17_Strobemer_syncmer_hash_9_"+str(w)+"_0_"+str(27-w)+"_counts.out" for w in [24,20,16,12]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_min3.png", 4)

minimiser = read_file([], ["Rep2_rand_3_1_17_Strobemer_minimiser_hash_9_"+str(w)+"_counts.out" for w in [i for i in range(29,44,4)]])
modmer = read_file([], ["Rep2_rand_3_1_17_Strobemer_modmer_hash_9_"+str(w)+"_counts.out" for w in [2,4,6,8]])
opensyncmer = read_file([], ["rand_3_1_17_Strobemer_syncmer_hash_9_"+str(w)+"_0_0_counts.out" for w in [26,24,22,20]])
closedsyncmer = read_file([], ["rand_3_1_17_Strobemer_syncmer_hash_9_"+str(w)+"_0_"+str(27-w)+"_counts.out" for w in [24,20,16,12]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_rand3.png", 4)

minimiser = read_file([], ["Rep2_hybrid_3_1_17_Strobemer_minimiser_hash_9_"+str(w)+"_counts.out" for w in [i for i in range(29,44,4)]])
modmer = read_file([], ["Rep2_hybrid_3_1_17_Strobemer_modmer_hash_9_"+str(w)+"_counts.out" for w in [2,4,6,8]])
opensyncmer = read_file([], ["hybrid_3_1_17_Strobemer_syncmer_hash_9_"+str(w)+"_0_0_counts.out" for w in [26,24,22,20]])
closedsyncmer = read_file([], ["hybrid_3_1_17_Strobemer_syncmer_hash_9_"+str(w)+"_0_"+str(27-w)+"_counts.out" for w in [24,20,16,12]])
create_plot(minimiser, modmer, opensyncmer, closedsyncmer, "../results/Count_representative_hybrid3.png",4)


# Plot Uniqueness
minimiser_hybrid1 = []
modmer_hybrid1 = []
opensyncmer_hybrid1 = []
closedsyncmer_hybrid1 = []

minimiser_hybrid = []
modmer_hybrid = []
opensyncmer_hybrid = []
closedsyncmer_hybrid = []

minimiser_hybrid31 = []
modmer_hybrid31 = []
opensyncmer_hybrid31 = []
closedsyncmer_hybrid31 = []

minimiser_hybrid3 = []
modmer_hybrid3 = []
opensyncmer_hybrid3 = []
closedsyncmer_hybrid3 = []

minimiser_rand1 = []
modmer_rand1 = []
opensyncmer_rand1 = []
closedsyncmer_rand1 = []

minimiser_rand = []
modmer_rand = []
opensyncmer_rand = []
closedsyncmer_rand = []

minimiser_rand31 = []
modmer_rand31 = []
opensyncmer_rand31 = []
closedsyncmer_rand31 = []

minimiser_rand3 = []
modmer_rand3 = []
opensyncmer_rand3 = []
closedsyncmer_rand3 = []

minimiser_min1 = []
modmer_min1 = []
opensyncmer_min1 = []
closedsyncmer_min1 = []

minimiser_min = []
modmer_min = []
opensyncmer_min = []
closedsyncmer_min = []

minimiser_min31 = []
modmer_min31 = []
opensyncmer_min31 = []
closedsyncmer_min31 = []

minimiser_min3 = []
modmer_min3 = []
opensyncmer_min3 = []
closedsyncmer_min3 = []

it = 0
with open("../results/Unique_representative2.out", 'r') as f:
    for line in f:
        number = float(line.split()[1])
        if ("2_1_18_" in line):
            if ("hybrid" in line):
                if ("minimiser" in line):
                    minimiser_hybrid.append(number)
                elif ("modmer" in line):
                    modmer_hybrid.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_hybrid.append(number)
                else:
                    closedsyncmer_hybrid.append(number)
            elif ("rand" in line):
                if ("minimiser" in line):
                    minimiser_rand.append(number)
                elif ("modmer" in line):
                    modmer_rand.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_rand.append(number)
                else:
                    closedsyncmer_rand.append(number)
            elif ("min_" in line):
                if ("minimiser" in line):
                    minimiser_min.append(number)
                elif ("modmer" in line):
                    modmer_min.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_min.append(number)
                else:
                    closedsyncmer_min.append(number)
        elif ("3_1_17_" in line):
            if ("hybrid" in line):
                if ("minimiser" in line):
                    minimiser_hybrid3.append(number)
                elif ("modmer" in line):
                    modmer_hybrid3.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_hybrid3.append(number)
                else:
                    closedsyncmer_hybrid3.append(number)
            elif ("rand" in line):
                if ("minimiser" in line):
                    minimiser_rand3.append(number)
                elif ("modmer" in line):
                    modmer_rand3.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_rand3.append(number)
                else:
                    closedsyncmer_rand3.append(number)
            elif ("min_" in line):
                if ("minimiser" in line):
                    minimiser_min3.append(number)
                elif ("modmer" in line):
                    modmer_min3.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_min3.append(number)
                else:
                    closedsyncmer_min3.append(number)
        elif ("3_1_13_" in line):
            if ("hybrid" in line):
                if ("minimiser" in line):
                    minimiser_hybrid31.append(number)
                elif ("modmer" in line):
                    modmer_hybrid31.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_hybrid31.append(number)
                else:
                    closedsyncmer_hybrid31.append(number)
            elif ("rand" in line):
                if ("minimiser" in line):
                    minimiser_rand31.append(number)
                elif ("modmer" in line):
                    modmer_rand31.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_rand31.append(number)
                else:
                    closedsyncmer_rand31.append(number)
            elif ("min_" in line):
                if ("minimiser" in line):
                    minimiser_min31.append(number)
                elif ("modmer" in line):
                    modmer_min31.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_min31.append(number)
                else:
                    closedsyncmer_min31.append(number)
        else:
            if ("hybrid" in line):
                if ("minimiser" in line):
                    minimiser_hybrid1.append(number)
                elif ("modmer" in line):
                    modmer_hybrid1.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_hybrid1.append(number)
                else:
                    closedsyncmer_hybrid1.append(number)
            elif ("rand" in line):
                if ("minimiser" in line):
                    minimiser_rand1.append(number)
                elif ("modmer" in line):
                    modmer_rand1.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_rand1.append(number)
                else:
                    closedsyncmer_rand1.append(number)
            elif ("min_" in line):
                if ("minimiser" in line):
                    minimiser_min1.append(number)
                elif ("modmer" in line):
                    modmer_min1.append(number)
                elif ("_0_0_" in line):
                    opensyncmer_min1.append(number)
                else:
                    closedsyncmer_min1.append(number)
        it += 1
print(modmer)
# Plot comparison between k-mers

def plot_unique(minimiser, modmer, opensyncmer, closedsyncmer, outfile, num_elem = 5):
    k_size = [i for i in range(num_elem)]
    pos = [x+0.25 for x in range(len(k_size))]

    fig = plt.figure()
    X = np.arange(len(k_size))
    colors = ["#890015","#5cffca","#a13ff0","#ff9ba0"]
    colors_error = ["#01d63a","#00e7e0","#fefea1","#748beb"]
    plt.xlabel("w,m or s")
    plt.xticks(pos, k_size)
    plt.ylabel("% of unique submers") # in microseconds

    if (num_elem == 5):
        plt.plot(pos, minimiser[:5], color = colors[0], label='(w,20)-minimizer',linewidth=3.0)
        plt.plot(pos, modmer[:5], color = colors[1], label='(20,m)-modmer',linewidth=3.0)
        plt.plot(pos, opensyncmer[:5], color = colors[2], label='(20,s,[0],1)-syncmer',linewidth=3.0)
        plt.plot(pos, closedsyncmer[:5], color = colors[3], label='(20,s,[0,20-s],1)-syncmer',linewidth=3.0)
    else:
        plt.plot(pos, minimiser[:5], color = colors[0], label='(w,27)-minimizer',linewidth=3.0)
        plt.plot(pos, modmer[:5], color = colors[1], label='(27,m)-modmer',linewidth=3.0)
        plt.plot(pos, opensyncmer[:5], color = colors[2], label='(27,s,[0],1)-syncmer',linewidth=3.0)
        plt.plot(pos, closedsyncmer[:5], color = colors[3], label='(27,s,[0,27-s],1)-syncmer',linewidth=3.0)

    plt.legend(title="Methods")
    plt.savefig(outfile, bbox_inches='tight')

plot_unique(minimiser_min, modmer_min,opensyncmer_min, closedsyncmer_min, "../results/Unique_Representative_min.png")
plot_unique(minimiser_min1, modmer_min1,opensyncmer_min1, closedsyncmer_min1, "../results/Unique_Representative_min1.png")
plot_unique(minimiser_rand, modmer_min,opensyncmer_rand, closedsyncmer_rand, "../results/Unique_Representative_rand.png")
plot_unique(minimiser_rand1, modmer_min1,opensyncmer_rand1, closedsyncmer_rand1, "../results/Unique_Representative_rand1.png")
plot_unique(minimiser_hybrid, modmer_min,opensyncmer_hybrid, closedsyncmer_hybrid, "../results/Unique_Representative_hybrid.png")
plot_unique(minimiser_hybrid1, modmer_min1,opensyncmer_hybrid1, closedsyncmer_hybrid1, "../results/Unique_Representative_hybrid1.png")

plot_unique(minimiser_min3, modmer_min3,opensyncmer_min3, closedsyncmer_min3, "../results/Unique_Representative_min3.png",4)
plot_unique(minimiser_min31, modmer_min31,opensyncmer_min31, closedsyncmer_min31, "../results/Unique_Representative_min31.png",4)
plot_unique(minimiser_rand3, modmer_min3,opensyncmer_rand3, closedsyncmer_rand3, "../results/Unique_Representative_rand3.png",4)
plot_unique(minimiser_rand31, modmer_min31,opensyncmer_rand31, closedsyncmer_rand31, "../results/Unique_Representative_rand31.png",4)
plot_unique(minimiser_hybrid3, modmer_min3,opensyncmer_hybrid3, closedsyncmer_hybrid3, "../results/Unique_Representative_hybrid3.png",4)
plot_unique(minimiser_hybrid31, modmer_min31,opensyncmer_hybrid31, closedsyncmer_hybrid31, "../results/Unique_Representative_hybrid31.png",4)


minimiser = []
modmer = []
opensyncmer = []
closedsyncmer = []
it = 0
with open("../results/Unique_representative.out", 'r') as f:
    for line in f:
        number = float(line.split()[1])
        if (it < 10):
            minimiser.append(number)
        elif (it < 20):
            modmer.append(number)
        elif (it < 30):
            opensyncmer.append(number)
        elif (it < 40):
            closedsyncmer.append(number)
        it += 1

minimiser_all = [np.mean(minimiser),np.mean(minimiser_min1),np.mean(minimiser_hybrid1),np.mean(minimiser_rand1)]
modmer_all = [np.mean(modmer),np.mean(modmer_min1),np.mean(modmer_hybrid1),np.mean(modmer_rand1)]
opensyncmer_all = [np.mean(opensyncmer),np.mean(opensyncmer_min1),np.mean(opensyncmer_hybrid1),np.mean(opensyncmer_rand1)]
closedsyncmer_all = [np.mean(closedsyncmer),np.mean(closedsyncmer_min1),np.mean(closedsyncmer_hybrid1),np.mean(closedsyncmer_rand1)]

def plot_bar(minimiser_all, modmer_all, opensyncmer_all, closedsyncmer_all, outfile):
    X = np.arange(4)
    colors = ["#890015","#5cffca","#a13ff0","#ff9ba0"]
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.bar(X + 0.00, minimiser_all, color = colors[0], label='minimizer', width = 0.2)
    ax.bar(X + 0.2, modmer_all, color = colors[1], label='modmer', width = 0.2)
    ax.bar(X + 0.4, opensyncmer_all, color = colors[2], label='syncmer', width = 0.2)
    ax.bar(X + 0.6, closedsyncmer_all, color = colors[3], label='syncmer', width = 0.2)
    ax.set_xticks([0.3,1.3,2.3,3.3])
    ax.set_xticklabels(["k-mer","min", "hybrid", "rand"])
    plt.ylabel("% of unique submers")

    plt.legend(title="Methods", bbox_to_anchor=(1.01, 0.65))
    plt.savefig(outfile,bbox_inches='tight')

minimiser_all = [np.mean(minimiser),np.mean(minimiser_min1),np.mean(minimiser_hybrid1),np.mean(minimiser_rand1)]
modmer_all = [np.mean(modmer),np.mean(modmer_min1),np.mean(modmer_hybrid1),np.mean(modmer_rand1)]
opensyncmer_all = [np.mean(opensyncmer),np.mean(opensyncmer_min1),np.mean(opensyncmer_hybrid1),np.mean(opensyncmer_rand1)]
closedsyncmer_all = [np.mean(closedsyncmer),np.mean(closedsyncmer_min1),np.mean(closedsyncmer_hybrid1),np.mean(closedsyncmer_rand1)]
plot_bar(minimiser_all, modmer_all, opensyncmer_all, closedsyncmer_all, "../results/Unique_representative_all_bar1.png")

minimiser_all = [np.mean(minimiser),np.mean(minimiser_min),np.mean(minimiser_hybrid),np.mean(minimiser_rand)]
modmer_all = [np.mean(modmer),np.mean(modmer_min),np.mean(modmer_hybrid),np.mean(modmer_rand)]
opensyncmer_all = [np.mean(opensyncmer),np.mean(opensyncmer_min),np.mean(opensyncmer_hybrid),np.mean(opensyncmer_rand)]
closedsyncmer_all = [np.mean(closedsyncmer),np.mean(closedsyncmer_min),np.mean(closedsyncmer_hybrid),np.mean(closedsyncmer_rand)]
plot_bar(minimiser_all, modmer_all, opensyncmer_all, closedsyncmer_all, "../results/Unique_representative_all_bar.png")

minimiser_all = [np.mean(minimiser),np.mean(minimiser_min3),np.mean(minimiser_hybrid3),np.mean(minimiser_rand3)]
modmer_all = [np.mean(modmer),np.mean(modmer_min3),np.mean(modmer_hybrid3),np.mean(modmer_rand3)]
opensyncmer_all = [np.mean(opensyncmer),np.mean(opensyncmer_min3),np.mean(opensyncmer_hybrid3),np.mean(opensyncmer_rand3)]
closedsyncmer_all = [np.mean(closedsyncmer),np.mean(closedsyncmer_min3),np.mean(closedsyncmer_hybrid3),np.mean(closedsyncmer_rand3)]
plot_bar(minimiser_all, modmer_all, opensyncmer_all, closedsyncmer_all, "../results/Unique_representative_all_bar3.png")

minimiser_all = [np.mean(minimiser),np.mean(minimiser_min31),np.mean(minimiser_hybrid31),np.mean(minimiser_rand31)]
modmer_all = [np.mean(modmer),np.mean(modmer_min31),np.mean(modmer_hybrid31),np.mean(modmer_rand31)]
opensyncmer_all = [np.mean(opensyncmer),np.mean(opensyncmer_min31),np.mean(opensyncmer_hybrid31),np.mean(opensyncmer_rand31)]
closedsyncmer_all = [np.mean(closedsyncmer),np.mean(closedsyncmer_min31),np.mean(closedsyncmer_hybrid31),np.mean(closedsyncmer_rand31)]
plot_bar(minimiser_all, modmer_all, opensyncmer_all, closedsyncmer_all, "../results/Unique_representative_all_bar31.png")
