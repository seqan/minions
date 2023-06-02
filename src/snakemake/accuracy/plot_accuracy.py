import sys

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

thresholds = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7]
pos = [x+0.25 for x in range(len(thresholds))]
strobe_range = [10]

def read_file(results, files):
    num_fp = 0
    fn_0 = 0.0
    i = 0
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                tp = int(line.split()[1])
                tn = int(line.split()[2])
                fp = int(line.split()[3])
                fn = int(line.split()[4])
                if (fn == 0):
                    fn_0 = thresholds[i]
                    num_fp = fp
                i += 1
    results.append([num_fp, fn_0])

# Read all files for an error
for error in [2,3,4,5]:
    results = []
    read_file(results, ["0_"+str(error)+"_"+str(threshold)+"_minimiser_hash_"+str(k)+"_"+str(k)+"_all_accuracy.out" for k in [20] for threshold in thresholds])
    read_file(results, ["16252901_"+str(error)+"_"+str(threshold)+"_minimiser_hash_"+str(k)+"_"+str(k)+"_all_accuracy.out" for k in [24] for threshold in thresholds])
    read_file(results, ["180082591_"+str(error)+"_"+str(threshold)+"_minimiser_hash_"+str(k)+"_"+str(k)+"_all_accuracy.out" for k in [28] for threshold in thresholds])

    read_file(results, [str(error)+"_"+str(threshold)+"_minstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(3+k)+"_all_accuracy.out" for k in strobe_range for threshold in thresholds])
    read_file(results,[str(error)+"_"+str(threshold)+"_hybridstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(4+k)+"_all_accuracy.out" for k in strobe_range for threshold in thresholds])
    read_file(results, [str(error)+"_"+str(threshold)+"_randstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(3+k)+"_all_accuracy.out" for k in strobe_range for threshold in thresholds])
    read_file(results, [str(error)+"_"+str(threshold)+"_minstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_all_accuracy.out" for k in strobe_range for threshold in thresholds])
    read_file(results,[str(error)+"_"+str(threshold)+"_hybridstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_all_accuracy.out" for k in strobe_range for threshold in thresholds])
    read_file(results, [str(error)+"_"+str(threshold)+"_randstrobemers_"+str(k)+"_2_"+str(0)+"_"+str(7+k)+"_all_accuracy.out" for k in strobe_range for threshold in thresholds])

    print("Error: ",error, "\n", results)

print("Representative:")
for error in [2,3,4,5]:
    results = []
    read_file(results, ["0_"+str(error)+"_"+str(threshold)+"_minimiser_hash_20_24_all_accuracy.out" for k in [20] for threshold in thresholds])
    #read_file(results, ["0_"+str(error)+"_"+str(threshold)+"_modmer_hash_20_3_all_accuracy.out" for k in [20] for threshold in thresholds])
    read_file(results, [str(error)+"_"+str(threshold)+"_syncmer_hash_20_18_0_0"+"_all_accuracy.out" for k in [20] for threshold in thresholds])
    read_file(results, [str(error)+"_"+str(threshold)+"_syncmer_hash_20_15_0_6"+"_all_accuracy.out" for k in [20] for threshold in thresholds])
    print("Error: ",error, "\n", results)

def fix():
    fig = plt.figure()
    labels = ['k-mer','4 k-mer','8 k-mer', 'minstrobemers','hybridstrobemers','randstrobemers']

    colors = ["#00ba32","#00d6e7","#fad100","#697ed5","#c76674","#9350a1"]
    #plt.xlabel("Threshold")
    plt.xticks(pos, thresholds)
    plt.ylabel("# False Positives")
    y_pos = np.arange(len(labels))

    plt.bar(y_pos, [x[0] for x in results[:6]], align='center')

    plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
    plt.savefig("../results/FPR_"+str(error)+".png",bbox_inches='tight')

    # Plot comparison between all match coverage
    fig = plt.figure()
    X = np.arange(len(thresholds))

    colors = ["#00ba32","#00d6e7","#fad100","#697ed5","#c76674","#9350a1"]
    plt.xlabel("Threshold")
    plt.xticks(pos, thresholds)
    plt.ylabel("False Negative Rate")

    plt.plot(pos, [x[1] for x in kmers], color = colors[0], label='k-mer', linewidth=3.0)
    plt.plot(pos, [x[1] for x in gapped4_kmers], color = colors[1], label='4 k-mer',linewidth=3.0)
    plt.plot(pos, [x[1] for x in gapped8_kmers], color = colors[2], label='8 k-mer',linewidth=3.0)
    plt.plot(pos, [x[1] for x in minstrobemers2], color = colors[3], label='minstrobemers',linewidth=3.0)
    plt.plot(pos, [x[1] for x in hybridstrobemers2], color = colors[4], label='hybridstrobemers',linewidth=3.0)
    plt.plot(pos, [x[1] for x in randstrobemers2], color = colors[5], label='randstrobemers',linewidth=3.0)

    plt.legend(bbox_to_anchor=(1.01, 0.75), title="Methods")
    plt.savefig("../results/FNR_"+str(error)+".png",bbox_inches='tight')
