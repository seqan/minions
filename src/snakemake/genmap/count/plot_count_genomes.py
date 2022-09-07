import os

import numpy as np
import matplotlib.pyplot as plt

# Read count files
def get_counts(infile):
    with open(infile, 'r') as f:
        for line in f:
            return int(line.split()[1])

# Different species
species = ["Human", "Mouse", "Wheat"]
ungapped = {"Human":[],"Mouse":[],"Wheat":[]}
gapped_1 = {"Human":[],"Mouse":[],"Wheat":[]} #[36607,933855,14548847,192920575,3758077695]
gapped_2 = {"Human":[],"Mouse":[],"Wheat":[]} #[51755,975475,13954519,241004285,3169577727]

k_mers = [16,20,24,28,32]
k_mers_1 = [36607,933855,14548847,192920575,3758077695]
k_mers_2 = [51755,975475,13954519,241004285,3169577727]
all_k_mers = [4**i for i in k_mers]
all_k_mers_1 = [4**(i-4) for i in k_mers]
all_k_mers_2 = [4**(i-8) for i in k_mers]

# Go over all count files
for i in range(len(k_mers)):
    for spec in species:
        ungapped[spec].append(get_counts("0_"+spec+"_kmer_hash_"+str(k_mers[i])+"_counts.out"))
        gapped_1[spec].append(get_counts(str(k_mers_1[i])+"_"+spec+"_kmer_hash_"+str(k_mers[i])+"_counts.out"))
        gapped_2[spec].append(get_counts(str(k_mers_2[i])+"_"+spec+"_kmer_hash_"+str(k_mers[i])+"_counts.out"))

# Create figures
for spec in species:
    print(spec)
    print(ungapped[spec])
    print(gapped_1[spec])
    print(gapped_2[spec])

    data = [ungapped[spec], gapped_1[spec], gapped_2[spec]]
    X = np.arange(5)
    fig = plt.figure()
    colors = ["#00ba32","#00d6e7","#fad100"]
    pos = [0.25,1.25,2.25,3.25,4.25]
    plt.xlabel("k")
    plt.xticks(pos, k_mers)
    plt.ylabel("# of k-mers")
    plt.bar(X + 0.00, data[0], color = colors[0], width = 0.25,label="0 gaps")
    plt.bar(X + 0.25, data[1], color = colors[1], width = 0.25, label = "4 gaps")
    plt.bar(X + 0.50, data[2], color = colors[2], width = 0.25, label = "8 gaps")
    plt.legend()
    plt.savefig("kmer_set_size_"+spec+".png")

    data = [[ungapped[spec][i]*100.0/all_k_mers[i] for i in range(len(all_k_mers))], [gapped_1[spec][i]*100.0/all_k_mers_1[i] for i in range(len(all_k_mers))], [gapped_2[spec][i]*100.0/all_k_mers_2[i] for i in range(len(all_k_mers))]]
    fig = plt.figure()
    plt.xlabel("k")
    plt.xticks(pos, k_mers)
    plt.ylabel("% of k-mers")
    plt.bar(X + 0.00, data[0], color = colors[0], width = 0.25,label="0 gaps")
    plt.bar(X + 0.25, data[1], color = colors[1], width = 0.25, label = "4 gaps")
    plt.bar(X + 0.50, data[2], color = colors[2], width = 0.25, label = "8 gaps")
    plt.legend()
    plt.savefig("percentage_kmer_set_size_"+spec+".png")
