import sys

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

k_mers = [16,20,24,28,32]
errors = [0,1,2]

def get_unique(in_file):
    frequencies = np.fromfile(in_file, dtype=np.uint16)

    number_elements = len (frequencies)
    number_unique = (frequencies == 1).sum()
    # print(number_elements, number_unique, (number_unique*100.0)/number_elements)
    return (number_unique*100.0)/number_elements

def get_results(species):
    results = []
    for e in errors:
        results.append([])
    for k in k_mers:
        for e in errors:
            genmap_file = 'genmap/'+species+'_K_'+str(k)+'_E_'+str(e)+'.freq16'
            results[errors.index(e)].append(get_unique(genmap_file))
    fig = plt.figure()
    X = np.arange(len(k_mers))
    colors = ["#00ba32","#00d6e7","#fad100"] # ["#1b9e77","#d95f02","#7570b3"] #["#00ba32","#00d6e7","#fad100"] 
    pos = [0.25,1.25,2.25,3.25,4.25]
    plt.xlabel("k")
    plt.xticks(pos, k_mers)
    plt.ylabel("percentage of unique k-mers")

    plt.bar(X + 0.00, results[0], color = colors[0], width = 0.25, label='0')
    plt.bar(X + 0.25, results[1], color = colors[1], width = 0.25, label='1')
    plt.bar(X + 0.50, results[2], color = colors[2], width = 0.25, label='2')
    plt.legend(title="# of errors")
    plt.savefig("Uniqueness_"+species+".png")

get_results("Human")
get_results("Mouse")
get_results("Wheat")
