import matplotlib.pyplot as plt
import numpy as np
import sys

input_file = sys.argv[1]
frequencies = []
# Read Data
with open(input_file, 'rb') as fin:
    byte = fin.read(8)
    while byte: # Ignore minimiser hash, get count value
        frequencies.append(int.from_bytes(fin.read(2), byteorder='big'))
        byte = fin.read(8)

frequencies.sort()
hist, bin_edges = np.histogram(frequencies)
print(hist)
plt.hist(x=frequencies, bins=200, color='#0504aa', alpha=0.7, rwidth=0.85)
plt.show()
