import os
import sys
import random

infile = sys.argv[1]
outfile = sys.argv[2]
error = int(sys.argv[3])

def mutated_seq(sequence):
    positions = random.sample(range(len(sequence)), error)
    seqs = list(sequence)
    for pos in positions:
        seqs[pos] = random.choice([i for i in "ACGT" if (i != sequence[pos])])
    return "".join(seqs)

with open(outfile, 'w') as o:
    with open(infile, 'r') as f:
        for line in f:
            if (line[0] == '>'):
                o.write(line)
            else:
                o.write(mutated_seq(line.strip()))
                o.write("\n")
