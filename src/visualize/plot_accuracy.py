import matplotlib.pyplot as plt
import sys


def get_labels(name, method):
    if (name == "kmer"):
        return method.split("_")[2]
    elif ((name == "minimiser") | (name == "modmer") | (name == "syncmer")):
        return "("+method.split("_")[2]+","+method.split("_")[3]+")"

name = sys.argv[1]
name2 = sys.argv[2]
all_fp = []
all_fn = []
labels = []
for i in range(3, len(sys.argv)):
    fp = 0
    fn = 0
    with open(sys.argv[i]) as f:
        for line in f:
            fp = int(line.split()[3])
            fn = int(line.split()[4])
            all_fp.append(fp)
            all_fn.append(fn)
            labels.append(get_labels(name, line.split()[0]))

fig = plt.figure()
ax = plt.axes()

plt.title("Accuracy")
plt.xlabel(name)
plt.ylabel("Hits")
ticks = [j for j in range(len(all_fp))]
plt.xticks(ticks, labels)
x = [j for j in range(len(all_fp))]
ax.plot(x, all_fp, color='blue', label = "False Positives")
ax.plot(x, all_fn, color='black', label = "False Negatives")
plt.legend()

plt.savefig("Plot_"+name+"_"+name2+".png")
