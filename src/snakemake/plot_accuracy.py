import matplotlib.pyplot as plt
import sys

name = sys.argv[1]
all_fp = []
all_fn = []
for i in range(2, len(sys.argv)):
    fp = 0
    fn = 0
    with open(sys.argv[i]) as f:
        for line in f:
            fp = int(line.split()[2])
            fn = int(line.split()[3])
            all_fp.append(fp)
            all_fn.append(fn)

fig = plt.figure()
ax = plt.axes()

plt.title(name)
plt.xlabel("Method")
plt.ylabel("Hits")
# plt.xticks(ticks, labels)
x = [j for j in range(len(all_fp))]
ax.plot(x, all_fp, color='blue', label = "False Positives")
ax.plot(x, all_fn, color='black', label = "False Negatives")
plt.legend()

plt.savefig("Plot_"+name+".png")
