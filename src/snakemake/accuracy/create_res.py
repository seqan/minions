import math

exp = 0
with open("search_results.out",'w') as o:
    for i in range(1048576):
        o.write(str(i))
        o.write('\t')
        o.write(str(math.floor(i/16384)))
        o.write('\n')
        if ((exp / 16384)  == 0):
            print(i)
            exp=exp+1
