# Simulated Data

The file "../results/simulated_10000.fa.gz" contains 1,000 random sequences of length 10,000 bp and was created via {mason](https://www.seqan.de/apps/mason.html) with the following command:

```
mason_genome -o simulated_10000.fa $(for i in {1..100}; do echo -l 10000; done;) -s 42
```
