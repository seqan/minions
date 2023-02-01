# Comparison between different methods to simplify sequence data

[![Build Status](https://github.com/seqan/app-template/workflows/App%20CI/badge.svg)](https://github.com/seqan/app-template/actions?query=branch%3Amaster+workflow%3A%22App+CI%22) [![codecov](https://codecov.io/gh/seqan/minions/branch/master/graph/badge.svg?token=SJVMYRUKW2)](https://codecov.io/gh/seqan/minions)

The aim of this repository is to compare different methods to simplify sequence data in a way that new methods can
easily be added and therefore included in the comparison. Furthermore, this repository is designed in such a way that
all methods can be found in header files in the include folder, so that this repository can be used as a library in
other projects.
Currently, the following methods are supported:

- k-mers
- minimizers
- modmers
- strobemers (integrated as submodule from [here](https://github.com/ksahlin/strobemers) and implemented as a view)
- syncmers

See Issue #1 for a list of methods that will be added in the future (see down below here for an example usage of each method).

The following evaluation metrics are implemented at the moment (see an example usage for each metric down below):

- speed

# Download
```
git clone --recurse-submodules https://github.com/seqan/minions.git
mkdir build-minions && cd build-minions
cmake ../minions
make
```
Run test to check, if Comparison is working as intended. All tests should pass.
```
make test
```

# Speed

Speeds creates a file called `{method}_speed.out` and returns the speed of processing a singular sequence in microseconds. As typical one sequence file contains multiple sequences the minimum speed, the mean, the variance and the maximum speed are returned. Speed can also handle multiple files and calculate the mean over all sequences found in all files. Speed considers for all supported methods the non-canonical version.

Example usage for calculating the k-mers of a given input file `in.fa`:
```
minions speed --method kmer -k 16 in.fasta
```

This results in the file `kmer_hash_16_speed.out`, which looks like:
```
kmer_hash_16	10	11.478	0.970317	21	-1590685541
```

The first number is the minimum, then follows the mean, the variance and the maximum. The last number in the row can be ignored as it's only used for internal purposes.

**Note:**
Currently, speed supports two implementation of the strobemers. The original one from [Kristoffer Sahlin](https://github.com/ksahlin/strobemers) and the one here presented. The one here presented is more comparable to the other methods used here, because they are based on the same hash functions. Therefore, these strobemers are used for every other evaluation metric.

For the original implementation, add the flag `--original` and note that for the original implementation, only randstrobemers are supported for order 2 and 3, minstrobemers and hybridstrobemers only support order 2. Furthermore, the flags `--w-min` and `--w-max` have different meanings between the original implementation and the implementation here.

`w-min` in the implementation from minions is the distance between the first strobe to second strobe. While for the original implementation, it is the starting position in the sequence of the window that is considered for the second strobe. Therefore, the call with original should always add (k+1) to `w-min` compared to the minion implementation.

`w-max` in the implementation from minions is the window length that should be considered for every strobe besides the first one. All strobes need to be completely inside this window length to be considered. While for the original implementation, it is the position in the sequence until which a strobe that is considered has to start. Therefore, for a strobemer with a strobe length of 8, `w-min` of 0 and `w-max` of 15 in the minion implementation would equal a `w-min` of 9 and `w-max` of 17. For more details, please read the documentation for both implementations.

# Methods

If a metric supports a method, pick it with the flag `--method`.

## k-mers

k-mers are defined by their value of k, which can be defined with `-k`. A gapped k-mer can be used by defining a shape with `--shape`. Shape expects a number, this number should be the decimal representation of a binary number with a starting and an ending 1, each 0 in the binary number will be considered a gap.  

## strobemers

Currently, minstrobemers with the flag `--min`, hybridstrobmers with the flag `--hybrid` and randstrobmers with the flag `--rand` are supported for order 2 and 3, which can be defined with `--order`.
With `-k` the length of a single strobe can be defined and with `--w-min` the distance between the first strobe to the next one and with `--w-max` the length of the window to pick the second and third strobe.

## minimizers

Minimizers support ungapped and gapped k-mers. A window size can be given with `-w`. The randomization of the order is achieved by XOR all k-mer hash values with a seed, if the lexicographical order is wanted `--seed` should be set to 0. For more information, see the [seqan tutorial](http://docs.seqan.de/seqan/3-master-user/tutorial_minimiser.html).

## modmers

Minimizers support ungapped and gapped k-mers. The mod value can be given with `-w`. The randomization of the order is achieved by XOR all k-mer hash values with a seed, if the lexicographical order is wanted `--seed` should be set to 0.

## syncmers

Syncmers support ungapped. The s-mer value can be given with `-w`. The positions of a s-mer that make a k-mer a syncmer can be given with `-p`. The randomization of the order is achieved by XOR all k-mer hash values with a seed, if the lexicographical order is wanted `--seed` should be set to 0.
