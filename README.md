# Comparison between different methods to simplify sequence data

[![Build Status](https://github.com/seqan/app-template/workflows/App%20CI/badge.svg)](https://github.com/seqan/app-template/actions?query=branch%3Amaster+workflow%3A%22App+CI%22)

The aim of this repository is to compare different methods to simplify sequence data in a way that new methods can
easily be added and therefore included in the comparison. Furthermore, this repository is designed in such a way that
all methods can be found in header files in the include folder, so that this repository can be used as a library in
other projects.
Currently, the following methods are supported:

- k-mers
- minimizers

See Issue #1 for a list of methods that will be added in the future.

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

Download the data that is used to perform the comparison.
```
wget https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
tar xfz example_data.tar.gz
```
