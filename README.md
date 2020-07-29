# Preliminaries
## Disclaimer and Licensing

DMFT++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
DMFT++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with DMFT++. If not, see <https://www.gnu.org/licenses/>.
The full software license for DMFT++ version 1.0.0
can be found in
file LICENSE.

## Please cite this work

DMFT++ is a free and open source
DMFT implementation with a DMRG impurity solver
for strongly correlated electrons.
The full software license for DMFT++ version 0.
can be found in
file LICENSE.
You are welcomed to use it and publish data
obtained with DMFT++. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.

## Mission Statement

DMFT++ is a C++ native application implementing the DMFT algorithm for strongly correlated electron models. 
Features and options are chosen from a user-friendly input file.
The implementation aims to be as fast as possible, so that it compiles natively, uses optimized linear algebra libraries,
the DMRG++ code base and its shared memory parallelization (pthreads), 
and MPI for running different frequencies in parallel in different compute nodes.

## Papers used

TBW

## Code Signature

TBW

# Building and Running DMFT++

## Required Software

* GNU C++
* PsimagLite (see below)
* DMRG++

## Optional Software

* make or gmake (only needed to use the Makefile)
* perl (may be needed to run some auxiliary script)

## Quick Start

1. Use your distribution repository tool to install gcc with support for C++,
make, perl, and git if you don't have them.

2. Issue

```BASH
    $ cd someDirectory/
    $ git clone https://github.com/g1257/PsimagLite.git
	$ git clone https://github.com/g1257/dmrgpp.git
    $ git clone https://github.com/g1257/dmft.git
```

3. Compile PsimagLite

```BASH
    $ cd PsimagLite
    $ git checkout features
    $ git pull origin features
    $ cd lib/
    $ ./configure.pl
    $ make -j something
    $ cd ../../
```

4. Compile libdmrgpp.a

```BASH
    $ cd dmrgpp/src
    $ git checkout features
    $ git pull origin features
    $ ./configure.pl
    $ make libdmrgpp.a libkronutil.a -j something
    $ cd ../../
```
4. Now issue

```BASH
    $ cd DMFTpp
    $ git checkout features
    $ git pull origin features
    $ cd src
    $ ./configure.pl
    $ make clean
    $ make all -j something
```

5. You can run it with

```BASH
    $ ./DMFTpp -f input.ain 
```


