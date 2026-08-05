# Preliminaries
## Disclaimer and Licensing

The full software license for DMRG++ version 3
can be found in file LICENSE in the root directory of the code.
DMRG++ is a free and open source implementation of the
DMRG algorithm. You are welcomed to use it and publish data
obtained with DMRG++. If you do,
<b>please cite this work</b> (see next subsection).

See DISCLAIMERS in the LICENSE file.

## Please cite this work

cincuenta is a free and open source
DMFT implementation with a DMRG impurity solver
for strongly correlated electrons.
You are welcomed to use it and publish data
obtained with cincuenta. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.

## Mission Statement

cincuenta is a C++ native application implementing the DMFT algorithm for strongly correlated electron models.
Features and options are chosen from a user-friendly input file.
The implementation aims to be as fast as possible, so that it compiles natively, uses optimized linear algebra libraries,
the DMRG++ code base and its shared memory parallelization (pthreads),
and MPI for running different frequencies in parallel in different compute nodes.

## Papers used

TBW

## Code Signature

TBW

# Building and Running cincuenta

## Required Software

* GNU C++
* PsimagLite (see below)
* DMRG++

## Optional Software

* make or gmake (only needed to use the Makefile)
* perl (may be needed to run some auxiliary script)

## Quick Start
Compile DMRG++; instructions are in the top directory README file.

5. You can run it with

```BASH
./cincuenta -f ../cincuenta/TestSuite/inputs/inputCincuenta0.ain
perl ../cincuenta/scripts/extract runForinputCincuenta0.cout
```


