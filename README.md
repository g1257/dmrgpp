
# Introduction

PsimagLite is free software (see file LICENSE)
Parts might have its own License. See Parts of PsimagLite below.

Please cite PsimagLite if you base any scientific
publication on this software. Citation should read:
G. Alvarez, (2011), PsimagLite (version 1.0)
[computer software], Oak Ridge National Laboratory.

-------------------------------------------------------------------------------

# DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

# Description

PsimagLite is a collection of C++ classes that are common to
 codes for the simulation of strongly correlated electrons.
PsimagLite is inspired in T.S.'s Psimag software (but PsimagLite is not a fork of Psimag).

The reason for PsimagLite is to share code among different applications.
Applications that depend on PsimagLite are:
SpinPhononFermion, DMRG++, Lanczos++, FreeFermions, GpusDoneRight, BetheAnsatz

# Code integrity

Hash of the latest commit is also posted at
https://g1257.github.com/hashes.html

Latest commit should always be signed.
Keys at https://g1257.github.com/keys.html

# Parts of PsimagLite

[DenseLinearAlgebra]
BLAS // wrapper
LAPACK // wrapper
Matrix /// a matrix class

[SparseLinearAlgebra]
CrsMatrix
TridiagonalMatrix
LanczosSolver
LanczosVectors
ChebyshevSolver
SparseRow       <-- slower, consumes less memory
SparseRowCached <-- faster, consumes more memory

[Io] Input output support
IoSimple
ChebyshevSerializer

[Concurrency] To write the same code for serial and MPI (and in
the future pthreads and GPUs)
Concurrency
ConcurrencyMpi
ConcurrencySerial
Pthreads
NoPthreads
PackIndices
Range

[SystemInfo] Basic time/date capability, os, hostname, compiler info
HostInfo
MemoryUsage
Rusage // Rusage class is deprecated, use MemoryUsage instead
Profiling // Profiling through constructor/destructor paradigm as done by M.S in DCA++
It should be called actually scope
ProgressIndicator
Tokenizer
TypeToString
LineMarker

[Math]
LinearPrediction
Minimizer
PlotParams
Sort
Random48
RandomForTests
AkimaSpline
GslWrapper
Fermi
AlmostEqual
BitManip
ContinuedFractionCollection
ContinuedFraction
ChebyshevFunction

[STL extensions] These add operations to std classes and put those
extensions in the std namespace, no new classes here!
Complex
Vector
Stack

[TestSuite]
See README under TestSuite

[PsimagDoc]
See README under PsimagDoc

[loki]
////////////////////////////////////////////////////////////////////////////////
// The Loki Library
// Copyright (c) 2001 by Andrei Alexandrescu
// This code accompanies the book:
// Alexandrescu, Andrei. "Modern C++ Design: Generic Programming and Design
//     Patterns Applied". Copyright (c) 2001. Addison-Wesley.
// Permission to use, copy, modify, distribute and sell this software for any
//     purpose is hereby granted without fee, provided that the above copyright
//     notice appear in all copies and that both that copyright notice and this
//     permission notice appear in supporting documentation.
// The author or Addison-Welsey Longman make no representations about the
//     suitability of this software for any purpose. It is provided "as is"
//     without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////

-------------------------------------------------------------------------------

