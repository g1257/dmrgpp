# Quick Start

## Licensing


 The full software license for DMRG++ version 2.0.0
 can be found in file LICENSE in the root directory of the code.
 DMRG++ is a free and open source implementation of the
 DMRG algorithm. You are welcomed to use it and publish data
 obtained with DMRG++. If you do,
<b>please cite this work</b> (see next subsection).

## DISCLAIMER

<small>
<pre>
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
</pre>
</small>

## How To Cite This Work

<small>
<pre>
\@article{re:alvarez0209,
author="G. Alvarez",
title="The Density Matrix Renormalization Group for
Strongly Correlated Electron Systems: A Generic Implementation",
journal="Computer Physics Communications",
volume="180",
pages="1572-1578",
year="2009"}

\@article{re:alvarez0310,
author="G. Alvarez",
title="Implementation of the SU(2) Hamiltonian Symmetry for the DMRG Algorithm",
journal="Computer Physics Communications",
volume="183",
pages="2226-2232",
year="2012"}


\@article{re:alvarez0311,
author="G. Alvarez and  L. G. G. V. Dias da Silva and
E. Ponce and  E. Dagotto",
title="Time Evolution with the DMRG Algorithm: A Generic Implementation
for Strongly Correlated Electronic Systems",
journal="Phys. Rev. E",
volume="84",
pages="056706",
year="2011"}

\@article{re:alvarez0713,
author="G. Alvarez",
title="Production of minimally entangled typical thermal states
with the Krylov-space approach",
journal="Phys. Rev. B",
volume="87",
pages="245130",
year="2013"}

And also:
\@article{re:alvarez08,
 re:webDmrgPlusPlus,
 Author = {G. Alvarez},
 Title = {DMRG++ Website},
 Publisher = {\\url{http://www.ornl.gov/~gz1/dmrgPlusPlus}} }
</pre>
</small>

## Code Integrity

Hash of the latest commit is also posted at

https://web.ornl.gov/~gz1/hashes.html

## Building and Running DMRG++

### Required Software

- Item GNU C++

- Item (required) The LAPACK library.

 The configure.pl script will ask for the LDFLAGS variable
 to pass to the compiler/linker. If the linux platform was
 chosen the default/suggested LDFLAGS will include -llapack.
 If the osx platform was chosen the default/suggested LDFLAGS will
 include  -framework Accelerate.
 For other platforms the appropriate linker flags must be given.
 More information on LAPACK is at http://netlib.org/lapack/

- Item (required) PsimagLite.

 This is here https://github.com/g1257/PsimagLite/.
 You can do <code>git clone https://github.com/g1257/PsimagLite.git</code> in a separate directory
 outside of the DMRG++ distribution. `configure.pl` will ask you where you put it.

- Item (optional) make or gmake
(only needed to use the Makefile)

- Item (optional) perl
(only needed to run the configure.pl script)

### Building DMRG++ 
<code>
<pre>
 cd src
 perl configure.pl
 (please answer questions regarding depedencies and libraries)
 make
</pre>
</code>

### Running DMRG++
 <code>./dmrg -f input.inp</code>

 Sample input files can be found under <code>TestSuite/inputs/</code>.

<code>configure.pl</code> creates the <code>Makefile</code>
 according to the answers to questions given.
 In the Makefile, LDFLAGS must contain the linker flags to
 link with the LAPACK library. Defaults provided
 automatically by configure.pl should work in most cases.
 If MPI is not selected (serial code) then the compiler will be chosen to be g++.
 Other compilers may work but only the GNU C++ compiler, g++, was tested.
 If MPI is selected then the compiler will be chosen to be mpicxx, which
 is usually a wrapper script around g++ to take care of linking with MPI libraries
 and to include MPI headers. Depending on your MPI installation you might need to
 change the name of this script.

