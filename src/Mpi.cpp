/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

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

*********************************************************

*/
/** \ingroup PsimagLite */
/*@{*/

/*! \file Mpi.cpp
 *
 */

#include "Mpi.h"

namespace PsimagLite {

namespace MPI {

#ifdef USE_MPI
CommType COMM_WORLD = MPI_COMM_WORLD;

template<>
const MPI_Datatype MpiData<unsigned int long>::Type = MPI_LONG;

template<>
const MPI_Datatype MpiData<unsigned int>::Type = MPI_INTEGER;

template<>
const MPI_Datatype MpiData<double>::Type = MPI_DOUBLE;

template<>
const MPI_Datatype MpiData<int>::Type = MPI_INTEGER;

void checkError(int errorCode,PsimagLite::String caller,CommType comm)
{
	if (errorCode == MPI_SUCCESS)
		return;

	char errorMessage[MPI_MAX_ERROR_STRING];
	int messageLength = 0;
	MPI_Error_string(errorCode,errorMessage,&messageLength);
	std::cerr<<"Error in call to "<<caller<<" ";
	std::cerr<<errorMessage<<"\n";
	MPI_Abort(comm, -1);
}

void init(int* argc, char **argv[])
{
	MPI_Init(argc,argv);
}

void finalize()
{
	MPI_Finalize();
}

SizeType commSize(CommType mpiComm)
{
	int tmp = 1;
	if (mpiComm != 0)
		MPI_Comm_size(mpiComm,&tmp);
	return tmp;
}

SizeType commRank(CommType mpiComm)
{
	int tmp = 0;
	if (mpiComm != 0)
		MPI_Comm_rank(mpiComm,&tmp);
	return tmp;
}

#else

int COMM_WORLD = 0;
int SUM = 0;

void init(int, char **) {}

void finalize() {}

SizeType commSize(CommType)
{
	return 1;
}

SizeType commRank(CommType)
{
	return 0;
}

#endif
} // namespace MPI

} // namespace PsimagLite

