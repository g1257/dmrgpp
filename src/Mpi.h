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

/*! \file Mpi.h
 *
 */
#ifndef MPI_HEADER_H
#define MPI_HEADER_H
#include <stdexcept>
#include "Vector.h"
#ifdef USE_MPI
#include <mpi.h>
#include <loki/TypeTraits.h>
#endif

namespace PsimagLite {

namespace MPI {

#ifdef USE_MPI
typedef MPI_Comm CommType;
CommType COMM_WORLD = MPI_COMM_WORLD;

template<typename T>
struct MpiData
{
	static const MPI_Datatype Type;
};

template<>
const MPI_Datatype MpiData<unsigned int long>::Type = MPI_LONG;

template<>
const MPI_Datatype MpiData<unsigned int>::Type = MPI_INTEGER;

template<>
const MPI_Datatype MpiData<double>::Type = MPI_DOUBLE;

void checkError(int errorCode,const PsimagLite::String& caller,CommType comm = COMM_WORLD)
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

void init(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
}

void finalize()
{
	MPI_Finalize();
}

SizeType commSize(CommType mpiComm)
{
	int tmp = 0;
	MPI_Comm_size(mpiComm,&tmp);
	return tmp;
}

SizeType commRank(CommType mpiComm)
{
	int tmp = 0;
	MPI_Comm_rank(mpiComm,&tmp);
	return tmp;
}

template<typename NumericType>
typename EnableIf<Loki::TypeTraits<NumericType>::isArith,
void>::Type bcast(NumericType& v,int root = 0, CommType mpiComm = COMM_WORLD)
{
	MPI_Datatype datatype = MpiData<NumericType>::Type;
	int errorCode = MPI_Bcast(&v,1,datatype,root,mpiComm);
	checkError(errorCode,"MPI_Bcast",mpiComm);
}

template<typename SomeVectorType>
typename EnableIf<IsVectorLike<SomeVectorType>::True &
Loki::TypeTraits<typename SomeVectorType::value_type>::IsArith,
void>::Type allGather(SomeVectorType& v,CommType mpiComm = COMM_WORLD)
{
	SomeVectorType recvbuf = v;
	MPI_Datatype datatype = MpiData<typename SomeVectorType::value_type>::Type;
	int errorCode = MPI_Allgather(&(v[0]),v.size(),datatype,&(recvbuf[0]),v.size(),datatype,mpiComm);
	checkError(errorCode,"MPI_Allgather",mpiComm);

	v = recvbuf;
}

template<typename SomeVectorType>
typename EnableIf<IsVectorLike<SomeVectorType>::True &
IsVectorLike<typename SomeVectorType::value_type>::True,
void>::Type allGather(SomeVectorType& v,CommType mpiComm = COMM_WORLD)
{
	typedef typename SomeVectorType::value_type DataType;
	for (SizeType i=0;i<v.size();i++) {
		DataType& vv = v[i];
		DataType recvbuf = vv;
		MPI_Datatype datatype = MpiData<typename DataType::value_type>::Type;
		int errorCode = MPI_Allgather(&(vv[0]),vv.size(),datatype,&(recvbuf[0]),vv.size(),datatype,mpiComm);
		checkError(errorCode,"MPI_Allgather",mpiComm);
		vv = recvbuf;
	}
}

template<typename SomeVectorType>
typename EnableIf<IsVectorLike<SomeVectorType>::True &
Loki::TypeTraits<typename SomeVectorType::value_type>::IsArith,
void>::Type gather(SomeVectorType& v,int root = 0, CommType mpiComm = COMM_WORLD)
{
	SomeVectorType recvbuf(v.size());
	MPI_Datatype datatype = MpiData<typename SomeVectorType::value_type>::Type;
	int errorCode = MPI_Gather(&(v[0]),v.size(),datatype,&(recvbuf[0]),v.size(),datatype,root,mpiComm);
	checkError(errorCode,"MPI_Gather",mpiComm);

	if (commRank(mpiComm) == root)
		v = recvbuf;
}

template<typename NumericType>
typename EnableIf<Loki::TypeTraits<NumericType>::isArith,
void>::Type gather(NumericType& v,int root = 0, CommType mpiComm = COMM_WORLD)
{
	int recvbuf = 0;
	MPI_Datatype datatype = MpiData<NumericType>::Type;
	int errorCode = MPI_Gather(&v,1,datatype,&recvbuf,1,datatype,root,mpiComm);
	checkError(errorCode,"MPI_Gather",mpiComm);

	if (commRank(mpiComm) == static_cast<SizeType>(root))
		v = recvbuf;
}

template<typename SomeVectorType>
typename EnableIf<IsVectorLike<SomeVectorType>::True &
Loki::TypeTraits<typename SomeVectorType::value_type>::isIntegral,
void>::Type allReduce(SomeVectorType& v,MPI_Op op = MPI_SUM, CommType mpiComm = COMM_WORLD)
{
	SomeVectorType recvbuf = v;
	MPI_Datatype datatype = MpiData<typename SomeVectorType::value_type>::Type;
	int errorCode = MPI_Allreduce(&(v[0]),&(recvbuf[0]),v.size(),datatype,op,mpiComm);
	checkError(errorCode,"MPI_Allreduce",mpiComm);
	v = recvbuf;
}

void allReduce(SizeType& v,MPI_Op op = MPI_SUM, CommType mpiComm = COMM_WORLD)
{
	int result = 0;
	int errorCode = MPI_Allreduce(&v,&result,1,MPI_INTEGER,op,mpiComm);
	checkError(errorCode,"MPI_Allreduce",mpiComm);
	v = result;
}

#else
typedef int CommType;
int COMM_WORLD = 0;

void init(int argc, char *argv[]) {}

void finalize() {}

SizeType commSize(CommType mpiComm)
{
	return 1;
}

SizeType commRank(CommType mpiComm)
{
	return 0;
}

template<typename T>
void bcast(T &t)
{}

template<typename T>
void allGather(T &t)
{}

template<typename T>
void gather(T &t)
{}

template<typename T>
void allReduce(T &t)
{}

#endif
} // namespace MPI

template<typename InstanceType>
class Mpi {
public:

	Mpi(SizeType npthreads=1,MPI::CommType comm = MPI::COMM_WORLD)
	    : comm_(comm)
	{
	}

	void loopCreate(SizeType total,
	                InstanceType& pfh)
	{
		SizeType procs = threads();
		SizeType rank = MPI::commRank(comm_);
		SizeType block = static_cast<SizeType>(total/procs);
		if (total % procs !=0) block++;
		pfh.thread_function_(rank,block,total,0);
	}

	String name() const { return "mpi"; }

	SizeType threads() const
	{
		return MPI::commSize(comm_);
	}

private:

	MPI::CommType comm_;

}; // Mpi

} // namespace PsimagLite 

/*@}*/	
#endif
