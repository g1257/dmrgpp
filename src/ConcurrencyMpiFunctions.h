// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup PsimagLite */
/*@{*/

/*! \file ConcurrencyMpiFunctions.h
 *
 * 
 */

#ifndef CONCURRENCY_MPIFUNCTIONS_H
#define CONCURRENCY_MPIFUNCTIONS_H

#include <mpi.h>
#include <complex>
#include <iostream>
#include <vector>
#include "Matrix.h" // under psimag
 
namespace PsimagLite {
	
	// MpiBroadcast non-template functions
	inline void MpiBroadcast(std::vector<double> &v,int iproc)
	{
		int x = v.size();
		int y=x;
		MPI_Bcast(&x,1,MPI_INTEGER,iproc,MPI_COMM_WORLD);
		if (x!=y) v.resize(x);
		MPI_Bcast(&(v[0]),x,MPI_DOUBLE,iproc,MPI_COMM_WORLD);
	}

	inline void MpiBroadcast(std::vector<std::complex<double> > &v,int iproc)
	{
		int x = v.size();
		int y=x;
		MPI_Bcast(&x,1,MPI_INTEGER,iproc,MPI_COMM_WORLD);
		if (x!=y) v.resize(x);
		MPI_Bcast(&(v[0]),2*x,MPI_DOUBLE,iproc,MPI_COMM_WORLD);
	}

	inline void MpiBroadcast(std::vector<int> &v,int iproc)
	{
		int x = v.size();
		int y=x;
		//std::cerr<<"before x="<<x<<" y="<<y<<" iproc="<<iproc<<"\n";
		MPI_Bcast(&x,1,MPI_INTEGER,iproc,MPI_COMM_WORLD);
		//std::cerr<<"after x="<<x<<" y="<<y<<"\n";
		//MPI_Barrier(MPI_COMM_WORLD);
		if (x!=y) v.resize(x);
		MPI_Bcast(&(v[0]),x,MPI_INTEGER,iproc,MPI_COMM_WORLD);
		//std::cerr<<"done with vector<int>\n";
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// MpiBroadcast template functions
	template<typename CrsMatrixType>
	void MpiBroadcast(CrsMatrixType *S,int rank)
	{
		//std::cerr<<"rank="<<rank<<" broadcasting S.rowptr.size()="<<S->rowptr.size()<<"\n";
		MpiBroadcast(S->rowptr_,rank);
		//std::cerr<<"rank="<<rank<<" broadcasting S.rcolind.size()="<<S->colind.size()<<"\n";
		MpiBroadcast(S->colind_,rank);
		//std::cerr<<"rank="<<rank<<" broadcasting S.size="<<S->size<<"\n";
		MpiBroadcast(&(S->size_),rank);
		//std::cerr<<"rank="<<rank<<" broadcasting S.values.size()="<<S->values.size()<<"\n";
		MpiBroadcast(S->values_,rank);
	}
	
	template<>
	void MpiBroadcast(double *v,int iproc)
	{
			MPI_Bcast(v,1,MPI_DOUBLE,iproc,MPI_COMM_WORLD);
	}
	
	template<>
	void MpiBroadcast(std::complex<double> *v,int iproc)
	{
			MPI_Bcast(v,2,MPI_DOUBLE,iproc,MPI_COMM_WORLD);
	}
	
	template<>
	void MpiBroadcast(int *v,int iproc)
	{
			MPI_Bcast(v,1,MPI_INTEGER,iproc,MPI_COMM_WORLD);
	}

	template<>
	void MpiBroadcast(psimag::Matrix<double> *m,int iproc)
	{
		int total = m->n_row()*m->n_col();
		MPI_Bcast(&((*m)(0,0)),total,MPI_DOUBLE,iproc,MPI_COMM_WORLD);
	}
	
	template<>
	void MpiBroadcast(std::vector<double> *v,int iproc)
	{
		int x = v->size();
		MPI_Bcast(&((*v)[0]),x,MPI_DOUBLE,iproc,MPI_COMM_WORLD);
	}

	
	

	// MpiRecv non-template functions
	inline void MpiRecv(std::vector<int> &v,int iproc,int tag)
	{
		int x;
		MPI_Status status;
		MPI_Recv(&x,1,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,&status);
		v.resize(x);
		MPI_Recv(&(v[0]),x,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,&status);
	}

	inline void MpiRecv(std::vector<std::complex<double> > &v,int iproc,int tag)
	{
		int x;
		MPI_Status status;
		MPI_Recv(&x,1,MPI_INTEGER,iproc,iproc+tag,MPI_COMM_WORLD,&status);
		v.resize(x);
		MPI_Recv(&(v[0]),2*x,MPI_DOUBLE,iproc,iproc+tag,MPI_COMM_WORLD,&status);
	}
	
	inline void MpiRecv(std::vector<double> &v,int iproc,int tag)
	{
		int x;
		MPI_Status status;
		MPI_Recv(&x,1,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,&status);
		v.resize(x);
		MPI_Recv(&(v[0]),x,MPI_DOUBLE,iproc,tag,MPI_COMM_WORLD,&status);
	}
	
	// MpiRecv template functions
	
	template<typename CrsMatrixType>
	inline void MpiRecv(CrsMatrixType *v,int iproc,int i)
	{
		int tag1=1024,tag2=2048,tag3=3096,tag4=4096;
		MPI_Status status;
		
		MpiRecv(v->rowptr_,iproc,i+tag1);
		MpiRecv(v->colind_,iproc,i+tag2);
		MpiRecv(v->values_,iproc,i+tag3);
		int x;
		MPI_Recv(&x,1,MPI_INTEGER,iproc,i+tag4,MPI_COMM_WORLD,&status);
		v->size_=x;
	}
	
	template<>
	void MpiRecv(psimag::Matrix<double> *m,int iproc,int tag)
	{
		int total = m->n_row()*m->n_col();
		MPI_Status status;
		MPI_Recv(&((*m)(0,0)),total,MPI_DOUBLE,iproc,tag,MPI_COMM_WORLD,&status);
	}
	
	template<>
	void MpiRecv(int *v,int iproc,int tag)
	{
		MPI_Status status;
		MPI_Recv(v,1,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,&status);
	}
	
	template<>
	void MpiRecv(double *v,int iproc,int tag)
	{
		MPI_Status status;
		MPI_Recv(v,1,MPI_DOUBLE,iproc,tag,MPI_COMM_WORLD,&status);
	}

	
	// send non-template
	inline void MpiSend(std::vector<std::complex<double> > &v,int iproc,int tag)
	{
		int x = v.size();
		MPI_Send(&x,1,MPI_INTEGER,0,iproc+tag,MPI_COMM_WORLD);
		MPI_Send(&(v[0]),2*x,MPI_DOUBLE,0,iproc+tag,MPI_COMM_WORLD);
	}
	
	inline void MpiSend(std::vector<int> &v,int iproc,int tag)
	{
		int x = v.size();
		MPI_Send(&x,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD);
		MPI_Send(&(v[0]),x,MPI_INTEGER,0,tag,MPI_COMM_WORLD);
	}

	inline void MpiSend(std::vector<double> &v,int iproc,int tag)
	{
		int x = v.size();
		MPI_Send(&x,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD);
		MPI_Send(&(v[0]),x,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	}


	inline void MpiSend(double *v,int iproc,int tag)
	{
		MPI_Send(v,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	}
	
	inline void MpiSend(int *v,int iproc,int tag)
	{
		MPI_Send(v,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD);
	}


	inline void MpiSend(psimag::Matrix<double> *m,int iproc,int tag)
	{
		int total = m->n_row()*m->n_col();
		MPI_Send(&((*m)(0,0)),total,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
	}

	template<typename CrsMatrixType>
	inline void MpiSend(CrsMatrixType  *v,int iproc,int i)
	{
		int tag1=1024,tag2=2048,tag3=3096,tag4=4096;
		int x;
		
		MpiSend(v->rowptr_,iproc,i+tag1);
		MpiSend(v->colind_,iproc,i+tag2);
		MpiSend(v->values_,iproc,i+tag3);
		x = v->size_;
		MPI_Send(&x,1,MPI_INTEGER,0,i+tag4,MPI_COMM_WORLD);
	}
	
} // namespace Dmrg

/*@}*/
#endif
