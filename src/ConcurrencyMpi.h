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

/*! \file ConcurrencyMpi.h
 *
 * Implements the Concurrency.h interface for MPI parallelization mode
 * Implements the Concurrency.h interface for MPI parallelization mode,
 * supports load-balancing
 */
#ifndef CONCURRENCY_MPI_HEADER_H
#define CONCURRENCY_MPI_HEADER_H
#include <mpi.h>
#include "Concurrency.h"
#include "TypeToString.h"

namespace PsimagLite {
	template<typename BogusType=int>
	class ConcurrencyMpi : public Concurrency<BogusType> {

	public:

		typedef MPI_Comm CommType;
		typedef std::pair<CommType,CommType> CommPairType;

		static const CommType  COMM_WORLD;

		ConcurrencyMpi(int argc, char *argv[])
		{
			MPI_Init(&argc,&argv);
		}

		std::string name() const { return "mpi"; }

		~ConcurrencyMpi()
		{
			MPI_Finalize();
			for (size_t i=0;i<garbage_.size();i++)
				delete garbage_[i];
		}

		int nprocs(CommType mpiComm=COMM_WORLD) const
		{ 
			int tmp;
			MPI_Comm_size(mpiComm,&tmp);
			return tmp;
		}

		int rank(CommType mpiComm=COMM_WORLD) const
		{
			int tmp;
			MPI_Comm_rank(mpiComm,&tmp);
			return tmp; 
		}

		bool root(CommType mpiComm=COMM_WORLD) const
		{
			if (rank(mpiComm)==0) return true;
			return false;
		}
		
		CommPairType newCommFromSegments(size_t numberOfSegments,CommType mpiComm=COMM_WORLD)
		{
			size_t procs = nprocs(mpiComm);
			if (procs%numberOfSegments !=0) {
				std::string s("Segment size must be a divisor of nprocs ");
				s += std::string("__FUNCTION__") + __FILE__+" : " + ttos(__LINE__);
				throw std::runtime_error(s.c_str());
			}
			/* Extract the original group handle */ 
			MPI_Group origGroup;
			MPI_Comm_group(mpiComm, &origGroup); 
			
			/* Divide tasks into procs/x distinct groups based upon rank */ 
			size_t segmentSize = size_t(procs/numberOfSegments);
			size_t r = rank(mpiComm);
			std::vector<int> rv;

			getSegmentsDirect(rv,numberOfSegments,segmentSize,r);
			CommType comm1 = getCommFromSegments(origGroup,rv);
			
			getSegmentsAdjuct(rv,numberOfSegments,segmentSize,r);
			CommType comm2 =getCommFromSegments(origGroup,rv);

			return CommPairType(comm1,comm2);
		}

		void reduce(double& v,CommType mpiComm=COMM_WORLD) const
		{
			double w = 0;
			int x = MPI_Reduce(&v,&w,1,MPI_DOUBLE,MPI_SUM,0,mpiComm);
			checkError(x,"MPI_Reduce");
			if (rank(mpiComm)==0) v = w;
		}

		void reduce(std::vector<double>& v,CommType mpiComm=COMM_WORLD) const
		{
			std::vector<double> w(v.size());
			
			int x = MPI_Reduce(&(v[0]),&(w[0]),v.size(),
			                   MPI_DOUBLE,MPI_SUM,0,mpiComm);
			checkError(x,"MPI_Reduce");
			
			if (root(mpiComm)) v = w;
		}

		void reduce(std::vector<std::complex<double> >& v,CommType mpiComm=COMM_WORLD) const
		{
			std::vector<std::complex<double> > w(v.size());
			
			int x = MPI_Reduce(&(v[0]),&(w[0]),2*v.size(),
			                   MPI_DOUBLE,MPI_SUM,0,mpiComm);
			checkError(x,"MPI_Reduce");
			
			if (rank(mpiComm)==0) v = w;
		}

		void reduce(PsimagLite::Matrix<double>& m,CommType mpiComm=COMM_WORLD) const
		{
			PsimagLite::Matrix<double> w(m.n_row(),m.n_col());
			int n = m.n_row()*m.n_col();
			int x = MPI_Reduce(&(m(0,0)),&(w(0,0)),n,MPI_DOUBLE,MPI_SUM,0,mpiComm);
			checkError(x,"MPI_Reduce");
			if (rank(mpiComm)==0) m = w;
		}

		void reduce(PsimagLite::Matrix<std::complex<double> >& m,CommType mpiComm=COMM_WORLD) const
		{
			PsimagLite::Matrix<std::complex<double> > w(m.n_row(),m.n_col());
			int n = 2*m.n_row()*m.n_col();
			int x = MPI_Reduce(&(m(0,0)),&(w(0,0)),n,MPI_DOUBLE,MPI_SUM,0,mpiComm);
			checkError(x,"MPI_Reduce");
			if (rank(mpiComm)==0) m = w;
		}

		void allReduce(std::vector<double>& v,CommType mpiComm=COMM_WORLD) const
		{
			std::vector<double> w(v.size(),0);
			int x = MPI_Allreduce(&(v[0]),&(w[0]),v.size(),MPI_DOUBLE,MPI_SUM,mpiComm);
			checkError(x,"MPI_Allreduce");
			v = w;
		}

		void allReduce(std::vector<std::complex<double> >& v,CommType mpiComm=COMM_WORLD) const
		{
			std::vector<std::complex<double> > w(v.size(),0);
			int x = MPI_Allreduce(&(v[0]),&(w[0]),2*v.size(),MPI_DOUBLE,MPI_SUM,mpiComm);
			checkError(x,"MPI_Allreduce");
			v = w;
		}
		
		void allReduce(PsimagLite::Matrix<double>& m,CommType mpiComm=COMM_WORLD) const
		{
			PsimagLite::Matrix<double> w(m.n_row(),m.n_col());
			int n = m.n_row()*m.n_col();
			int x = MPI_Allreduce(&(m(0,0)),&(w(0,0)),n,MPI_DOUBLE,MPI_SUM,mpiComm);
			checkError(x,"MPI_Allreduce");
			//if (rank(mpiComm)==0) std::cout<<__FUNCTION__<<" has been called\n";
			m = w;
		}

		void allReduce(PsimagLite::Matrix<std::complex<double> >& m,CommType mpiComm=COMM_WORLD) const
		{
			PsimagLite::Matrix<std::complex<double> > w(m.n_row(),m.n_col());
			int n = 2*m.n_row()*m.n_col();
			int x = MPI_Allreduce(&(m(0,0)),&(w(0,0)),n,MPI_DOUBLE,MPI_SUM,mpiComm);
			checkError(x,"MPI_Reduce");
			//if (rank(mpiComm)==0) std::cout<<__FUNCTION__<<" has been called\n";
			m = w;
		}

		void gather(std::vector<PsimagLite::Matrix<double> > &v,CommType mpiComm=COMM_WORLD)
		{
			std::string s = "You hit an unimplemented function.\n";
			s += "Contribute to PsimagLite development and make a difference!\n";
			s += "Implement this function!\n";
			s += ttos(__FUNCTION__) + __FILE__ + " : " + ttos(__LINE__) + "\n"; 
			throw std::runtime_error(s.c_str());
		}

		void gather(std::vector<std::complex<double> >& v,CommType mpiComm=COMM_WORLD)
		{
			std::vector<std::complex<double> > w(v.size(),0);
			int x = MPI_Gather(&(v[0]),2*v.size(), MPI_DOUBLE,&(w[0]),2*w.size(),MPI_DOUBLE,0,mpiComm); 
			checkError(x,"MPI_Gather");
			v = w;
		}

		void gather(std::vector<double>& v,CommType mpiComm=COMM_WORLD)
		{
			std::vector<double> w(v.size(),0);
			assert(v.size()>0);
			assert(w.size()>0);
			assert(v.size()==w.size());
			int x = MPI_Gather(&(v[0]),v.size(), MPI_DOUBLE,&(w[0]),w.size(),MPI_DOUBLE,0,mpiComm);
			checkError(x,"MPI_Gather");
			v = w;
		}


// 		template<typename T>
// 		void broadcast(std::vector<std::vector<T> > &v,CommType mpiComm=COMM_WORLD) 
// 		{ 
// 			for (size_t i=0;i<v.size();i++) MpiBroadcast(&(v[i]),0);
// 		}
// 
// 		template<typename DataType>
// 		void broadcast(std::vector<DataType> &v,CommType mpiComm=COMM_WORLD) 
// 		{ 
// 			for (size_t i=0;i<v.size();i++) MpiBroadcast(&(v[i]),0);
// 		}
// 		
// 		template<typename DataType>
// 		void broadcast(std::vector<DataType*> &v,CommType mpiComm=COMM_WORLD) 
// 		{ 
// 			for (size_t i=0;i<v.size();i++) MpiBroadcast(v[i],0);
// 		}

		void broadcast(bool& b,CommType mpiComm=COMM_WORLD)
		{
			int val = (b) ? 1 : 0;
			int x = MPI_Bcast(&val,1, MPI_INTEGER,0,mpiComm); 
			checkError(x,"MPI_Bcast");
			b = (val==1) ? true : false;
		}

		void broadcast(std::vector<double>& v,CommType mpiComm=COMM_WORLD)
		{
			int x = MPI_Bcast(&(v[0]),v.size(), MPI_DOUBLE,0,mpiComm);
			checkError(x,"MPI_Bcast");
		}

		void broadcast(double& v,CommType mpiComm=COMM_WORLD)
		{
			int x = MPI_Bcast(&v,1, MPI_DOUBLE,0,mpiComm);
			checkError(x,"MPI_Bcast");
		}

		void barrier(CommType mpiComm=COMM_WORLD)
		{
			MPI_Barrier(mpiComm);
		}

	private:
		std::vector<CommType*> garbage_;

		void checkError(int x,const std::string& s) const
		{
			if (x!=MPI_SUCCESS) {
				std::string s2("ConcurrencyMpi::");
				s2 += s + " failed\n";
				throw std::runtime_error(s2.c_str());
			}
		}

		void getSegmentsDirect(std::vector<int>& rv,size_t numberOfSegments,size_t segmentSize,size_t r)
		{
			std::vector<std::vector<int> > ranks;
			size_t thisSegment = 0;
			for (size_t i=0;i<segmentSize;i++) {
				std::vector<int> tmp(numberOfSegments);
				for (size_t j=0;j<numberOfSegments;j++) {
					tmp[j] = j*segmentSize+i;
					if (r==size_t(tmp[j])) thisSegment=i;
				}
				ranks.push_back(tmp);
			}
			rv = ranks[thisSegment];
		}

		void getSegmentsAdjuct(std::vector<int>& rv,size_t numberOfSegments,size_t segmentSize,size_t r)
		{
			
			std::vector<std::vector<int> > ranks;
			size_t thisSegment = 0;
			for (size_t i=0;i<numberOfSegments;i++) {
				std::vector<int> tmp;
				size_t start = i*segmentSize;
				size_t end = (i+1)*segmentSize;
				if (r>=start && r<end) thisSegment = i;
				for (size_t j=start;j<end;j++) tmp.push_back(j);
				ranks.push_back(tmp);
			}
			rv = ranks[thisSegment];
		}

		CommType getCommFromSegments(MPI_Group origGroup,std::vector<int>& rv)
		{
			MPI_Group newGroup;
			int status = MPI_Group_incl(origGroup,rv.size(),&(rv[0]),&newGroup);
			if (status!=MPI_SUCCESS) throw std::runtime_error("getCommFromSegments\n");
			CommType* newComm = new CommType;
			garbage_.push_back(newComm);
			MPI_Comm_create(COMM_WORLD, newGroup, newComm);
			return *newComm;
		}
	}; // class ConcurrencyMpi

	template<typename BogusType>
	const typename ConcurrencyMpi<BogusType>::CommType
	ConcurrencyMpi<BogusType>::COMM_WORLD = MPI_COMM_WORLD;
	
} // namespace Dmrg

/*@}*/
#endif
