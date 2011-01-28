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
 */
#ifndef CONCURRENCY_MPI_HEADER_H
#define CONCURRENCY_MPI_HEADER_H
#include "ConcurrencyMpiFunctions.h"
#include "Concurrency.h"

namespace PsimagLite {
	//! Implements the Concurrency.h interface for MPI parallelization mode, supports load-balancing
	template<typename FieldType>
	class ConcurrencyMpi : public Concurrency<FieldType> {
	public:
		ConcurrencyMpi(int argc, char *argv[])
		{
			MPI_Init(&argc,&argv);
			step_= -1;
			total_=0;
			rank_ = rank();
			nprocs_=nprocs();
			
		}

		std::string name() const { return "mpi"; }

		~ConcurrencyMpi()
		{
			MPI_Finalize();
		}

		int nprocs(MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{ 
			int tmp;
			MPI_Comm_size(mpiComm,&tmp);
			return tmp;
		}

		int rank(MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{
			int tmp;
			MPI_Comm_rank(mpiComm,&tmp);
			return tmp; 
		}

		bool root() 
		{
			if (rank_==0) return true;
			return false;
		}

		void loopCreate(size_t total,std::vector<size_t> const &weights,MPI_Comm mpiComm=MPI_COMM_WORLD)
		{
			nprocs_=nprocs(mpiComm);
			rank_=rank(mpiComm);
			total_ = total;
			step_=0;

			// distribute the load among the processors
			std::vector<int> loads(nprocs_,0);
			indicesOfThisProc_.resize(nprocs_);
			int r;

			for (r=0;r<nprocs_;r++) indicesOfThisProc_[r].clear();

			assigned_=false;
			for (size_t i=0;i<total_;i++) {
				r = findLowestLoad(loads);
				indicesOfThisProc_[r].push_back(i);
				loads[r] += weights[i];
				if (r==rank_) assigned_=true;
			}
			// set myIndices_
			myIndices_=indicesOfThisProc_[rank_];
			MPI_Barrier(mpiComm);
		}

		void loopCreate(size_t total,MPI_Comm mpiComm=MPI_COMM_WORLD)
		{
			std::vector<size_t> weights(total,1);
			loopCreate(total,weights,mpiComm);
		}

		bool loop(size_t &i)
		{
			if (!assigned_) return false;
			
			if (step_<0 || total_==0) throw std::runtime_error("ConcurrencySerial::loop() loopCreate() must be called before.\n"); 
			 //i = rank_ + step_*nprocs_;
			i=myIndices_[step_];
			
			if (i>=total_ || step_>=int(myIndices_.size())) {
				step_= -1; 
				return false;
			}
			step_++;
			return true;
		}

		void loopReset()
		{
			step_=0;
		}

		void gather(std::vector<std::vector<std::complex<double> > > &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{
			int i,x;
			std::vector<std::complex<double> > tmpVec;
			MPI_Status status;
			int tag=999;

			if (!assigned_) return;

			if (total_>v.size()  || myIndices_.size()<=0) { 
				std::cerr<<"total_="<<total_<<" v.size()="<<v.size()<<" myindices.size="<<myIndices_.size()<<" line="<<__LINE__<<"\n";
				throw std::runtime_error("ConcurrencyMpi::gather() loopCreate() must be called before.\n"); 
			}
			if (rank_>0) {
				for (step_=0;step_<myIndices_.size();step_++) {
					//i = rank_ + step_*nprocs_;
					i=myIndices_[step_];
					if (i>=total_) break;
					//std::cerr<<"sent rank="<<rank_<<" step="<<step_<<" tag="<<i<<"\n";
					x=v[i].size();
					MPI_Send(&x,1,MPI_INTEGER,0,i,mpiComm);
					//std::cerr<<"sent rank="<<rank_<<" step="<<step_<<" tag="<<tag<<"\n";
					MPI_Send(&(v[i][0]),2*x,MPI_DOUBLE,0,i,mpiComm);
				
				}
			} else {
				for (int r=1;r<nprocs_;r++) {
					for (step_=0;step_<indicesOfThisProc_[r].size();step_++) {
						//i = r + step_*nprocs_;
						i = indicesOfThisProc_[r][step_];
						if (i>=total_) continue;
						
						MPI_Recv(&x,1,MPI_INTEGER,r,i,mpiComm,&status);
						tmpVec.resize(x);
						MPI_Recv(&(tmpVec[0]),2*x,MPI_DOUBLE,r,i,mpiComm,&status);
						v[i]=tmpVec;
					}
				}
			}
			step_= -1;
		}

		void gather(std::vector<std::vector<double> > &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{
			int x;
			size_t i;
			std::vector<double> tmpVec;
			MPI_Status status;

			if (!assigned_) return;

			if (total_>v.size() ||  myIndices_.size()<=0) { 
				std::cerr<<"total_="<<total_<<" v.size()="<<v.size()<<" myindices.size="<<myIndices_.size()<<" line="<<__LINE__<<"\n";
				throw std::runtime_error("ConcurrencyMpi::gather() loopCreate() must be called before.\n"); 
			}
			if (rank_>0) {
				for (step_=0;step_<int(myIndices_.size());step_++) {
					//i = rank_ + step_*nprocs_;
					i=myIndices_[step_];
					if (i>=total_) break;
					//std::cerr<<"sent rank="<<rank_<<" step="<<step_<<" tag="<<i<<"\n";
					x=v[i].size();
					MPI_Send(&x,1,MPI_INTEGER,0,i,mpiComm);
					//std::cerr<<"sent rank="<<rank_<<" step="<<step_<<" tag="<<tag<<"\n";
					MPI_Send(&(v[i][0]),x,MPI_DOUBLE,0,i,mpiComm);
				}
			} else {
				//std::cerr<<"rank_="<<rank_<<"nprocs="<<nprocs_<<"mpirun="<<mpirun<<"\n";
				for (int r=1;r<nprocs_;r++) {
					for (step_=0;step_<int(indicesOfThisProc_[r].size());step_++) {
						//i = r + step_*nprocs_;
						i=indicesOfThisProc_[r][step_];
						if (i>=total_) continue;
						
						MPI_Recv(&x,1,MPI_INTEGER,r,i,mpiComm,&status);
						tmpVec.resize(x);
						MPI_Recv(&(tmpVec[0]),x,MPI_DOUBLE,r,i,mpiComm,&status);
						v[i]=tmpVec;
					}
				}
			}
			step_= -1;
		}
		
		template<typename T>
		void gather(std::vector<T> &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{
			size_t i;

			if (!assigned_) return;

			if (total_>v.size() ||  myIndices_.size()<=0) { 
				std::cerr<<"total_="<<total_<<" v.size()="<<v.size()<<" myIndices_.size()="<<myIndices_.size()<<"\n";
				throw std::runtime_error("ConcurrencyMpi::broadcast() loopCreate() must be called before.\n"); 
			}

			if (rank_>0) {
				for (step_=0;step_<int(myIndices_.size());step_++) {
					//i = rank_ + step_*nprocs_;
					i=myIndices_[step_];
					if (i>=total_) break;
					MpiSend(&(v[i]),rank_,i);
				
				}
			} else {
				//std::cerr<<"rank_="<<rank_<<"nprocs="<<nprocs_<<"mpirun="<<mpirun<<"\n";
				for (int r=1;r<nprocs_;r++) {
					for (step_=0;step_<int(indicesOfThisProc_[r].size());step_++) {
						i=indicesOfThisProc_[r][step_];
						//i = r + step_*nprocs_;
						if (i>=total_) continue;
						MpiRecv(&(v[i]),r,i);
					}
				}
			}
			step_= -1;
		}
		
		template<typename T>
		void gather(std::vector<PsimagLite::Matrix<T> > &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{
			size_t i;

			if (!assigned_) return;

			if (total_>v.size() ||  myIndices_.size()<=0) { 
				std::cerr<<"total_="<<total_<<" v.size()="<<v.size()<<" myIndices_.size()="<<myIndices_.size()<<"\n";
				throw std::runtime_error("ConcurrencyMpi::broadcast() loopCreate() must be called before.\n"); 
			}

			if (rank_>0) {
				for (step_=0;step_<int(myIndices_.size());step_++) {
					//i = rank_ + step_*nprocs_;
					i=myIndices_[step_];
					if (i>=total_) break;
					int nrow = v[i].n_row();
					int ncol = v[i].n_col();
					MpiSend(&nrow,rank_,i);
					MpiSend(&ncol,rank_,i);
					MpiSend(&(v[i]),rank_,i);
				
				}
			} else {
				//std::cerr<<"rank_="<<rank_<<"nprocs="<<nprocs_<<"mpirun="<<mpirun<<"\n";
				for (int r=1;r<nprocs_;r++) {
					for (step_=0;step_<int(indicesOfThisProc_[r].size());step_++) {
						i=indicesOfThisProc_[r][step_];
						//i = r + step_*nprocs_;
						if (i>=total_) continue;
						int nrow,ncol;
						MpiRecv(&nrow,r,i);
						MpiRecv(&ncol,r,i);
						//std::cerr<<"r="<<r<<" nrow="<<nrow<<" ncol="<<ncol<<" i="<<i<<"\n";
						v[i].resize(nrow,ncol);
						MpiRecv(&(v[i]),r,i);
					}
				}
			}
			step_= -1;
		}
		
		
		template<typename T>
		void gather(std::vector<T*> &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{
			size_t i;

			if (!assigned_) return;

			if (total_>v.size() ||  myIndices_.size()<=0) { 
				std::cerr<<"total_="<<total_<<" v.size()="<<v.size()<<" myIndices_.size()="<<myIndices_.size()<<"\n";
				throw std::runtime_error("ConcurrencyMpi::broadcast() loopCreate() must be called before.\n"); 
			}

			if (rank_>0) {
				for (step_=0;step_<int(myIndices_.size());step_++) {
					//i = rank_ + step_*nprocs_;
					i=myIndices_[step_];
					if (i>=total_) break;
					MpiSend(v[i],rank_,i);
				
				}
			} else {
				//std::cerr<<"rank_="<<rank_<<"nprocs="<<nprocs_<<"mpirun="<<mpirun<<"\n";
				for (int r=1;r<nprocs_;r++) {
					for (step_=0;step_<int(indicesOfThisProc_[r].size());step_++) {
						i=indicesOfThisProc_[r][step_];
						//i = r + step_*nprocs_;
						if (i>=total_) continue;
						MpiRecv(v[i],r,i);
					}
				}
			}
			step_= -1;
		}

		template<typename T>
		void broadcast(std::vector<std::vector<T> > &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{ 
			for (size_t i=0;i<v.size();i++) MpiBroadcast(&(v[i]),0);
		}

		template<typename DataType>
		void broadcast(std::vector<DataType> &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{ 
			for (size_t i=0;i<v.size();i++) MpiBroadcast(&(v[i]),0);
		}
		
		template<typename DataType>
		void broadcast(std::vector<DataType*> &v,MPI_Comm mpiComm=MPI_COMM_WORLD) 
		{ 
			for (size_t i=0;i<v.size();i++) MpiBroadcast(v[i],0);
		}

		void barrier(MPI_Comm mpiComm=MPI_COMM_WORLD)
		{
			MPI_Barrier(mpiComm);
		
		}

	private:
		int nprocs_; // total number of processors
		std::vector<int> myIndices_; // indices assigned to this processor
		int rank_; // rank of this processor
		int step_; // step within this processor
		size_t total_; // total number of indices
		std::vector<std::vector<int> > indicesOfThisProc_; // given rank and step it maps the index
		bool assigned_;

		void MpiGather(std::vector<double> &vrec,double vsend,int iproc,MPI_Comm mpiComm)
		{
			MPI_Gather(&vsend,1,MPI_DOUBLE,&(vrec[0]),1,MPI_DOUBLE,iproc,mpiComm);
		}

		void MpiGather(std::vector<FieldType> &vrec,FieldType &vsend,int iproc,MPI_Comm mpiComm)
		{
			MPI_Gather(&vsend,2,MPI_DOUBLE,&(vrec[0]),2,MPI_DOUBLE,iproc,mpiComm);
		}

		void MpiGather(std::vector<std::vector<FieldType> > &vrec,std::vector<FieldType> &vsend,int iproc,MPI_Comm mpiComm)
		{
			int x = vsend.size();
			MPI_Gather(&(vsend[0]),2*x,MPI_DOUBLE,&(vrec[0][0]),2*x,MPI_DOUBLE,iproc,mpiComm);
		}

		int findLowestLoad(std::vector<int> const &loads)
		{
			int x= 1000000;
			int ret=0;
			for (size_t i=0;i<loads.size();i++) {
				if (loads[i]<x) {
					x=loads[i];
					ret =i;
				}
			}
			return ret;
		}
	}; // class ConcurrencyMpi

} // namespace Dmrg

/*@}*/
#endif
