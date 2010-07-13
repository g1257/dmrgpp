// BEGIN LICENSE BLOCK
/*
Copyright © 2008 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 1.0.0]
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
/** \ingroup DMRG */
/*@{*/

/*! \file ObserverHelper.h
 *
 *  A class to read and serve precomputed data to the observer
 *
 */
#ifndef PRECOMPUTED_H
#define PRECOMPUTED_H
#include "SparseVector.h"
#include "ProgramGlobals.h"
#include "FermionSign.h"
#include "TimeSerializer.h"
#include "DmrgSerializer.h"
#include "VectorWithOffsets.h" // to include norm
#include "VectorWithOffset.h" // to include norm

namespace Dmrg {
	template<typename IoType,typename MatrixType,typename VectorType_,typename VectorWithOffsetType,typename BasisType>
	class ObserverHelper {
	public:
		typedef VectorType_ VectorType;
		typedef size_t IndexType;
		typedef typename VectorType::value_type FieldType;
		typedef typename BasisType::RealType RealType;
		typedef TimeSerializer<RealType,VectorWithOffsetType> TimeSerializerType;
		typedef DmrgSerializer<RealType,VectorWithOffsetType,MatrixType,BasisType,FermionSign> DmrgSerializerType;
		
		enum {NOTIMEVECTOR=0,USETIMEVECTOR=1};
		
		ObserverHelper(const std::string& filename,size_t nf,bool verbose=true) 
			:	filename_(filename),
				io_(filename),
				dSerializerV_(1,DmrgSerializerType(io_,true)),
				currentPos_(0),
				verbose_(verbose)
		{
			std::cerr<<"Observer will use file: "<<filename<<" for core DMRG data\n";
			init(nf);
		}
		
		ObserverHelper(const std::string& filename,const std::string& timeFilename,size_t nf=0,bool verbose=true) 
			:	filename_(filename),
				io_(filename),
				io2_(timeFilename),
				dSerializerV_(1,DmrgSerializerType(io_,true)),
				timeSerializerV_(nf),
				currentPos_(0),
				verbose_(verbose)
		{
			std::cerr<<"Observer will use file: "<<filename<<" for core DMRG data\n";
			std::cerr<<"Observer will use file: "<<timeFilename<<" for time DMRG data\n";
			init(nf);
			integrityChecks();
		}
		
		void setPointer(size_t pos)
		{
			//std::cerr<<"POS="<<pos<<"\n";
			currentPos_=pos;
		}

		void transform(MatrixType& ret,const MatrixType& O2) const
		{
			return dSerializerV_[currentPos_].transform(ret,O2);
		}
		
		size_t columns() const
		{
			return dSerializerV_[currentPos_].columns();
		}
		
		size_t rows() const
		{
			return dSerializerV_[currentPos_].rows();
		}
		
		const FermionSign& fermionicSign() const
		{
			return dSerializerV_[currentPos_].fermionicSign();
		}
		
		
		const BasisType& basisS() const 
		{
			return dSerializerV_[currentPos_].basisS();
		}

		const BasisType& basisE() const 
		{
			return dSerializerV_[currentPos_].basisE();
		}
		
		const BasisType& basisSE() const 
		{
			return dSerializerV_[currentPos_].basisSE();
		}
		
		size_t direction() const
		{
			return dSerializerV_[currentPos_].direction();
		}

		const VectorWithOffsetType& wavefunction() const
		{
			return dSerializerV_[currentPos_].wavefunction();
		}
		
		RealType time() const
		{
			return timeSerializerV_[currentPos_].time();	
		}
		
		//! This applies more generally (ie. not only to time)
		size_t site() const
		{
			return timeSerializerV_[currentPos_].site();
		}
		
		size_t size() const
		{
			return dSerializerV_.size()-1;
		}
		
		const VectorWithOffsetType& timeVector() const
		{
			if (currentPos_>=timeSerializerV_.size() || 
						 timeSerializerV_[currentPos_].size()==0)
				throw std::runtime_error("timeVector has a problem\n");
			return timeSerializerV_[currentPos_].vector();
		}
		
		template<typename IoType1,typename MatrixType1,typename VectorType1,typename VectorWithOffsetType1,typename BasisType1>
		friend std::ostream& operator<<(std::ostream& os,
			ObserverHelper<IoType1,MatrixType1,VectorType1,VectorWithOffsetType1,BasisType1>& precomp);

	private:
		void init(size_t nf)
		{
			rewind(true);
			//std::vector<size_t> el0; // not really needed, but needs to read to keep in sync
			//getElectronsOneSite(el0);
			//for (size_t i=0;i<nf-1;i++) {
			dSerializerV_.clear();
			while(true) {
				if (nf>0 && dSerializerV_.size()>=nf) break;
				if (verbose_) std::cerr<<"ObserverHelper "<<dSerializerV_.size()<<"\n";
				try {
					DmrgSerializerType dSerializer(io_);
					dSerializerV_.push_back(dSerializer);
				} catch (std::exception& e)
				{
					if (dSerializerV_.size()==0) {
						std::cerr<<e.what()<<" rethrowing...\n";
						throw e;
					}
					break;
				}
			}
			
			FieldType dummy = 0;
			initTimeVectors(dSerializerV_.size(),dummy);
			// Line below might cause trouble under gcc v3
			//if (verbose_) std::cerr<<(*this);	
		}
		
		void integrityChecks()
		{
			if (dSerializerV_.size()!=timeSerializerV_.size()) throw std::runtime_error("Error 1\n");
			if (dSerializerV_.size()==0) return;
			for (size_t x=0;x<dSerializerV_.size()-1;x++) {
				if (dSerializerV_[x].basisSE().size()==0) continue;
				if (dSerializerV_[x].basisSE().size()!=timeSerializerV_[x].size()) throw std::runtime_error("Error 2\n");
			}
			
		}
		
		void initTimeVectors(size_t nf,RealType dummy)
		{
			
		}
		
		void initTimeVectors(size_t nf,std::complex<RealType> dummy)
		{
			if (nf!=timeSerializerV_.size()) timeSerializerV_.resize(nf);
			for (size_t i=0;i<timeSerializerV_.size();i++) { // up to i<nf-1 FIXME
				TimeSerializerType ts(io2_);
				timeSerializerV_[i] = ts;
				//std::cerr<<"time vector "<<i<<" has size "<<psiTimeVector_[i].size()<<"\n";
				//RealType tmp = std::norm(psiTimeVector_[i]);
				//std::cerr<<"Mod="<<tmp<<"\n";
				//std::cerr<<"----------------------------------\n";
			}
		}
		
		void rewind(bool doIt=false) 
		{
			if (doIt) {
				io_.rewind();
				io2_.rewind();
			}
		}
		
		// Not needed, but if you remove this, also remove in DmrgSolver the corresponding
		// printing of the first basis to keep everythign in sync
		/*void getElectronsOneSite(std::vector<size_t>& el0)
		{
			BasisType b(io_,"one site");
			if (b.block().size()!=1) throw std::runtime_error("getElectronsOneSite\n");
			el0 = b.electronsVector();
		}*/
		

		void getTransform(MatrixType& transform,int ns)
		{
			io_.readMatrix(transform,"#TRANSFORM_sites",ns);
			rewind();
		}

		void getWaveFunction(VectorType& wavefunction,size_t ns)
		{
			VectorWithOffsetType tmpV;
			tmpV.load(io_,"#WAVEFUNCTION_sites=",ns);
			tmpV.toSparse(wavefunction);
			rewind();
		}
		
		void getDirection(int& x,int ns)
		{
			io_.readline(x,"#DIRECTION=",ns);
			rewind();
		}
		
		void getTimeVector(VectorWithOffsetType& timeVector,size_t ns)
		{
			timeVector.load(io2_,"targetVector0",ns);
			rewind();
		}

		std::string filename_; 
		typename IoType::In io_;
		typename IoType::In io2_;
		std::vector<DmrgSerializerType> dSerializerV_;
		std::vector<TimeSerializerType> timeSerializerV_;
		size_t currentPos_;
		bool verbose_;
	};  //ObserverHelper

	template<typename IoType1,typename MatrixType1,typename VectorType1,typename VectorWithOffsetType1,typename BasisType1>
	std::ostream& operator<<(std::ostream& os,
		ObserverHelper<IoType1,MatrixType1,VectorType1,VectorWithOffsetType1,BasisType1>& p)
	{
		for (size_t i=0;i<p.SpermutationInverse_.size();i++) {
			os<<"i="<<i<<"\n";
			os<<"\tS.size="<<p.SpermutationInverse_[i].size()<<" "<<p.Spermutation_[i].size()<<"\n";
			os<<"\tSE.size="<<p.SEpermutationInverse_[i].size()<<" "<<p.SEpermutation_[i].size()<<"\n";
			os<<"\tElectrons.size="<<p.electrons_[i].size()<<"\n";
			os<<"\tTransform="<<p.transform_[i].n_row()<<"x"<<p.transform_[i].n_col()<<"\n";
			os<<"\tWF.size="<<p.wavefunction_[i].size()<<"\n";
		}
		return os;
	}
} // namespace Dmrg

/*@}*/
#endif
