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
		enum {NOTIMEVECTOR=0,USETIMEVECTOR=1};
		
		ObserverHelper(const std::string& filename,size_t nf,bool verbose=true) 
			:	filename_(filename),
				io_(filename),
				bogusBasis_("Bogus"),
				fermionSigns_(nf,bogusBasis_.electronsVector()),
				basisS_(nf,bogusBasis_),
				basisE_(nf,bogusBasis_),
				basisSE_(nf,bogusBasis_),
				transform_(nf),
				directions_(nf),
				wavefunction_(nf),
				psiTimeVector_(nf),
				currentPos_(0),
				verbose_(verbose),
				nf_(nf),
				stepTimes_(0)
		{
			init(nf);
		}
		
		//! stepTimes must always be equal to nf
		ObserverHelper(const std::string& filename,const std::string& timeFilename,size_t nf,size_t stepTimes,bool verbose=true) 
			:	filename_(filename),
				io_(filename),
				io2_(timeFilename),
				bogusBasis_("Bogus"),
				fermionSigns_(nf,bogusBasis_.electronsVector()),
				basisS_(nf,bogusBasis_),
				basisE_(nf,bogusBasis_),
				basisSE_(nf,bogusBasis_),
				transform_(nf),
				directions_(nf),
				wavefunction_(nf),
				psiTimeVector_(nf),
				currentPos_(0),
				verbose_(verbose),
				nf_(nf),
				stepTimes_(stepTimes)
		{
			init(nf);
			integrityChecks();
		}
		
		void setPointer(size_t pos)
		{
			//std::cerr<<"POS="<<pos<<"\n";
			currentPos_=pos;
		}

		const MatrixType& transform() const
		{
			return transform_[currentPos_];
		}

		const FermionSign& fermionicSign() const
		{
			return fermionSigns_[currentPos_];
		}
		
		
		const BasisType& basisS() const 
		{
			return basisS_[currentPos_];
		}

		const BasisType& basisE() const 
		{
			return basisE_[currentPos_];
		}
		
		const BasisType& basisSE() const 
		{
			return basisSE_[currentPos_];
		}

		void transform(MatrixType& ret,MatrixType& O)
		{
			//typedef typename MatrixType::value_type FieldType;
			int nBig = O.n_row();
			int nSmall = transform_[currentPos_].n_col();
			MatrixType fmTmp(nSmall,nBig);
			FieldType alpha=1.0,beta=0.0;
		
			psimag::BLAS::GEMM('N','N',nBig,nSmall,nBig,alpha,
					   &(O(0,0)),nBig,&(transform_[currentPos_](0,0)),nBig,beta,&(fmTmp(0,0)),nBig);
			
			psimag::BLAS::GEMM('C','N',nSmall,nSmall,nBig,alpha,
					   &(transform_[currentPos_](0,0)),nBig,&(fmTmp(0,0)),nBig,beta,&(ret(0,0)),nSmall);
			
		}
		
		size_t direction() const
		{
			return directions_[currentPos_];
		}

		const VectorType& wavefunction() const
		{
			return wavefunction_[currentPos_];
		}
		
		const VectorWithOffsetType& timeVector() const
		{
			if (currentPos_>=psiTimeVector_.size() || 
						 psiTimeVector_[currentPos_].size()==0)
				throw std::runtime_error("timeVector has a problem\n");
			return psiTimeVector_[currentPos_]; //-nf_+1+stepTimes_];	
		}
		
		template<typename IoType1,typename MatrixType1,typename VectorType1,typename VectorWithOffsetType1,typename BasisType1>
		friend std::ostream& operator<<(std::ostream& os,
			ObserverHelper<IoType1,MatrixType1,VectorType1,VectorWithOffsetType1,BasisType1>& precomp);

	private:
		void init(size_t nf)
		{
			rewind(true);
			std::vector<size_t> el0; // not really needed, but needs to read to keep in sync
			getElectronsOneSite(el0);
			for (size_t i=0;i<nf-1;i++) {
				if (verbose_) std::cerr<<"ObserverHelper "<<i<<" out of "<<(nf-1)<<"\n";
				size_t j = 0; // = i;
				
				/*getPermutation(Spermutation_[i],SpermutationInverse_[i],
						"#pSprime.permutationInverse_sites",j);
				getPermutation(SEpermutation_[i],SEpermutationInverse_[i],
						"#pSE.permutationInverse_sites=",j);
				*/
				fermionSigns_[i].load(io_);
				basisS_[i].load(io_);
				basisE_[i].load(io_);
				basisSE_[i].load(io_);
				getWaveFunction(wavefunction_[i],j);
				//getElectrons(electrons_[i],j);
				getTransform(transform_[i],j);
				int x = 0;
				getDirection(x,j);
				if (x<0) throw std::runtime_error("OBserverHelper:: direction must be non-negative\n");
				directions_[i] = x;
			}
			
			FieldType dummy = 0;
			initTimeVectors(dummy);
			// Line below might cause trouble under gcc v3
			//if (verbose_) std::cerr<<(*this);	
		}
		
		void integrityChecks()
		{
			if (basisSE_.size()!=psiTimeVector_.size()) throw std::runtime_error("Error 1\n");
			for (size_t x=0;x<basisSE_.size();x++) {
				if (basisSE_[x].size()==0) continue;
				if (basisSE_[x].size()!=psiTimeVector_[x].size()) throw std::runtime_error("Error 2\n");
			}
			
		}
		
		void initTimeVectors(size_t nf,RealType dummy)
		{
			
		}
		
		void initTimeVectors(std::complex<RealType> dummy)
		{
			std::cerr<<"steptimes = "<<stepTimes_<<"\n";
			for (size_t i=0;i<stepTimes_;i++) { // up to i<nf-1 FIXME
				size_t j = 0;
				getTimeVector(psiTimeVector_[i],j);
				std::cerr<<"time vector "<<i<<" has size "<<psiTimeVector_[i].size()<<"\n";
				RealType tmp = std::norm(psiTimeVector_[i]);
				std::cerr<<"Mod="<<tmp<<"\n";
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
		void getElectronsOneSite(std::vector<size_t>& el0)
		{
			BasisType b("one site");
			b.load(io_);
			if (b.block().size()!=1) throw std::runtime_error("getElectronsOneSite\n");
			el0 = b.electronsVector();
		}
		

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
		BasisType bogusBasis_;
		std::vector<FermionSign> fermionSigns_;
		std::vector<BasisType> basisS_,basisE_,basisSE_;
		std::vector<MatrixType> transform_;
		std::vector<size_t> directions_;
		std::vector<VectorType> wavefunction_;
		std::vector<VectorWithOffsetType> psiTimeVector_;
		size_t currentPos_;
		bool verbose_;
		size_t nf_;
		size_t stepTimes_;
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
