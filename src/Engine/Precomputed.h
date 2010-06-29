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

/*! \file Precomputed.h
 *
 *  A class to read and serve precomputed data to the observer
 *
 */
#ifndef PRECOMPUTED_H
#define PRECOMPUTED_H
#include "SparseVector.h"
#include "ProgramLimits.h"

namespace Dmrg {
	template<typename RealType,typename FieldType,typename IoType,
 		typename MatrixType,template<typename> class VectorTemplate>
	class Precomputed {
	public:
		typedef size_t IndexType;
		typedef VectorTemplate<FieldType> VectorType;
		enum {NOTIMEVECTOR=0,USETIMEVECTOR=1};
		
		
		
		Precomputed(const std::string& filename,size_t nf,size_t stepTimes,bool verbose=true) 
			:	filename_(filename),
				io_(filename),
				SpermutationInverse_(nf),Spermutation_(nf),
				SEpermutationInverse_(nf),SEpermutation_(nf),
				electrons_(nf),
				transform_(nf),
				wavefunction_(nf),
				psiTimeVector_(nf),
				currentPos_(0),
				verbose_(verbose),
				nf_(nf),
				stepTimes_(stepTimes)
		{
			rewind(true);
			for (size_t i=0;i<nf-1;i++) {
				if (verbose_) std::cerr<<"Precomputed "<<i<<" out of "<<(nf-1)<<"\n";
				size_t j = 0; // = i;
				
				getPermutation(Spermutation_[i],SpermutationInverse_[i],
						"#pSprime.permutationInverse_sites",j);
				getPermutation(SEpermutation_[i],SEpermutationInverse_[i],
						"#pSE.permutationInverse_sites=",j);
				getWaveFunction(wavefunction_[i],j);
				getElectrons(electrons_[i],j);
				getTransform(transform_[i],j);
				
			}
			FieldType dummy = 0;
			initTimeVectors(dummy);
			// Line below might cause trouble under gcc v3
			//if (verbose_) std::cerr<<(*this);
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

		const std::vector<IndexType>& electrons() const
		{
			return electrons_[currentPos_];
		}

		IndexType electrons(size_t i) const
		{
			return electrons_[currentPos_][i];
		}

		IndexType SEpermutation(size_t i) const
		{
			return SEpermutation_[currentPos_][i];
		}

		IndexType SEpermutation() const
		{
			return SEpermutation_[currentPos_].size();
		}

		IndexType SEpermutationInverse(size_t i) const
		{
			return SEpermutationInverse_[currentPos_][i];
		}

		IndexType Spermutation(size_t i) const
		{
			return 	Spermutation_[currentPos_][i];
		}

		const std::vector<IndexType>&   Spermutation() const
		{
			return 	Spermutation_[currentPos_];
		}

		const std::vector<IndexType>&  SpermutationInverse() const
		{
			return SpermutationInverse_[currentPos_];
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

		const VectorType& wavefunction() const
		{
			return wavefunction_[currentPos_];
		}
		
		const VectorType& timeVector() const
		{
			if (currentPos_>=psiTimeVector_.size() || 
						 psiTimeVector_[currentPos_].size()==0)
				throw std::runtime_error("timeVector has a problem\n");
			return psiTimeVector_[currentPos_]; //-nf_+1+stepTimes_];	
		}
		
		template<typename RealType1,
  			typename FieldType1,typename IoType1,typename MatrixType1,template<typename> class VectorTemplate1>
		friend std::ostream& operator<<(std::ostream& os,
			Precomputed<RealType1,FieldType1,IoType1,MatrixType1,VectorTemplate1>& precomp);

	private:
		
		void initTimeVectors(size_t nf,RealType dummy)
		{
			
		}
		
		void initTimeVectors(std::complex<RealType> dummy)
		{
			for (size_t i=0;i<stepTimes_;i++) {
				size_t j = 0;
				getTimeVector(psiTimeVector_[i],j);
				std::cerr<<psiTimeVector_[i];
				std::cerr<<"----------------------------------\n";
			}
		}
		
		void rewind(bool doIt=false) 
		{
			if (doIt) io_.rewind();
		}

		void  getPermutation(std::vector<IndexType>& pS,std::vector<IndexType>& pSi,
				     const std::string& label,int ns)
		{
			io_.read(pSi,label,ns);
			rewind();
			pS.resize(pSi.size());
			for (size_t i=0;i<pSi.size();i++) pS[pSi[i]]=i;
			
		}

		void getTransform(MatrixType& transform,int ns)
		{
			io_.readMatrix(transform,"#TRANSFORM_sites",ns);
			rewind();
		}

		void getElectrons(std::vector<IndexType>& electrons,int ns)
		{
			io_.read(electrons,"#ELECTRONS_sites=",ns);
			rewind();
		}

		void getWaveFunction(VectorType& wavefunction,size_t ns)
		{
			io_.readSparseVector(wavefunction,"#WAVEFUNCTION_sites=",ns);
			rewind();
		}
		
		void getTimeVector(VectorType& wavefunction,size_t ns)
		{
			io_.readSparseVector(wavefunction,"targetVector0",ns);
			rewind();
		}

		std::string filename_; 
		typename IoType::In io_;
		std::vector<std::vector<IndexType> >	SpermutationInverse_,Spermutation_,
  							SEpermutationInverse_,SEpermutation_,
	 						electrons_;
		std::vector<MatrixType> transform_;
		std::vector<VectorType> wavefunction_;
		std::vector<VectorType> psiTimeVector_;
		size_t currentPos_;
		bool verbose_;
		size_t nf_;
		size_t stepTimes_;
	};  //Precomputed

	template<typename RealType1,typename FieldType1,typename IoType1,
 		typename MatrixType1,template<typename> class VectorTemplate1>
	std::ostream& operator<<(std::ostream& os,
		Precomputed<RealType1,FieldType1,IoType1,MatrixType1,VectorTemplate1>& p)
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
