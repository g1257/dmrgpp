
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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

/*! \file CommonTargetting.h
 *
 * Functionality used by many targetting classes
 *
 */

#ifndef COMMON_TARGETTING_H
#define COMMON_TARGETTING_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "DynamicSerializer.h"
#include "DynamicDmrgParams.h"
#include "VectorWithOffsets.h"
#include "ContinuedFraction.h"
#include <cassert>

namespace Dmrg {

	template<typename ModelType,
	         typename TargettingParamsType,
	         typename WaveFunctionTransfType,
	         typename VectorWithOffsetType,
			 typename LanczosSolverType>
	class CommonTargetting  {

	public:

		typedef typename ModelType::RealType RealType;
		typedef typename ModelType::ModelHelperType ModelHelperType;
		typedef typename ModelHelperType::LeftRightSuperType
		                                  LeftRightSuperType;
		typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
		typedef typename LanczosSolverType::PostProcType PostProcType;
		typedef typename VectorWithOffsetType::VectorType VectorType;
		typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
		typedef DynamicSerializer<VectorWithOffsetType,PostProcType> DynamicSerializerType;
		typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
		typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename BasisWithOperatorsType::OperatorType OperatorType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename BasisWithOperatorsType::BasisDataType BasisDataType;
		typedef typename BasisType::BlockType BlockType;
		typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,VectorType> ApplyOperatorType;

		enum {DISABLED,OPERATOR,CONVERGING};
		enum {
			EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE
		};

		CommonTargetting(const LeftRightSuperType& lrs,
		                 const ModelType& model,
						 const TargettingParamsType& tstStruct)
		: lrs_(lrs),
		  model_(model),
		  tstStruct_(tstStruct),
		  applyOpLocal_(lrs),
		  progress_("CommonTargetting")
		{
		}

		RealType normSquared(const VectorWithOffsetType& v) const
		{
			// call to mult will conjugate one of the vector
			return std::real(multiply(v,v));
		}
		
		template<typename IoOutputType>
		void save(const typename PsimagLite::Vector<SizeType>::Type& block,
		          IoOutputType& io,
				  const PostProcType& cf,
				  const typename PsimagLite::Vector<VectorWithOffsetType>::Type& targetVectors) const
		{
			DynamicSerializerType dynS(cf,block[0],targetVectors);
			dynS.save(io);
		}

		template<typename IoInputType>
		void load(IoInputType& io,typename PsimagLite::Vector<VectorWithOffsetType>::Type& targetVectors)
		{
			DynamicSerializerType dynS(io,IoInputType::In::LAST_INSTANCE);
			for (SizeType i=0;i<targetVectors.size();i++)
				targetVectors[i] = dynS.vector(i);
			
		}

		void checkOrder(SizeType i,const typename PsimagLite::Vector<SizeType>::Type& stage) const
		{
			if (i==0) return;
			for (SizeType j=0;j<i;j++) {
				if (stage[j] == DISABLED) {
					PsimagLite::String s ="TST:: Seeing dynamic site "+ttos(tstStruct_.sites[i]);
					s =s + " before having seen";
					s = s + " site "+ttos(j);
					s = s +". Please order your dynamic sites in order of appearance.\n";
					throw PsimagLite::RuntimeError(s);
				}
			}
		}

		bool allStages(SizeType x,const typename PsimagLite::Vector<SizeType>::Type& stage) const
		{
			for (SizeType i=0;i<stage.size();i++)
				if (stage[i]!=x) return false;
			return true;
		}

		bool noStageIs(SizeType x,const typename PsimagLite::Vector<SizeType>::Type& stage) const
		{
			for (SizeType i=0;i<stage.size();i++)
				if (stage[i]==x) return false;
			return true;
		}

		PsimagLite::String getStage(SizeType i,const typename PsimagLite::Vector<SizeType>::Type& stage) const
		{
			switch (stage[i]) {
			case DISABLED:
				return "Disabled";
				break;
			case OPERATOR:
				return "Applying operator for the first time";
				break;
			case CONVERGING:
				return "Converging DDMRG";
				break;
			}
			return "undefined";
		}

		void initialGuess(VectorWithOffsetType& v,
		                  const WaveFunctionTransfType& wft,
		                  const VectorWithOffsetType& psi,
		                  const typename PsimagLite::Vector<SizeType>::Type& stage,
		                  const typename PsimagLite::Vector<RealType>::Type& weights,
						  const typename PsimagLite::Vector<SizeType>::Type& block,
						  const typename PsimagLite::Vector<VectorWithOffsetType>::Type& targetVectors) const
		{
			typename PsimagLite::Vector<SizeType>::Type nk;
			setNk(nk,block);
			wft.setInitialVector(v,psi,lrs_,nk);
			if (!allStages(CONVERGING,stage)) return;
			typename PsimagLite::Vector<VectorWithOffsetType>::Type vv(targetVectors.size());
			for (SizeType i=0;i<targetVectors.size();i++) {
				wft.setInitialVector(vv[i],targetVectors[i],lrs_,nk);
				if (std::norm(vv[i])<1e-6) continue;
				VectorWithOffsetType w= weights[i]*vv[i];
				v += w;
			}
		}

		void findElectronsOfOneSite(typename PsimagLite::Vector<SizeType>::Type& electrons,SizeType site) const
		{
			typename PsimagLite::Vector<SizeType>::Type block(1,site);
			typename ModelType::HilbertBasisType basis;
			typename PsimagLite::Vector<SizeType>::Type quantumNumbs;
			model_.setNaturalBasis(basis,quantumNumbs,block);
			model_.findElectrons(electrons,basis,site);
		}

		void noCocoon(const PsimagLite::String& msg) const
		{
			std::cout<<"-------------&*&*&* In-situ measurements start\n";
			std::cout<<"----- NO IN-SITU MEAS. POSSIBLE, reason="<<msg<<"\n";
			std::cout<<"-------------&*&*&* In-situ measurements end\n";
		}

		// in situ computation:
		void cocoon(SizeType direction,SizeType site,const VectorWithOffsetType& psi) const
		{
			int fermionSign1 = 1;
			const std::pair<SizeType,SizeType> jm1(0,0);
			RealType angularFactor1 = 1.0;
			typename OperatorType::Su2RelatedType su2Related1;

			std::cout<<"-------------&*&*&* In-situ measurements start\n";

			typename PsimagLite::Vector<PsimagLite::String>::Type vecStr;
			PsimagLite::tokenizer(model_.params().insitu,vecStr,",");
			for (SizeType i=0;i<vecStr.size();i++) {
				const PsimagLite::String& opLabel = vecStr[i];
				OperatorType nup;
				if (!fillOperatorFromFile(nup,opLabel)) {
					PsimagLite::CrsMatrix<RealType> tmpC(model_.naturalOperator(opLabel,site,0));
					nup = OperatorType(tmpC,fermionSign1,jm1,angularFactor1,su2Related1);
				}
				PsimagLite::String tmpStr = "<PSI|" + opLabel + "|PSI>";
				test(psi,psi,direction,tmpStr,site,nup);
			}

			std::cout<<"-------------&*&*&* In-situ measurements end\n";
		}

		void computeCorrection(VectorWithOffsetType& v,
							   SizeType direction,
							   const BlockType& block1,
							   const VectorWithOffsetType& psi) const
		{
			// operators in the one-site basis:
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			SparseMatrixType hmatrix;
			BasisDataType q;

			RealType time = 0;
			model_.setNaturalBasis(creationMatrix,hmatrix,q,block1,time);

			FermionSign fs(lrs_.left(),q.electrons);
			for (SizeType j=0;j<creationMatrix.size();j++) {
				VectorWithOffsetType phiTemp;
				applyOpLocal_(phiTemp,psi,creationMatrix[j],
							  fs,direction);
				if (j==0) v = phiTemp;
				else v += phiTemp;
			}
		}

		void setNk(typename PsimagLite::Vector<SizeType>::Type& nk,const typename PsimagLite::Vector<SizeType>::Type& block) const
		{
			for (SizeType i=0;i<block.size();i++)
				nk.push_back(model_.hilbertSize(block[i]));
		}

	private:

		void test(const VectorWithOffsetType& src1,
				  const VectorWithOffsetType& src2,
				  SizeType systemOrEnviron,
				  const PsimagLite::String& label,
				  SizeType site,
				  const OperatorType& A) const
		{
			typename PsimagLite::Vector<SizeType>::Type electrons;
			model_.findElectronsOfOneSite(electrons,site);
			FermionSign fs(lrs_.left(),electrons);
			VectorWithOffsetType dest;
			applyOpLocal_(dest,src1,A,fs,systemOrEnviron);

			RealType sum = 0;
			for (SizeType ii=0;ii<dest.sectors();ii++) {
				SizeType i = dest.sector(ii);
				SizeType offset1 = dest.offset(i);
				for (SizeType jj=0;jj<src2.sectors();jj++) {
					SizeType j = src2.sector(jj);
					SizeType offset2 = src2.offset(j);
					if (i!=j) continue; //throw PsimagLite::RuntimeError("Not same sector\n");
					for (SizeType k=0;k<dest.effectiveSize(i);k++)
						sum+= dest[k+offset1] * std::conj(src2[k+offset2]);
				}
			}
			std::cout<<site<<" "<<sum<<" "<<" 0";
			std::cout<<" "<<label<<" "<<(src1*src2)<<"\n";
		}

		bool fillOperatorFromFile(OperatorType& nup,const PsimagLite::String& label2) const
		{
			if (label2.length()<2 || label2[0]!=':') return false;
			PsimagLite::String label = label2.substr(1,label2.length()-1);

			std::ifstream fin(label.c_str());
			if (!fin || fin.bad() || !fin.good()) return false;

			PsimagLite::String line1("");
			fin>>line1;
			if (fin.eof() || line1!="TSPOperator=raw") return false;

			line1="";
			fin>>line1;
			if (fin.eof() || line1!="RAW_MATRIX") return false;

			line1="";
			fin>>line1;
			if (fin.eof()) return false;

			PsimagLite::String line2("");
			fin>>line2;
			if (fin.eof()) return false;

			int n = atoi(line1.c_str());
			if (n!=atoi(line2.c_str()) || n<=0) return false;

			PsimagLite::Matrix<RealType> m(n,n);
			for (int i=0;i<n;i++) {
				for (int j=0;j<n;j++) {
					line1="";
					fin>>line1;
					if (fin.eof()) return false;
					m(i,j) = atof(line1.c_str());
				}
			}

			line1="";
			fin>>line1;
			if (fin.eof()) return false;
			int fermionicSign = 0;
			if (line1=="FERMIONSIGN=-1") fermionicSign = -1;
			if (line1!="FERMIONSIGN=1") fermionicSign = 1;

			if (fermionicSign==0) return false;

			line1="";
			fin>>line1;
			if (fin.eof() || line1!="JMVALUES") return false;

			const std::pair<SizeType,SizeType> jm1(0,0);
			line1="";
			fin>>line1;
			if (fin.eof()) return false;


			line1="";
			fin>>line1;
			if (fin.eof()) return false;

			typename OperatorType::Su2RelatedType su2Related1;
			RealType angularFactor1 = 1.0;
			line1="";
			fin>>line1;
			if (fin.eof() || line1.substr(0,14)!="AngularFactor=") return false;

			PsimagLite::CrsMatrix<RealType> msparse(m);
			nup = OperatorType(msparse,fermionicSign,jm1,angularFactor1,su2Related1);
			return true;
		}


		const LeftRightSuperType& lrs_;
		const ModelType& model_;
		const TargettingParamsType& tstStruct_;
		ApplyOperatorType applyOpLocal_;
		PsimagLite::ProgressIndicator progress_;
	}; // class CommonTargetting

	template<typename ModelType,
	         typename TargettingParamsType,
	         typename WaveFunctionTransfType,
	         typename VectorWithOffsetType,
	         typename LanczosSolverType
	>
	std::ostream& operator<<(std::ostream& os,
	                         const CommonTargetting<ModelType,TargettingParamsType,WaveFunctionTransfType,VectorWithOffsetType,LanczosSolverType>& tst)
	{
		os<<"DT=NothingToSeeHereYet\n";
		return os;
	}

} // namespace
/*@}*/
#endif // COMMON_TARGETTING_H

