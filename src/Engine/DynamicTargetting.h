
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

#ifndef DYNAMICTARGETTING_H
#define DYNAMICTARGETTING_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "TimeSerializer.h"
#include "DynamicDmrgParams.h"
#include "DynamicFunctional.h"
#include "Minimizer.h"
#include "VectorWithOffsets.h"

namespace Dmrg {
	
	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
		template<typename,typename> class InternalProductTemplate,
		typename WaveFunctionTransformationType_,
		typename ModelType_,
		typename ConcurrencyType_,
		typename IoType_,
		template<typename> class VectorWithOffsetTemplate>
	class DynamicTargetting  {
	public:
		
		typedef WaveFunctionTransformationType_ WaveFunctionTransformationType;
		typedef ModelType_ ModelType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef IoType_ IoType;
		typedef typename ModelType::RealType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef InternalProductTemplate<RealType,ModelType> InternalProductType;
		typedef typename ModelType::OperatorsType OperatorsType;
		//typedef typename OperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename ModelType::MyBasisWithOperators BasisWithOperatorsType;
		//typedef std::vector<RealType> VectorType;
		//typedef PsimagLite::Matrix<ComplexType> ComplexMatrixType;
		//typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
		typedef typename BasisWithOperatorsType::OperatorType OperatorType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef DynamicDmrgParams<ModelType> TargettingParamsType;
		typedef typename BasisType::BlockType BlockType;
		typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
		typedef typename VectorWithOffsetType::VectorType VectorType;
		typedef LanczosSolverTemplate<RealType,InternalProductType,VectorType> LanczosSolverType;
		typedef VectorType TargetVectorType;
		typedef ApplyOperatorLocal<BasisWithOperatorsType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
		typedef TimeSerializer<RealType,VectorWithOffsetType> TimeSerializerType;
		
		
		enum {DISABLED,OPERATOR,CONVERGING};
		enum {	EXPAND_ENVIRON=WaveFunctionTransformationType::EXPAND_ENVIRON,
				EXPAND_SYSTEM=WaveFunctionTransformationType::EXPAND_SYSTEM,
				INFINITE=WaveFunctionTransformationType::INFINITE};

		static const size_t parallelRank_ = 0; // DYNT needs to support concurrency FIXME
		
		
		DynamicTargetting(
				const BasisWithOperatorsType& basisS,
				const BasisWithOperatorsType& basisE,
				const BasisType& basisSE,
				const ModelType& model,
				const TargettingParamsType& tstStruct,
				const WaveFunctionTransformationType& wft)
		:	
		 	stage_(tstStruct.sites.size(),DISABLED),
		 	basisS_(basisS),
		 	basisE_(basisE),
		 	basisSE_(basisSE),
		 	model_(model),
		 	tstStruct_(tstStruct),
		 	waveFunctionTransformation_(wft),
		 	progress_("DynamicTargetting",0),
		 	currentOmega_(tstStruct_.omega),
		 	targetVectors_(3),
		 	weight_(targetVectors_.size()),
		 	io_(tstStruct_.filename,parallelRank_),
		 	applyOpLocal_(basisS,basisE,basisSE)
		 	

		{
			if (!wft.isEnabled()) throw std::runtime_error(" DynamicTargetting "
					"needs an enabled wft\n");
			RealType sum = 0;
			size_t n = weight_.size();
			for (size_t i=0;i<n;i++) {
				weight_[i] = 1.0/(n+1);
				sum += weight_[i];
			}

			gsWeight_=1.0-sum;
			sum += gsWeight_;
			if (fabs(sum-1.0)>1e-5) throw std::runtime_error("Weights don't amount to one\n");
			printHeader();
		}
		
		
		

		RealType weight(size_t i) const
		{
			if (allStages(DISABLED)) throw std::runtime_error("TST: What are you doing here?\n");
			return weight_[i];
			//return 1.0;
		}
		

		RealType gsWeight() const
		{
			if (allStages(DISABLED)) return 1.0;
			return gsWeight_;
		}
		

		RealType normSquared(size_t i) const
		{
			// call to mult will conjugate one of the vector
			return std::real(multiply(targetVectors_[i],targetVectors_[i]));
		}
		

		template<typename SomeBasisType>
		void setGs(const std::vector<TargetVectorType>& v,
				const SomeBasisType& someBasis)
		{
			psi_.set(v,someBasis);
		}
		

		const RealType& operator[](size_t i) const { return psi_[i]; }
					
		RealType& operator[](size_t i) { return psi_[i]; }
		

		const VectorWithOffsetType& gs() const { return psi_; }
		

		bool includeGroundStage() const {return true; }
		

		size_t size() const
		{
			if (allStages(DISABLED)) return 0;
			return targetVectors_.size();
		}
		

		const VectorWithOffsetType& operator()(size_t i) const
		{
			return targetVectors_[i];
		}
		

		void evolve(RealType Eg,size_t direction,const BlockType& block,
				size_t loopNumber, bool needsPrinting)
		{
			size_t count =0;
			VectorWithOffsetType phiOld = psi_;
			VectorWithOffsetType phiNew;
			size_t max = tstStruct_.sites.size();

			if (noStageIs(DISABLED)) max = 1;

			// Loop over each operator that needs to be applied
			// in turn to the g.s.
			for (size_t i=0;i<max;i++) {
				count += evolve(i,phiNew,phiOld,Eg,direction,block,loopNumber,max-1);
				phiOld = phiNew;
			}

			if (count==0) {
				// always print to keep observer driver in sync
				if (needsPrinting) {
					zeroOutVectors();
					printVectors(block);
				}
				return;
			}

			ComplexType val = calcDynVectors(Eg,phiNew,direction);

			cocoon(val,direction,block); // in-situ

			if (needsPrinting) printVectors(block); // for post-processing
		}
		

		void initialGuess(VectorWithOffsetType& v) const
		{
			waveFunctionTransformation_.setInitialVector(v,psi_,basisS_,basisE_,basisSE_);
			if (!allStages(CONVERGING)) return;
			std::vector<VectorWithOffsetType> vv(targetVectors_.size());
			for (size_t i=0;i<targetVectors_.size();i++) {
				waveFunctionTransformation_.setInitialVector(vv[i],
						targetVectors_[i],basisS_,basisE_,basisSE_);
				if (std::norm(vv[i])<1e-6) continue;
				VectorWithOffsetType w= weight_[i]*vv[i];
				v += w;
			}
		}
		

		const BasisType& basisSE() const { return basisSE_; }

		const BasisWithOperatorsType& basisS() const { return basisS_; }

		const BasisWithOperatorsType& basisE() const { return basisE_; }
		
		
	private:
		

		size_t evolve(
				size_t i,
				VectorWithOffsetType& phiNew,
				VectorWithOffsetType& phiOld,
				RealType Eg,
				size_t direction,
				const BlockType& block,
				size_t loopNumber,
				size_t lastI)
		{
			
			if (tstStruct_.startingLoops[i]>loopNumber || direction==INFINITE) return 0;
			
			
			if (block.size()!=1) throw
			std::runtime_error("DynamicTargetting::evolve(...):"
				" blocks of size != 1 are unsupported (sorry)\n");
			size_t site = block[0];
			
			
			if (site != tstStruct_.sites[i] && stage_[i]==DISABLED) return 0;
			
			
			if (site == tstStruct_.sites[i] && stage_[i]==DISABLED) stage_[i]=OPERATOR;
			else stage_[i]=CONVERGING;
			if (stage_[i] == OPERATOR) checkOrder(i);
			
			
			std::ostringstream msg;
			msg<<"Evolving, stage="<<getStage(i)<<" site="<<site<<" loopNumber="<<loopNumber;
			msg<<" Eg="<<Eg;
			progress_.printline(msg,std::cout);
							
			// phi = A|psi>
			computePhi(i,phiNew,phiOld,direction);
			
			return 1;
		}
		

		void computePhi(size_t i,VectorWithOffsetType& phiNew,
			VectorWithOffsetType& phiOld,size_t systemOrEnviron)
		{
			if (stage_[i]==OPERATOR) {
				
				std::ostringstream msg;
				msg<<"I'm applying a local operator now";
				progress_.printline(msg,std::cout);
				FermionSign fs(basisS_,tstStruct_.electrons);
				applyOpLocal_(phiNew,phiOld,tstStruct_.aOperators[i],fs,systemOrEnviron);
				RealType norma = std::norm(phiNew);
				if (norma==0) throw std::runtime_error("Norm of phi is zero\n");
				//std::cerr<<"Norm of phi="<<norma<<" when i="<<i<<"\n";
				
			} else if (stage_[i]== CONVERGING) {
				
				std::ostringstream msg;
				msg<<"I'm calling the WFT now";
				progress_.printline(msg,std::cout);

				if (tstStruct_.aOperators.size()==1) guessPhiSectors(phiNew,i,systemOrEnviron);
				else phiNew.populateSectors(basisSE_);

				// OK, now that we got the partition number right, let's wft:
				waveFunctionTransformation_.setInitialVector(phiNew,targetVectors_[0],
					basisS_,basisE_,basisSE_); // generalize for su(2)
				phiNew.collapseSectors();
				
			} else {
				throw std::runtime_error("It's 5 am, do you know what line "
					" your code is exec-ing?\n");
			}
		}
		

		void cocoon(ComplexType& val,size_t direction,const BlockType& block) const
		{
			size_t site = block[0];
			std::cerr<<"-------------&*&*&* Cocoon output starts\n";
			test(psi_,psi_,direction,"<PSI|A|PSI>",site);
			std::cerr<<"OMEGA "<<currentOmega_<<" "<<std::imag(val)<<" "<<std::real(val)<<" "<<site<<"\n";
			for (size_t j=0;j<targetVectors_.size();j++) {
				std::string s = "<P"+utils::ttos(j)+"|A|P"+utils::ttos(j)+">";
				test(targetVectors_[j],targetVectors_[0],direction,s,site);
			}
			std::cerr<<"-------------&*&*&* Cocoon output ends\n";
		}
		

		void checkOrder(size_t i) const
		{
			if (i==0) return;
			for (size_t j=0;j<i;j++) {
				if (stage_[j] == DISABLED) {
					std::string s ="TST:: Seeing dynamic site "+utils::ttos(tstStruct_.sites[i]);
					s =s + " before having seen";
					s = s + " site "+utils::ttos(j);
					s = s +". Please order your dynamic sites in order of appearance.\n";
					throw std::runtime_error(s);
				}
			}
		}
		

		bool allStages(size_t x) const
		{
			for (size_t i=0;i<stage_.size();i++)
				if (stage_[i]!=x) return false;
			return true;
		}
		

		bool noStageIs(size_t x) const
		{
			for (size_t i=0;i<stage_.size();i++)
				if (stage_[i]==x) return false;
			return true;
		}
		

		std::string getStage(size_t i) const
		{
			switch (stage_[i]) {
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
		

		ComplexType calcDynVectors(
				RealType Eg,
				const VectorWithOffsetType& phi,
				size_t systemOrEnviron)
		{
			RealType retIm = minimizeFunctional(targetVectors_[1],Eg,phi,systemOrEnviron);
			obtainXA(targetVectors_[2],targetVectors_[1],Eg);
			RealType retRe = 0; //-std::real(targetVectors_[2]*phi)/M_PI; // Eq.~(12a)
			targetVectors_[0] = phi;
			areAllTargetsSensible();
			return ComplexType(retRe,retIm);
		}
		

		RealType minimizeFunctional(
				VectorWithOffsetType& psiMin,
				RealType Eg,
				const VectorWithOffsetType&phi,
				size_t systemOrEnviron)
		{
			VectorWithOffsetType phiCopy = phi;
			psiMin = phi;
			RealType ret = 0;
			for (size_t i=0;i<phiCopy.sectors();i++) {
				VectorType sv;
				size_t ii = phiCopy.sector(i);
				psiMin.extract(sv,ii);
				if (sv.size()==0) throw std::runtime_error("Non-zero sector is zero!\n");
				ret += minimizeFunctional(sv,Eg,phi,ii);
				psiMin.setDataInSector(sv,ii);
			}
			return ret;
		}
		

		RealType minimizeFunctional(VectorType& sv,RealType Eg,const VectorWithOffsetType& phi,size_t ind)
		{
			size_t p = basisSE_.findPartitionNumber(phi.offset(ind));
			typename ModelType::ModelHelperType modelHelper(p,basisSE_,basisS_,basisE_,model_.orbitals());
			typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
			LanczosMatrixType h(&model_,&modelHelper);
			typedef DynamicFunctional<RealType,LanczosMatrixType,VectorType> DynamicFunctionalType;
			VectorType aVector;
			phi.extract(aVector,ind);
			DynamicFunctionalType wFunctional(h,aVector,currentOmega_,Eg,tstStruct_.eta);
			size_t maxIter = 1000;

			PsimagLite::Minimizer<RealType,DynamicFunctionalType> min(wFunctional,maxIter);
			std::vector<RealType> svReal(sv.size());
			//wFunctional.packComplexToReal(svReal,sv);
			for (size_t i=0;i<svReal.size();i++) svReal[i]=drand48();
			RealType norma = std::norm(svReal);
			for (size_t i=0;i<svReal.size();i++) svReal[i]/=norma;

			int iter = -1;
			RealType delta = 0.01;
			RealType tolerance = 1e-5;
			size_t counter = 0;
			while (iter<0 && counter<100) {
				iter = min.simplex(svReal,delta,tolerance);
				if (iter>=0) {
					std::cerr<<"delta="<<delta<<" tolerance="<<tolerance<<" counter="<<counter<<"\n";
				}
				//delta /= 2;
				tolerance *= 1.5;
				counter++;
			}
			if (iter<0) {
				std::cerr<<"delta="<<delta<<" tol="<<tolerance<<"\n";
				throw std::runtime_error
					("DynTargetting::minimizeFunctional(...):No minimum found\n");
			}
			//wFunctional.packRealToComplex(sv,svReal);
			sv = svReal;
			return  -wFunctional(svReal)/(M_PI*tstStruct_.eta); // Eq.~(16)
		}
		

		void obtainXA(
				VectorWithOffsetType& xa,
				const VectorWithOffsetType& ya,
				RealType Eg)
		{
			xa = ya;
			for (size_t i=0;i<ya.sectors();i++) {
				size_t ii = ya.sector(i);
				obtainXA(xa,Eg,ya,ii);
			}
		}
		

		void obtainXA(VectorWithOffsetType& xa,RealType Eg,const VectorWithOffsetType& ya,size_t i)
		{
			size_t p = basisSE_.findPartitionNumber(ya.offset(i));
			typename ModelType::ModelHelperType modelHelper(p,basisSE_,basisS_,basisE_,model_.orbitals());
			typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
			LanczosMatrixType h(&model_,&modelHelper);
			VectorType yaThisSector;
			ya.extract(yaThisSector,i);
			VectorType sv(yaThisSector.size(),0.0);
			h.matrixVectorProduct(sv,yaThisSector); // sv = H * yaThisSector
			RealType factor =  (Eg+currentOmega_);
			sv -= (yaThisSector * factor);
			sv *= (1/tstStruct_.eta);
			xa.setDataInSector(sv,i);
		}
		

		void guessPhiSectors(VectorWithOffsetType& phi,size_t i,size_t systemOrEnviron)
		{
			FermionSign fs(basisS_,tstStruct_.electrons);
			if (allStages(CONVERGING)) {
				VectorWithOffsetType tmpVector = psi_;
				for (size_t j=0;j<tstStruct_.aOperators.size();j++) {
					applyOpLocal_(phi,tmpVector,tstStruct_.aOperators[j],fs,
							systemOrEnviron);
					tmpVector = phi;
				}
				return;
			}
			applyOpLocal_(phi,psi_,tstStruct_.aOperators[i],fs,
					systemOrEnviron);
		}
		

		void zeroOutVectors()
		{
			for (size_t i=0;i<targetVectors_.size();i++)
				targetVectors_[i].resize(basisSE_.size());
		}
		

		void printVectors(const std::vector<size_t>& block)
		{
			if (block.size()!=1) throw std::runtime_error(
					"DynamicTargetting only supports blocks of size 1\n");

			TimeSerializerType ts(currentOmega_,block[0],targetVectors_);
			ts.save(io_);
		}
		

		void printHeader()
		{
			io_.print(tstStruct_);
			std::string label = "omega";
			std::string s = "Omega=" + utils::ttos(currentOmega_);
			io_.printline(s);
			label = "weights";
			io_.printVector(weight_,label);
			s = "GsWeight="+utils::ttos(gsWeight_);
			io_.printline(s);
		}
		

		void test(
				const VectorWithOffsetType& src1,
				const VectorWithOffsetType& src2,
				size_t systemOrEnviron,
				const std::string& label,
				size_t site) const
		{
			VectorWithOffsetType dest;
			OperatorType A = tstStruct_.aOperators[0];
			CrsMatrix<RealType> tmpC(model_.getOperator("c",0,0));
			/*CrsMatrix<ComplexType> tmpCt;
						transposeConjugate(tmpCt,tmpC);
						multiply(A.data,tmpCt,tmpC);*/
			A.fermionSign = 1;
			A.data.makeDiagonal(tmpC.rank(),1.0);
			FermionSign fs(basisS_,tstStruct_.electrons);
			applyOpLocal_(dest,src1,A,fs,systemOrEnviron);

			RealType sum = 0;
			for (size_t ii=0;ii<dest.sectors();ii++) {
				size_t i = dest.sector(ii);
				size_t offset1 = dest.offset(i);
				for (size_t jj=0;jj<src2.sectors();jj++) {
					size_t j = src2.sector(jj);
					size_t offset2 = src2.offset(j);
					if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
					for (size_t k=0;k<dest.effectiveSize(i);k++)
						sum+= dest[k+offset1] * std::conj(src2[k+offset2]);
				}
			}
			std::cerr<<site<<" "<<sum<<" "<<" "<<currentOmega_;
			std::cerr<<" "<<label<<std::norm(src1)<<" "<<std::norm(src2)<<" "<<std::norm(dest)<<"\n";
		}
		

		void areAllTargetsSensible() const
		{
			for (size_t i=0;i<targetVectors_.size();i++)
				isThisTargetSensible(i);
		}
		

		void isThisTargetSensible(size_t i) const
		{
			RealType norma = std::norm(targetVectors_[i]);
			if (norma<1e-6) throw std::runtime_error("Norma is zero\n");
		}
		
		
		
		std::vector<size_t> stage_;
		VectorWithOffsetType psi_;
		const BasisWithOperatorsType& basisS_;
		const BasisWithOperatorsType& basisE_;
		const BasisType& basisSE_;
		const ModelType& model_;
		const TargettingParamsType& tstStruct_;
		const WaveFunctionTransformationType& waveFunctionTransformation_;
		PsimagLite::ProgressIndicator progress_;
		RealType currentOmega_;
		std::vector<VectorWithOffsetType> targetVectors_;
		std::vector<RealType> weight_;
		RealType gsWeight_;
		typename IoType::Out io_;
		ApplyOperatorType applyOpLocal_;
		
	}; // class DynamicTargetting
	
} // namespace
#endif // DYNAMICTARGETTING_H

