
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
#include "VectorWithOffsets.h"

namespace Dmrg {
	
	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
		template<typename,typename> class InternalProductTemplate,
		template<typename,typename> class WaveFunctionTransfTemplate,
		typename ModelType_,
		typename ConcurrencyType_,
		typename IoType_,
		template<typename> class VectorWithOffsetTemplate>
	class DynamicTargetting  {
	public:
		
		typedef ModelType_ ModelType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef IoType_ IoType;
		typedef typename ModelType::RealType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef InternalProductTemplate<RealType,ModelType> InternalProductType;
		typedef typename ModelType::OperatorsType OperatorsType;
		//typedef typename OperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename ModelType::ModelHelperType ModelHelperType;
		typedef typename ModelHelperType::LeftRightSuperType
			LeftRightSuperType;
		typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
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
		typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
		typedef TimeSerializer<RealType,VectorWithOffsetType> TimeSerializerType;
		typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
		typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
		typedef typename LanczosSolverType::DenseMatrixType DenseMatrixType;
		
		enum {DISABLED,OPERATOR,CONVERGING};
		enum {	EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
				EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
				INFINITE=WaveFunctionTransfType::INFINITE};

		static const size_t parallelRank_ = 0; // DYNT needs to support concurrency FIXME
		
		
		DynamicTargetting(
				const LeftRightSuperType& lrs,
				const ModelType& model,
				const TargettingParamsType& tstStruct,
				const WaveFunctionTransfType& wft)
		:	
		 	stage_(tstStruct.sites.size(),DISABLED),
		 	lrs_(lrs),
		 	model_(model),
		 	tstStruct_(tstStruct),
		 	waveFunctionTransformation_(wft),
		 	progress_("DynamicTargetting",0),
		 	applyOpLocal_(lrs),
		 	gsWeight_(1.0)

		{
			if (!wft.isEnabled()) throw std::runtime_error(" DynamicTargetting "
					"needs an enabled wft\n");
		}

		RealType weight(size_t i) const
		{
			if (allStages(DISABLED)) throw std::runtime_error(
					"DynTarget: What are you doing here?\n");
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
				size_t loopNumber)
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
		//		// always print to keep observer driver in sync
		//		if (needsPrinting) {
		//			zeroOutVectors();
		//			printVectors(block);
		//		}
				return;
			}

			calcDynVectors(Eg,phiNew,direction);

			//cocoon(direction,block); // in-situ

			//if (needsPrinting) printVectors(block); // for post-processing
		}
		

		void initialGuess(VectorWithOffsetType& v) const
		{
			waveFunctionTransformation_.setInitialVector(v,psi_,lrs_);
			if (!allStages(CONVERGING)) return;
			std::vector<VectorWithOffsetType> vv(targetVectors_.size());
			for (size_t i=0;i<targetVectors_.size();i++) {
				waveFunctionTransformation_.setInitialVector(vv[i],
						targetVectors_[i],lrs_);
				if (std::norm(vv[i])<1e-6) continue;
				VectorWithOffsetType w= weight_[i]*vv[i];
				v += w;
			}
		}
		
		const LeftRightSuperType& leftRightSuper() const { return lrs_; }

		template<typename IoOutputType>
		void save(const std::vector<size_t>& block,IoOutputType& io) const
		{
			std::cerr<<"WARNING: DynTarget won't save\n";
			return;
			if (block.size()!=1) throw std::runtime_error(
					"DynamicTargetting only supports blocks of size 1\n");

			TimeSerializerType ts(0.0,block[0],targetVectors_);
			ts.save(io);
		}
		

		void load(const std::string& f)
		{
			std::cerr<<"WARNING: No load implemented for DynamicTargetting\n";
		}
		
		
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
				FermionSign fs(lrs_.left(),tstStruct_.electrons);
				applyOpLocal_(phiNew,phiOld,tstStruct_.aOperators[i],
						fs,systemOrEnviron);
				RealType norma = std::norm(phiNew);
				if (norma==0) throw std::runtime_error(
						"Norm of phi is zero\n");
				//std::cerr<<"Norm of phi="<<norma<<" when i="<<i<<"\n";
				
			} else if (stage_[i]== CONVERGING) {
				
				std::ostringstream msg;
				msg<<"I'm calling the WFT now";
				progress_.printline(msg,std::cout);

				if (tstStruct_.aOperators.size()==1)
					guessPhiSectors(phiNew,i,systemOrEnviron);
				else phiNew.populateSectors(lrs_.super());

				// OK, now that we got the partition number right, let's wft:
				waveFunctionTransformation_.setInitialVector(
						phiNew,targetVectors_[0],lrs_); // generalize for su(2)
				phiNew.collapseSectors();
				
			} else {
				throw std::runtime_error("It's 5 am, do you know what line "
					" your code is exec-ing?\n");
			}
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
		
		void calcDynVectors(
				RealType Eg,
				const VectorWithOffsetType&phi,
				size_t systemOrEnviron)
		{
			for (size_t i=0;i<phi.sectors();i++) {
				VectorType sv;
				size_t i0 = phi.sector(i);
				phi.extract(sv,i0);
				DenseMatrixType V;
				size_t p = lrs_.super().findPartitionNumber(phi.offset(i0));
				getLanczosVectors(V,Eg,sv,p);
				if (i==0) {
					targetVectors_.resize(V.n_col());
					for (size_t j=0;j<targetVectors_.size();j++)
						targetVectors_[j] = phi;
				}
				setLanczosVectors(V,i0);
			}
			setWeights();
		}

		void getLanczosVectors(
				DenseMatrixType& V,
				const RealType& Eg,
				const VectorType& sv,
				size_t p) const
		{
			typename ModelType::ModelHelperType modelHelper(
					p,lrs_,model_.orbitals());
			typedef typename LanczosSolverType::LanczosMatrixType
					LanczosMatrixType;
			LanczosMatrixType h(&model_,&modelHelper);

			RealType eps= ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;

			//srand48(3243447);
			LanczosSolverType lanczosSolver(h,iter,eps,parallelRank_);

			TridiagonalMatrixType ab;

			lanczosSolver.tridiagonalDecomposition(sv,ab,V);
			calcIntensity(Eg,sv,V,ab);
		}
		
		void calcIntensity(
				const RealType& Eg,
				const VectorType& sv,
				const DenseMatrixType& V,
				const TridiagonalMatrixType& ab) const
		{
			PsimagLite::Matrix<RealType> S;
			ab.buildDenseMatrix(S);
			std::vector<RealType> eigs(S.n_row());
			diag(S,eigs,'V');
			RealType delta= 0.04;
			for (RealType omega = 0;omega <3;omega+=0.02) {
				ComplexType z(omega,delta);
				ComplexType res = calcIntensity(Eg,sv,V,eigs,S,z);
				std::cout<<omega<<" "<<res<<"\n";
			}
		}

		// FIXME: Needs optimization
		ComplexType calcIntensity(
				const RealType& Eg,
				const VectorType& sv,
				const DenseMatrixType& V,
				const std::vector<RealType>& eigs,
				const PsimagLite::Matrix<RealType>& S,
				const ComplexType& z) const

		{
			RealType tmp1 = 0;
			for (size_t m=0;m<sv.size();m++)
				tmp1 += V(m,0)*sv[m];

			RealType tmp2 = 0;
			for (size_t k2=0;k2<sv.size();k2++)
				tmp2 += std::conj(sv[k2]*V(k2,0));

			ComplexType sum = 0;
			for (size_t l=0;l<S.n_row();l++)
				sum += std::conj(S(0,l))*S(0,l)/(z-eigs[l]+Eg);

			return sum*tmp1*tmp2;
		}

		void setLanczosVectors(
				const DenseMatrixType& V,
				size_t i0)
		{
			for (size_t i=0;i<targetVectors_.size();i++) {
				VectorType tmp(V.n_row());
				for (size_t j=0;j<tmp.size();j++) tmp[j] = V(j,i);
				targetVectors_[i].setDataInSector(tmp,i0);
			}
		}


		void guessPhiSectors(VectorWithOffsetType& phi,size_t i,size_t systemOrEnviron)
		{
			FermionSign fs(lrs_.left(),tstStruct_.electrons);
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
		
		void setWeights()
		{
			RealType sum  = 0;
			weight_.resize(targetVectors_.size());
			for (size_t r=0;r<weight_.size();r++) {
				weight_[r] =0;
				for (size_t i=0;i<targetVectors_[0].sectors();i++) {
					VectorType v,w;
					size_t i0 = targetVectors_[0].sector(i);
					targetVectors_[0].extract(v,i0);
					targetVectors_[r].extract(w,i0);
					weight_[r] += dynWeightOf(v,w);
				}
				sum += weight_[r];
			}
			for (size_t r=0;r<weight_.size();r++) weight_[r] = 0.5/sum;
			gsWeight_ = 0.5;

		}

		RealType dynWeightOf(VectorType& v,const VectorType& w) const
		{
			RealType sum = 0;
			for (size_t i=0;i<v.size();i++)
				sum += utils::square(std::real(v[i]*w[i]));
			return sum;
		}

		void zeroOutVectors()
		{
			for (size_t i=0;i<targetVectors_.size();i++)
				targetVectors_[i].resize(lrs_.super().size());
		}
		

		//void printHeader()
		//{
		//	io_.print(tstStruct_);
		//	std::string label = "omega";
		//	std::string s = "Omega=" + utils::ttos(currentOmega_);
		//	io_.printline(s);
		//	label = "weights";
		//	io_.printVector(weight_,label);
		//	s = "GsWeight="+utils::ttos(gsWeight_);
		//	io_.printline(s);
		//}
		
		std::vector<size_t> stage_;
		VectorWithOffsetType psi_;
		const LeftRightSuperType& lrs_;
		const ModelType& model_;
		const TargettingParamsType& tstStruct_;
		const WaveFunctionTransfType& waveFunctionTransformation_;
		PsimagLite::ProgressIndicator progress_;
		ApplyOperatorType applyOpLocal_;
		RealType gsWeight_;
		std::vector<VectorWithOffsetType> targetVectors_;
		std::vector<RealType> weight_;
		//typename IoType::Out io_;
		

	}; // class DynamicTargetting
	
	template<
	template<typename,typename,typename> class LanczosSolverTemplate,
	template<typename,typename> class InternalProductTemplate,
	template<typename,typename> class WaveFunctionTransfTemplate,
	typename ModelType_,
	typename ConcurrencyType_,
	typename IoType_,
	template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const DynamicTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,ConcurrencyType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		os<<"DT=NothingToSeeHereYet\n";
		return os;
	}
	
} // namespace
#endif // DYNAMICTARGETTING_H

