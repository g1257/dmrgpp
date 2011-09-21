
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

#ifndef METTS_TARGETTING_H
#define METTS_TARGETTING_H
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "MettsSerializer.h"
#include "MettsParams.h"

namespace Dmrg {
	template<
			template<typename,typename,typename> class LanczosSolverTemplate,
   			template<typename,typename> class InternalProductTemplate,
	 		template<typename,typename> class WaveFunctionTransfTemplate,
    			typename ModelType_,
    			typename ConcurrencyType_,
    			typename IoType_,
       			template<typename> class VectorWithOffsetTemplate>
	class MettsTargetting  {

		public:
			typedef ModelType_ ModelType;
			typedef ConcurrencyType_ ConcurrencyType;
			typedef IoType_ IoType;
			typedef typename ModelType::RealType RealType;
			typedef InternalProductTemplate<ComplexType,ModelType>
			   InternalProductType;
			typedef typename ModelType::OperatorsType OperatorsType;
			typedef typename ModelType::ModelHelperType ModelHelperType;
			typedef typename ModelHelperType::LeftRightSuperType
			    LeftRightSuperType;
			typedef typename LeftRightSuperType::BasisWithOperatorsType
			    BasisWithOperatorsType;
			typedef std::vector<RealType> VectorType;
			typedef LanczosSolverTemplate<RealType,InternalProductType,VectorType>
			    LanczosSolverType;
			typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef MettsParams<ModelType> TargettingParamsType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
			typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
			typedef VectorType TargetVectorType;
			typedef PsimagLite::Matrix<RealType> MatrixType;
			typedef BlockMatrix<RealType,MatrixType> BlockMatrixType;
			typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
			typedef MettsSerializer<RealType,VectorWithOffsetType> MettsSerializerType;

			enum {DISABLED,OPERATOR,WFT_NOADVANCE,WFT_ADVANCE};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};

			static size_t const PRODUCT = TargettingParamsType::PRODUCT;
			static size_t const SUM = TargettingParamsType::SUM;
			static const size_t parallelRank_ = 0; // Metts needs to support concurrency FIXME

			MettsTargetting(const LeftRightSuperType& lrs,
	 		                const ModelType& model,
			                const TargettingParamsType& mettsStruct,
			                const WaveFunctionTransfType& wft)
			: stage_(mettsStruct.sites.size(),DISABLED),
			  lrs_(lrs),
			  model_(model),
			  mettsStruct_(mettsStruct),
			  waveFunctionTransformation_(wft),
			  progress_("MettsTargetting",0),
			  currentBeta_(0),
			  betas_(mettsStruct_.betaSteps),
			  weight_(mettsStruct_.betaSteps),
			  targetVectors_(mettsStruct_.betaSteps),
			  applyOpLocal_(lrs)
			{
				if (!wft.isEnabled()) throw std::runtime_error(" MettsTargetting "
							"needs an enabled wft\n");
				
				RealType tau =mettsStruct_.tau;
				RealType sum = 0;
				RealType factor = 1.0;
				size_t n = betas_.size();
				for (size_t i=0;i<n;i++) {
					betas_[i] = i*tau/(n-1);
					weight_[i] = factor/(n+4);
					sum += weight_[i];
				}
				sum -= weight_[0];
				sum -= weight_[n-1];
				weight_[0] = weight_[n-1] = 2*factor/(n+4);
				sum += weight_[n-1];
				sum += weight_[0];
				
				gsWeight_=1.0-sum;
				sum += gsWeight_;
				//for (size_t i=0;i<weight_.size();i++) sum += weight_[i];
				if (fabs(sum-1.0)>1e-5)
					throw std::runtime_error("Weights don't amount to one\n");
				//printHeader();
			}

			RealType weight(size_t i) const
			{
				if (allStages(DISABLED)) 
					throw std::runtime_error("Metts: What are you doing here?\n");
				return weight_[i];
			}

			RealType gsWeight() const
			{
				return (allStages(DISABLED)) ?  1.0 : gsWeight_;
			}

			RealType normSquared(size_t i) const
			{
				// call to mult will conjugate one of the vector
				return real(multiply(targetVectors_[i],targetVectors_[i])); 
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
				return (allStages(DISABLED)) ? 0 : targetVectors_.size();
			}

			const VectorWithOffsetType& operator()(size_t i) const
			{
				return targetVectors_[i];
			}

			const LeftRightSuperType& leftRightSuper() const
			{
				return lrs_;
			}

			void evolve(RealType Eg,
			            size_t direction,
			            const BlockType& block,
			            size_t loopNumber)
			{
				size_t count =0;
				VectorWithOffsetType phiOld = psi_;
				VectorWithOffsetType phiNew;
				VectorWithOffsetType vectorSum;
				size_t max = mettsStruct_.sites.size();
				
				if (noStageIs(DISABLED)) max = 1;
				
				// Loop over each operator that needs to be applied 
				// in turn to the g.s.
				for (size_t i=0;i<max;i++) {
					count += evolve(i,phiNew,phiOld,Eg,direction,block,loopNumber,max-1);
					if (mettsStruct_.concatenation==PRODUCT) {
						phiOld = phiNew;
					} else {
						vectorSum += phiNew;
					}
				}
				if (mettsStruct_.concatenation==SUM) phiNew = vectorSum;
				
				if (count==0) {
					// always print to keep observer driver in sync
//					if (needsPrinting) {
//						zeroOutVectors();
//						printVectors(block);
//					}
//					return;
				}
				
				calcTimeVectors(Eg,phiNew,direction);
				
				cocoon(direction,block); // in-situ
				
				//if (needsPrinting) printVectors(block); // for post-processing
			}

			
			void load(const std::string& f)
			{
				throw std::runtime_error("Metts: load() unimplemented\n");
// 				for (size_t i=0;i<stage_.size();i++) stage_[i] = WFT_NOADVANCE;
// 
// 				typename IoType::In io(f);
// 
// 				TimeSerializerType ts(io,IoType::In::LAST_INSTANCE);
// 				for (size_t i=0;i<targetVectors_.size();i++) targetVectors_[i] = ts.vector(i);
// 				currentBeta_ = ts.time();
// 
// 				psi_.load(io,"PSI");
			}
			
			
			void print(std::ostream& os) const
			{
				os<<"MettsWeightsTimeVectors=";
				for (size_t i=0;i<weight_.size();i++)
					os<<weight_[i]<<" ";
				os<<"\n";
				os<<"MettsWeightGroundState="<<gsWeight_<<"\n";
			}
			

			size_t evolve(size_t i,
			              VectorWithOffsetType& phiNew,
			              const VectorWithOffsetType& phiOld,
			              RealType Eg,
			              size_t direction,
			              const BlockType& block,
			              size_t loopNumber,
			              size_t lastI)
			{
				static size_t  timesWithoutAdvancement=0;

				if (mettsStruct_.startingLoops[i]>loopNumber || direction==INFINITE) return 0;

				if (block.size()!=1) throw 
					std::runtime_error("MettsTargetting::evolve(...):"
							" blocks of size != 1 are unsupported (sorry)\n");
				size_t site = block[0];

				if (site != mettsStruct_.sites[i] && stage_[i]==DISABLED) return 0;

				if (site == mettsStruct_.sites[i] && stage_[i]==DISABLED) stage_[i]=OPERATOR;
				else stage_[i]=WFT_NOADVANCE;
				if (stage_[i] == OPERATOR) checkOrder(i);

				if (timesWithoutAdvancement >= mettsStruct_.advanceEach) {
					stage_[i] = WFT_ADVANCE;
					if (i==lastI) {
						currentBeta_ += mettsStruct_.tau;
						timesWithoutAdvancement=0;
					}
				} else {
					if (i==lastI && stage_[i]==WFT_NOADVANCE) 
						timesWithoutAdvancement++;
				}

				std::ostringstream msg2;
				msg2<<"Steps without advance: "<<timesWithoutAdvancement;
				if (timesWithoutAdvancement>0) progress_.printline(msg2,std::cout);
				
				std::ostringstream msg;
				msg<<"Evolving, stage="<<getStage(i)<<" site="<<site<<" loopNumber="<<loopNumber;
				msg<<" Eg="<<Eg;
				progress_.printline(msg,std::cout);
				
				// phi = A|psi>
				computePhi(i,phiNew,phiOld,direction);
				
				return 1;
			}

			void initialGuess(VectorWithOffsetType& v) const
			{
				waveFunctionTransformation_.setInitialVector(v,psi_,lrs_);
				bool b = allStages(WFT_ADVANCE) || allStages(WFT_NOADVANCE);
				if (!b) return;
				std::vector<VectorWithOffsetType> vv(targetVectors_.size());
				for (size_t i=0;i<targetVectors_.size();i++) {
					waveFunctionTransformation_.setInitialVector(vv[i],
						targetVectors_[i],lrs_);
					if (norm(vv[i])<1e-6) continue;
					VectorWithOffsetType w= weight_[i]*vv[i];
					v += w;
				}
			}


			template<typename IoOutputType>
			void save(const std::vector<size_t>& block,IoOutputType& io) const
			{
				std::cerr<<"METTS: WARNING: save(...) unimplemented\n";
// 				std::ostringstream msg;
// 				msg<<"Saving state...";
// 				progress_.printline(msg,std::cout);
// 
// 				TimeSerializerType ts(currentBeta_,block[0],targetVectors_);
// 				ts.save(io);
// 				psi_.save(io,"PSI");
			}

		private:
			
			// in situ computation:
			void cocoon(size_t direction,const BlockType& block) const
			{
				size_t site = block[0];
				std::cerr<<"-------------&*&*&* In-situ measurements start\n";
				test(psi_,psi_,direction,"<PSI|A|PSI>",site);
				
				for (size_t j=0;j<targetVectors_.size();j++) {
					std::string s = "<P"+ttos(j)+"|A|P"+ttos(j)+">";
					test(targetVectors_[j],targetVectors_[j],direction,s,site);
				}
				std::cerr<<"-------------&*&*&* In-situ measurements end\n";
			}

			void checkOrder(size_t i) const
			{
				if (i==0) return;
				for (size_t j=0;j<i;j++) {
					if (stage_[j] == DISABLED) {
						std::string s ="TST:: Seeing tst site "+ttos(mettsStruct_.sites[i]);
						s =s + " before having seen";
						s = s + " site "+ttos(j);
						s = s +". Please order your tst sites in order of appearance.\n";
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
					case WFT_ADVANCE:
						return "WFT with time stepping";
						break;
					case WFT_NOADVANCE:
						return "WFT without time change";
						break;
				}
				return "undefined";
			}


			void computePhi(size_t i,VectorWithOffsetType& phiNew,
					const VectorWithOffsetType& phiOld,size_t systemOrEnviron)
			{
				size_t indexAdvance = betas_.size()-1;
				size_t indexNoAdvance = 0;
				if (stage_[i]==OPERATOR) {
					std::ostringstream msg;
					msg<<"I'm applying a local operator now";
					progress_.printline(msg,std::cout);
					FermionSign fs(lrs_.left(),mettsStruct_.electrons);
					applyOpLocal_(phiNew,phiOld,mettsStruct_.aOperators[i],fs,systemOrEnviron);
					RealType norma = norm(phiNew);
					if (norma==0) throw std::runtime_error("Norm of phi is zero\n");
					//std::cerr<<"Norm of phi="<<norma<<" when i="<<i<<"\n";


				} else if (stage_[i]== WFT_NOADVANCE || stage_[i]== WFT_ADVANCE) {
					size_t advance = indexNoAdvance;
					if (stage_[i] == WFT_ADVANCE) advance = indexAdvance;
					std::ostringstream msg;
					msg<<"I'm calling the WFT now";
					progress_.printline(msg,std::cout);
					
					if (mettsStruct_.aOperators.size()==1) guessPhiSectors(phiNew,i,systemOrEnviron);
					else phiNew.populateSectors(lrs_.super());
					
					// OK, now that we got the partition number right, let's wft:
					waveFunctionTransformation_.setInitialVector(phiNew,targetVectors_[advance],
							lrs_); // generalize for su(2)
					phiNew.collapseSectors();
					
				} else {
					throw std::runtime_error("It's 5 am, do you know what line "
						" your code is exec-ing?\n");
				}
			}		


			void calcTimeVectors(
						RealType Eg,
      						const VectorWithOffsetType& phi,
						size_t systemOrEnviron)
			{
				std::vector<ComplexMatrixType> V(phi.sectors());
				std::vector<ComplexMatrixType> T(phi.sectors());
				
				std::vector<size_t> steps(phi.sectors());
				
				triDiag(phi,T,V,steps);
				
				std::vector<std::vector<RealType> > eigs(phi.sectors());
						
				for (size_t ii=0;ii<phi.sectors();ii++) 
					PsimagLite::diag(T[ii],eigs[ii],'V');
				
				calcTargetVectors(phi,T,V,Eg,eigs,steps,systemOrEnviron);
			}

			void calcTargetVectors(
						const VectorWithOffsetType& phi,
						const std::vector<ComplexMatrixType>& T,
						const std::vector<ComplexMatrixType>& V,
						RealType Eg,
      						const std::vector<VectorType>& eigs,
	    					std::vector<size_t> steps,
					      	size_t systemOrEnviron)
			{
				targetVectors_[0] = phi;
				for (size_t i=1;i<betas_.size();i++) {
					// Only time differences here (i.e. betas_[i] not betas_[i]+currentBeta_)
					calcTargetVector(targetVectors_[i],phi,T,V,Eg,eigs,betas_[i],steps);
					//normalize(targetVectors_[i]);
				}
			}

			void calcTargetVector(
						VectorWithOffsetType& v,
      						const VectorWithOffsetType& phi,
						const std::vector<ComplexMatrixType>& T,
						const std::vector<ComplexMatrixType>& V,
						RealType Eg,
      						const std::vector<VectorType>& eigs,
	    					RealType t,
						std::vector<size_t> steps)
			{
				v = phi;
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i0 = phi.sector(ii);
					ComplexVectorType r;
					calcTargetVector(r,phi,T[ii],V[ii],Eg,eigs[ii],t,steps[ii],i0);
					v.setDataInSector(r,i0);
				}
			}

			void calcTargetVector(
						ComplexVectorType& r,
      						const VectorWithOffsetType& phi,
						const ComplexMatrixType& T,
						const ComplexMatrixType& V,
						RealType Eg,
      						const VectorType& eigs,
	    					RealType t,
	  					size_t steps,
						size_t i0)
			{
				size_t n2 = steps;
				size_t n = V.n_row();
				if (T.n_col()!=T.n_row()) throw std::runtime_error("T is not square\n");
				if (V.n_col()!=T.n_col()) throw std::runtime_error("V is not nxn2\n");
				// for (size_t j=0;j<v.size();j++) v[j] = 0; <-- harmful if v is sparse
				ComplexType zone = 1.0;
				ComplexType zzero = 0.0;
				
				ComplexVectorType tmp(n2);
				r.resize(n2);
				calcR(r,T,V,phi,Eg,eigs,t,steps,i0);
				psimag::BLAS::GEMV('N', n2, n2, zone, &(T(0,0)), n2, &(r[0]), 1, zzero, &(tmp[0]), 1 );
				r.resize(n);
				psimag::BLAS::GEMV('N', n,  n2, zone, &(V(0,0)), n, &(tmp[0]),1, zzero, &(r[0]),   1 );
			}

			void calcR(
				ComplexVectorType& r,
    				const ComplexMatrixType& T,
				const ComplexMatrixType& V,
    				const VectorWithOffsetType& phi,
    				RealType Eg,
				const VectorType& eigs,
    				RealType t,
				size_t n2,
				size_t i0)
			{
				for (size_t k=0;k<n2;k++) {
					ComplexType sum = 0.0;
					for (size_t kprime=0;kprime<n2;kprime++) {
						ComplexType tmpV = calcVTimesPhi(kprime,V,phi,i0);
						sum += conj(T(kprime,k))*tmpV;
					}
					RealType tmp = (eigs[k]-Eg)*t;
					ComplexType c(cos(tmp),sin(tmp));
					r[k] = c * sum;
				}
			}

			ComplexType calcVTimesPhi(size_t kprime,const ComplexMatrixType& V,const VectorWithOffsetType& phi,
						 size_t i0)
			{
				ComplexType ret = 0;
				size_t total = phi.effectiveSize(i0);
				
				for (size_t j=0;j<total;j++)
					ret += conj(V(j,kprime))*phi.fastAccess(i0,j);
				return ret;
			}

			void triDiag(
					const VectorWithOffsetType& phi,
					std::vector<ComplexMatrixType>& T,
	 				std::vector<ComplexMatrixType>& V,
					std::vector<size_t>& steps)
			{
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i = phi.sector(ii);
					steps[ii] = triDiag(phi,T[ii],V[ii],i);
				}
			}

			size_t triDiag(const VectorWithOffsetType& phi,ComplexMatrixType& T,ComplexMatrixType& V,size_t i0)
			{
				size_t p = lrs_.super().findPartitionNumber(phi.offset(i0));
				typename ModelType::ModelHelperType modelHelper(p,lrs_,model_.orbitals());
				 		//,useReflection_);
				typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);
			
				RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
				size_t iter= ProgramGlobals::LanczosSteps;

				//srand48(3243447);
				LanczosSolverType lanczosSolver(lanczosHelper,iter,eps,parallelRank_);
				
				TridiagonalMatrixType ab;
				size_t total = phi.effectiveSize(i0);
				TargetVectorType phi2(total);
				phi.extract(phi2,i0);
				/* std::ostringstream msg;
				msg<<"Calling tridiagonalDecomposition...\n";
				progress_.printline(msg,std::cerr);*/
				lanczosSolver.tridiagonalDecomposition(phi2,ab,V);
				ab.buildDenseMatrix(T);
				//check1(V,phi2);
				return lanczosSolver.steps();
			}

			//! This check is invalid if there are more than one sector
			void check1(const ComplexMatrixType& V,const TargetVectorType& phi2)
			{
				if (V.n_col()>V.n_row()) throw std::runtime_error("cols > rows\n");
				TargetVectorType r(V.n_col());
				for (size_t k=0;k<V.n_col();k++) {
					r[k] = 0.0;
					for (size_t j=0;j<V.n_row();j++) 
						r[k] += conj(V(j,k))*phi2[j];
					// is r(k) == \delta(k,0)
					if (k==0 && std::norm(r[k]-1.0)>1e-5) 
						std::cerr<<"WARNING: r[0]="<<r[0]<<" != 1\n";
					if (k>0 && std::norm(r[k])>1e-5) 
						std::cerr<<"WARNING: r["<<k<<"]="<<r[k]<<" !=0\n";
				}
			}

			void guessPhiSectors(VectorWithOffsetType& phi,size_t i,size_t systemOrEnviron)
			{
				FermionSign fs(lrs_.left(),mettsStruct_.electrons);
				if (allStages(WFT_NOADVANCE)) {
					VectorWithOffsetType tmpVector = psi_;
					for (size_t j=0;j<mettsStruct_.aOperators.size();j++) {
						applyOpLocal_(phi,tmpVector,mettsStruct_.aOperators[j],fs,
							systemOrEnviron);
						tmpVector = phi;
					}
					return;
				}
				applyOpLocal_(phi,psi_,mettsStruct_.aOperators[i],fs,
								systemOrEnviron);
			}

			void zeroOutVectors()
			{
				for (size_t i=0;i<targetVectors_.size();i++) 
					targetVectors_[i].resize(lrs_.super().size());
			}

//			void printHeader()
//			{
//				io_.print(mettsStruct_);
//				std::string label = "times";
//				io_.printVector(betas_,label);
//				label = "weights";
//				io_.printVector(weight_,label);
//				std::string s = "GsWeight="+ttos(gsWeight_);
//				io_.printline(s);
//			}

			void test(	
					const VectorWithOffsetType& src1,
					const VectorWithOffsetType& src2,
					size_t systemOrEnviron,
				 	const std::string& label,
					size_t site) const
			{
				VectorWithOffsetType dest;
				OperatorType A = mettsStruct_.aOperators[0];
				PsimagLite::CrsMatrix<ComplexType> tmpC(model_.getOperator("c",0,0));
				PsimagLite::CrsMatrix<ComplexType> tmpCt;
				transposeConjugate(tmpCt,tmpC);
				multiply(A.data,tmpCt,tmpC);
				A.fermionSign = 1;
				//A.data = tmpC;
				FermionSign fs(lrs_.left(),mettsStruct_.electrons);
				applyOpLocal_(dest,src1,A,fs,systemOrEnviron);

				ComplexType sum = 0;
				for (size_t ii=0;ii<dest.sectors();ii++) {
					size_t i = dest.sector(ii);
					size_t offset1 = dest.offset(i);
					for (size_t jj=0;jj<src2.sectors();jj++) {
						size_t j = src2.sector(jj);
						size_t offset2 = src2.offset(j);
						if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
						for (size_t k=0;k<dest.effectiveSize(i);k++) 
							sum+= dest[k+offset1] * conj(src2[k+offset2]);
					}
				}
				std::cerr<<site<<" "<<sum<<" "<<" "<<currentBeta_;
				std::cerr<<" "<<label<<std::norm(src1)<<" "<<std::norm(src2)<<" "<<std::norm(dest)<<"\n";
			}


			std::vector<size_t> stage_;
			VectorWithOffsetType psi_;
			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& mettsStruct_;
			const WaveFunctionTransfType& waveFunctionTransformation_;
			PsimagLite::ProgressIndicator progress_;
			RealType currentBeta_;
			std::vector<RealType> betas_,weight_;
			std::vector<VectorWithOffsetType> targetVectors_;
			RealType gsWeight_;
			//typename IoType::Out io_;
			ApplyOperatorType applyOpLocal_;

	};     //class MettsTargetting

	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
 		template<typename,typename> class WaveFunctionTransfTemplate,
			typename ModelType_,
 		typename ConcurrencyType_,
			typename IoType_,
   			template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const MettsTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,ConcurrencyType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		tst.print(os);
		return os;
	}
} // namespace Dmrg

#endif //METTS_TARGETTING_H
