// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

/*! \file TimeStepTargetting.h
 *
 *  see http://de.arxiv.org/abs/cond-mat/0606018v1
 *  Eq. (55) to (60)  page 26.
 *
 */
 
#ifndef TIMESTEPTARGETTING_H
#define TIMESTEPTARGETTING_H
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TimeStepStructure.h"
#include "ApplyOperatorLocal.h"

namespace Dmrg {
	template<
			template<typename,typename> class LanczosSolverTemplate,
   			template<typename,typename> class InternalProductTemplate,
	 		typename WaveFunctionTransformationType_,
    			typename ModelType_,
	 		typename ConcurrencyType_,
    			typename IoType_,
       			template<typename> class VectorWithOffsetTemplate>
	class TimeStepTargetting  {
		public:
			typedef WaveFunctionTransformationType_ WaveFunctionTransformationType;
			typedef ModelType_ ModelType;
			typedef ConcurrencyType_ ConcurrencyType;
			typedef IoType_ IoType;
			typedef typename ModelType::RealType RealType;
			typedef std::complex<RealType> ComplexType;
			typedef InternalProductTemplate<ComplexType,ModelType> InternalProductType;
			typedef typename ModelType::OperatorsType OperatorsType;
			typedef typename ModelType::MyBasisWithOperators BasisWithOperatorsType;
			//typedef BasisWithOperators<OperatorsType,ConcurrencyType> BasisWithOperatorsType;
			typedef std::vector<ComplexType> ComplexVectorType;
			//typedef std::VectorWithOffset<ComplexType> VectorWithOffsetType;
			typedef LanczosSolverTemplate<InternalProductType,ComplexVectorType> LanczosSolverType;
			typedef std::vector<RealType> VectorType;
			//typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
			typedef psimag::Matrix<ComplexType> ComplexMatrixType;
			typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef TimeStepStructure<OperatorType> TargettingStructureType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<ComplexType> VectorWithOffsetType;
			typedef ComplexVectorType TargetVectorType;
			typedef BlockMatrix<ComplexType,ComplexMatrixType> ComplexBlockMatrixType;
			typedef ApplyOperatorLocal<BasisWithOperatorsType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
			
			enum {DISABLED,OPERATOR,WFT_NOADVANCE,WFT_ADVANCE};
			enum {SHRINK_SYSTEM=WaveFunctionTransformationType::SHRINK_SYSTEM,
			SHRINK_ENVIRON=WaveFunctionTransformationType::SHRINK_ENVIRON,
			INFINITE=WaveFunctionTransformationType::INFINITE};

			static const size_t parallelRank_ = 0; // TST needs to support concurrency FIXME

			//enum {INDEX_NOADVANCE=0,INDEX_ADVANCE=1};

			TimeStepTargetting(
	  				const BasisWithOperatorsType& basisS,
       					const BasisWithOperatorsType& basisE,
	    				const BasisType& basisSE,
	 				const ModelType& model,
					const TargettingStructureType& tstStruct,
					const WaveFunctionTransformationType& wft) 
				: stage_(tstStruct.sites.size(),DISABLED),basisS_(basisS),basisE_(basisE),basisSE_(basisSE),
					model_(model),tstStruct_(tstStruct),waveFunctionTransformation_(wft),
					progress_("TimeStepTargetting",0),currentTime_(0),
							times_(tstStruct_.timeSteps),weight_(tstStruct_.timeSteps),targetVectors_(tstStruct_.timeSteps),
						io_(tstStruct_.filename,parallelRank_),applyOpLocal_(basisS,basisE,basisSE,SHRINK_ENVIRON)
			{
				if (!wft.isEnabled()) throw std::runtime_error(" TimeStepTargetting "
							"needs an enabled wft\n");
				
				RealType tau =tstStruct_.tau;
				RealType sum = 0;
				size_t n = times_.size();
				for (size_t i=0;i<n;i++) {
					times_[i] = i*tau/(n-1);
					weight_[i] = 1.0/(n+2);
					sum += weight_[i];
				}
				sum -= weight_[0];
				weight_[0] = weight_[n-1] = 2.0/(n+2);
				sum += weight_[0];
				
				gsWeight_=1.0-sum;
				sum += gsWeight_;
				//for (size_t i=0;i<weight_.size();i++) sum += weight_[i];
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
				//return 1.0;
			}
			
			RealType normSquared(size_t i) const
			{
				return real(multiply(targetVectors_[i],targetVectors_[i])); // call to mult will conj()
			}

			template<typename SomeBasisType>
			void setGs(const std::vector<TargetVectorType>& v,//const std::vector<size_t>& weights,
				   const SomeBasisType& someBasis)
			{
				psi_.set(v,someBasis);
			}

			const ComplexType& operator[](size_t i) const { return psi_[i]; }

			ComplexType& operator[](size_t i) { return psi_[i]; }

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

			void evolve(RealType Eg,size_t direction,const BlockType& block,size_t loopNumber,
				   bool needsPrinting)
			{
				size_t count =0;
				VectorWithOffsetType phiOld = psi_;
				VectorWithOffsetType phiNew;
				size_t max = tstStruct_.sites.size();
				
				if (noStageIs(DISABLED)) max = 1;
				
				for (size_t i=0;i<max;i++) {
					count += evolve(i,phiNew,phiOld,Eg,direction,block,loopNumber,max-1);
					phiOld = phiNew;
				}
				
				if (count==0) {
					// always print to keep observer driver sync
					if (needsPrinting) {
						zeroOutVectors();
						printVectors();
					}
					return;
				}
				
				calcTimeVectors(Eg,phiNew,direction);
				
				cocoon(direction,block); // in-situ
				
				if (needsPrinting) printVectors(); // for post-processing
			}

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
				
				static size_t  timesWithoutAdvancement=0;
				if (tstStruct_.startingLoops[i]>loopNumber || direction==INFINITE) return 0;
				if (block.size()!=1) throw 
					std::runtime_error("TimeStepTargetting::evolve(...):"
							" blocks of size != 1 are unsupported (sorry)\n");
				size_t site = block[0];
				//size_t max = tstStruct_.aOperators.size()-1;
				if (site != tstStruct_.sites[i] && stage_[i]==DISABLED) return 0;
				
				
				if (site == tstStruct_.sites[i] && stage_[i]==DISABLED) stage_[i]=OPERATOR;
				else stage_[i]=WFT_NOADVANCE;
				
				if (timesWithoutAdvancement >= tstStruct_.advanceEach) {
					stage_[i] = WFT_ADVANCE;
					if (i==lastI) {
						currentTime_ += tstStruct_.tau;
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

			const BasisType& basisSE() const { return basisSE_; }

			const BasisWithOperatorsType& basisS() const { return basisS_; }

			const BasisWithOperatorsType& basisE() const { return basisE_; }

			void initialGuess(VectorWithOffsetType& v) const
			{
				waveFunctionTransformation_.createRandomVector(v);
			}

		private:
			
			// in situ computation:
			void cocoon(size_t direction,const BlockType& block) const
			{
				size_t site = block[0];
				std::cerr<<"-------------&*&*&* Cocoon output starts\n";
				test(psi_,direction,"<PSI|A|PSI>",site);
				
				for (size_t j=0;j<targetVectors_.size();j++) {
					std::string s = "<P"+utils::ttos(j)+"|A|P"+utils::ttos(j)+">";
					test(targetVectors_[j],direction,s,site);
				}
				std::cerr<<"-------------&*&*&* Cocoon output ends\n";
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
					VectorWithOffsetType& phiOld,size_t systemOrEnviron)
			{
				size_t indexAdvance = times_.size()-1;
				size_t indexNoAdvance = 0;
				if (stage_[i]==OPERATOR) {
					std::ostringstream msg;
					msg<<"I'm applying a local operator now";
					progress_.printline(msg,std::cout);
					FermionSign fs(basisS_,tstStruct_.electrons);
					applyOpLocal_(phiNew,phiOld,tstStruct_.aOperators[i],fs,systemOrEnviron);
					RealType norma = norm(phiNew);
					if (norma==0) throw std::runtime_error("Norm of phi is zero\n");
					std::cerr<<"Norm of phi="<<norma<<" when i="<<i<<"\n";
				} else if (stage_[i]== WFT_NOADVANCE || stage_[i]== WFT_ADVANCE) {
					size_t advance = indexNoAdvance;
					if (stage_[i] == WFT_ADVANCE) advance = indexAdvance;
					std::ostringstream msg;
					msg<<"I'm calling the WFT now";
					progress_.printline(msg,std::cout);
					
					if (tstStruct_.aOperators.size()==1) guessPhiSectors(phiNew,i,systemOrEnviron);
					else phiNew.populateSectors(basisSE_);
					
					// OK, now that we got the partition number right, let's wft:
					waveFunctionTransformation_.setInitialVector(phiNew,targetVectors_[advance],
							basisS_,basisE_,basisSE_); // generalize for su(2)
					phiNew.collapseSectors();
					
				} else {
					throw std::runtime_error("It's 5 am, do you know what line your code is exec-ing?\n");
				}
				
				//normalize(phi);
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
					utils::diag(T[ii],eigs[ii],'V');
				
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
				for (size_t i=1;i<times_.size();i++) {
					// Only time differences here (i.e. times_[i] not times_[i]+currentTime_)
					calcTargetVector(targetVectors_[i],phi,T,V,Eg,eigs,times_[i],steps);
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
				
				/*std::cerr<<"Trying T\n";
				RealType eps = 1e-5;
				if (!isUnitary(T,eps)) throw std::runtime_error("Oops: T\n");;
				std::cerr<<"Passed T\n";*/
				//RealType norma = std::norm(r);
				//std::cerr<<"r norma="<<norma<<"\n";
				//for (size_t i=0;i<n;i++) 
				//	std::cout<<r[i]<<" "<<phi[i+phi.offset(i0)]<<"\n";
				
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
						//if (k==0) std::cerr<<"HERE:"<<tmpV<<" "<<V(0,kprime)<<"\n";
						sum += conj(T(kprime,k))*tmpV;
					}
					//if (t!=0) throw std::runtime_error("time is not zero\n");
					RealType tmp = (eigs[k]-Eg)*t;
					ComplexType c(cos(tmp),sin(tmp));
					r[k] = c * sum;
				}
			}

			ComplexType calcVTimesPhi(size_t kprime,const ComplexMatrixType& V,const VectorWithOffsetType& phi,
						 size_t i0)
			{
				ComplexType ret = 0;
				//for (size_t ii=0;ii<phi.sectors();ii++) {
					//size_t i = phi.sector(ii);
					size_t total = phi.effectiveSize(i0);
					//size_t offset = phi.offset(i);
				
					for (size_t j=0;j<total;j++)
						ret += conj(V(j,kprime))*phi.fastAccess(i0,j);
				//}
				
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
				size_t p = basisSE_.findPartitionNumber(phi.offset(i0));
				typename ModelType::ModelHelperType modelHelper(p,basisSE_,basisS_,basisE_,model_.orbitals());
				 		//,useReflection_);
				typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);
			
				size_t mode = LanczosSolverType::WITH_INFO;
				RealType eps= 0.01*ProgramLimits::LanczosTolerance;
				size_t iter= ProgramLimits::LanczosSteps;

				//srand48(3243447);
				LanczosSolverType lanczosSolver(lanczosHelper,iter,eps,parallelRank_,mode);
				
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

			// hack to get the current partition, note that this:
			// size_t partition = targetVectors_[0].findPartition(basisSE_);
			// doesn't work, since targetVectors_[0] is stale at this point
			// This hack could cause performance problems, must be checked
			// any other ideas?
			// Update: When more than one op. needs to apply all of 'em:
			void guessPhiSectors(VectorWithOffsetType& phi,size_t i,size_t systemOrEnviron)
			{
				FermionSign fs(basisS_,tstStruct_.electrons);
				if (allStages(WFT_NOADVANCE)) {
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

			void printVectors()
			{
				for (size_t i=0;i<targetVectors_.size();i++) {
					std::string label = "targetVector"+utils::ttos(i)+"_"+utils::ttos(currentTime_);
					targetVectors_[i].save(io_,label);
				}
			}

			void printHeader()
			{
				io_.print(tstStruct_);
				std::string label = "times";
				io_.printVector(times_,label);
				label = "weights";
				io_.printVector(weight_,label);
				std::string s = "GsWeight="+utils::ttos(gsWeight_);
				io_.printline(s);
			}

			void test(	
					const VectorWithOffsetType& src,
					size_t systemOrEnviron,
				 	const std::string& label,
					size_t site) const
			{
				VectorWithOffsetType dest;
				OperatorType A = tstStruct_.aOperators[0];
				CrsMatrix<ComplexType> tmpC(model_.getOperator("c",0,0));
				CrsMatrix<ComplexType> tmpCt;
				transposeConjugate(tmpCt,tmpC);
				multiply(A.data,tmpCt,tmpC);
				A.fermionSign = 1;
				//A.data = tmpC;
				FermionSign fs(basisS_,tstStruct_.electrons);
				applyOpLocal_(dest,src,A,fs,systemOrEnviron);
				
				
				ComplexType sum = 0;
				for (size_t ii=0;ii<dest.sectors();ii++) {
					size_t i = dest.sector(ii);
					for (size_t jj=0;jj<dest.sectors();jj++) {
						size_t j = src.sector(jj);
						if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
						size_t offset = dest.offset(i);
						for (size_t k=0;k<dest.effectiveSize(i);k++) 
							sum+= dest[k+offset] * conj(src[k+offset]);
					}
				}
				std::cerr<<site<<" "<<sum<<" "<<" "<<currentTime_<<" "<<label<<std::norm(src)<<" "<<std::norm(dest)<<"\n";
			}

			std::vector<size_t> stage_;
			VectorWithOffsetType psi_;
			const BasisWithOperatorsType& basisS_;
			const BasisWithOperatorsType& basisE_;
			const BasisType& basisSE_;
			const ModelType& model_;
			const TargettingStructureType& tstStruct_;
			const WaveFunctionTransformationType& waveFunctionTransformation_;
			ProgressIndicator progress_;
			RealType currentTime_;
			std::vector<RealType> times_,weight_;
			std::vector<VectorWithOffsetType> targetVectors_;
			RealType gsWeight_;
			typename IoType::Out io_;
			ApplyOperatorType applyOpLocal_;
	};     //class TimeStepTargetting
} // namespace Dmrg
/*@}*/
#endif
