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
#include "SimpleReader.h"
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TargetStructureParams.h"

namespace Dmrg {
	
	template<typename OperatorType>
	struct TimeStepStructure {
		std::string filename;
		typename OperatorType::RealType tau;
		size_t timeSteps;
		size_t advanceEach;
		std::vector<size_t> sites;
		std::vector<size_t> startingLoops;
		std::vector<OperatorType> aOperators;
		std::vector<size_t> electrons;
	};
	
	/*template<typename OperatorType>
	inline TimeStepStructure<OperatorType>&
	operator<=(TimeStepStructure<OperatorType>& t,SimpleReader& reader)
	{
		reader.read(t.aOperator);
		reader.read(t.electrons);
		reader.read(t.site);
		reader.read(t.startingLoop);
		reader.read(t.filename);
		return t;
	}*/
	
	template<typename OperatorType,typename ModelType>
	inline TargetStructureParams<TimeStepStructure<OperatorType>,ModelType>&
	operator<=(TargetStructureParams<TimeStepStructure<OperatorType>,ModelType>& tsp,SimpleReader& reader)
	{
		typedef typename ModelType::RealType RealType;
		std::vector<size_t> sites,loops;
		std::string s;
		reader.read(s); // filename
		RealType tau=0;
		reader.read(tau);
		size_t timeSteps=0;
		reader.read(timeSteps);
		size_t advanceEach=0;
		reader.read(advanceEach);
		reader.read(sites);
		reader.read(loops);
		
		tsp.init(s,tau,timeSteps,advanceEach,sites,loops);
		
		for (size_t i=0;i<sites.size();i++) {
			//std::string s;
			reader.read(s);
			if (s == "cooked") {
				reader.read(s);
				std::vector<size_t> v;
				reader.read(v);
				tsp.setCookedData(i,s,v);
			} else {
				psimag::Matrix<RealType> m;
				reader.read(m);
				tsp.setRawData(i,m);
			}
			int fermiSign=0;
			reader.read(fermiSign);
			std::pair<size_t,size_t> jmValues;
			reader.read(jmValues);
			RealType angularFactor;
			reader.read(angularFactor);
			tsp.set(i,fermiSign,jmValues,angularFactor);
		}
		
		return tsp;
	}
	
	template<typename OperatorType>
	inline std::ostream&
	operator<<(std::ostream& os,const TimeStepStructure<OperatorType>& t)
	{
		for (size_t i=0;i<t.aOperators.size();i++)
			os<<t.aOperators[i];
		os<<t.electrons;
		os<<"TimeStepStructure.site="<<t.sites<<"\n";
		os<<"TimeStepStructure.startingLoop="<<t.startingLoops<<"\n";
		os<<"TimeStepStructure.filename="<<t.filename<<"\n";
		os<<"TimeVectorsfilename.tau="<<t.tau<<"\n";
		os<<"TimeVectorsfilename.timeSteps="<<t.timeSteps<<"\n";
		os<<"TimeVectorsfilename.advanceEach="<<t.advanceEach<<"\n";
		return os;
	}
	
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
			
			enum {DISABLED,OPERATOR,WFT_NOADVANCE,WFT_ADVANCE};
			enum {SHRINK_SYSTEM=WaveFunctionTransformationType::SHRINK_SYSTEM,
			SHRINK_ENVIRON=WaveFunctionTransformationType::SHRINK_ENVIRON,
			INFINITE=WaveFunctionTransformationType::INFINITE};
			
			static const size_t parallelRank_ = 0; // TST needs to support concurrency FIXME
			
			enum {INDEX_NOADVANCE=0,INDEX_ADVANCE=1};
			
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
						io_(tstStruct_.filename,parallelRank_)
			{
				if (!wft.isEnabled()) throw std::runtime_error(" TimeStepTargetting "
							"needs an enabled wft\n");
				
				RealType tau =tstStruct_.tau;
				RealType sum = 0;
				times_[INDEX_NOADVANCE]=0;
				weight_[INDEX_NOADVANCE]=1.0/3.0;
				sum += weight_[INDEX_NOADVANCE];
				if (times_.size()>INDEX_ADVANCE) {
					times_[INDEX_ADVANCE]=tau;
					weight_[INDEX_ADVANCE]=1.0/3.0;
					sum += weight_[INDEX_ADVANCE];
				}
				if (times_.size()>2) {
					times_[2]=tau/3.;
					weight_[2]=1.5/12.0;
					sum += weight_[2];
				}
				
				if (times_.size()>3) {
					times_[3]=2.0*tau/3.;
					weight_[3]=1.5/12.0;
					sum += weight_[3];
				}
				
				
				for (size_t i=0;i<times_.size();i++) {
					times_[i] = i*tau/times_.size();
					weight_[i] = 1.0/(times_.size()+2);
				}
				weight_[0] = weight_[9] = 2.0/(times_.size()+2);
				sum = 1.0;
				
				gsWeight_=1.0-sum;
				sum = gsWeight_;
				for (size_t i=0;i<weight_.size();i++) sum += weight_[i];
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

			void evolve(RealType Eg,size_t direction,const BlockType& block,size_t loopNumber)
			{
				for (size_t i=0;i<tstStruct_.sites.size();i++) 
					evolve(i,Eg,direction,block,loopNumber);	
			}

			void evolve(size_t i,RealType Eg,size_t direction,const BlockType& block,size_t loopNumber)
			{
				static size_t  timesWithoutAdvancement=0;
				if (tstStruct_.startingLoops[i]>loopNumber || direction==INFINITE) return;
				if (block.size()!=1) throw 
					std::runtime_error("TimeStepTargetting::evolve(...):"
							" blocks of size != 1 are unsupported (sorry)\n");
				size_t site = block[0];
				size_t max = tstStruct_.aOperators.size()-1;
				if (site != tstStruct_.sites[i] && stage_[i]==DISABLED) return;
				
				
				if (site == tstStruct_.sites[i] && stage_[i]==DISABLED) stage_[i]=OPERATOR;
				else stage_[i]=WFT_NOADVANCE;
				
				if (timesWithoutAdvancement >= tstStruct_.advanceEach) {
					stage_[i] = WFT_ADVANCE;
					if (i==max) {
						currentTime_ += times_[INDEX_ADVANCE];
						timesWithoutAdvancement=0;
					}
				} else {
					if (i==max && stage_[i]==WFT_NOADVANCE) 
						timesWithoutAdvancement++;
				}
				
				std::ostringstream msg;
				msg<<"Evolving, stage="<<getStage(i)<<" site="<<site<<" loopNumber="<<loopNumber;
				progress_.printline(msg,std::cout);
				
				VectorWithOffsetType phi; // phi = A|psi>
				computePhi(i,phi,direction);
				if (norm(phi)==0) throw std::runtime_error("Norm of phi is zero\n");
				test(i,psi_,direction,"<PSI|A|PSI>",site);
				test(i,phi,direction,"<PHI|A|PHI>",site);

				calcTimeVectors(Eg,phi);
				test(i,targetVectors_[INDEX_ADVANCE],direction,"<PTAU|A|PTAU>",site);
				printVectors();
			}

			const BasisType& basisSE() const { return basisSE_; }

			const BasisWithOperatorsType& basisS() const { return basisS_; }

			const BasisWithOperatorsType& basisE() const { return basisE_; }

			void initialGuess(VectorWithOffsetType& v) const
			{
				waveFunctionTransformation_.createRandomVector(v);
			}

		private:
			
			bool allStages(size_t x) const
			{
				for (size_t i=0;i<stage_.size();i++)
					if (stage_[i]!=x) return false;
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
			
			void computePhi(size_t i,VectorWithOffsetType& phi,size_t systemOrEnviron)
			{
				if (stage_[i]==OPERATOR) {
					std::ostringstream msg;
					msg<<"Applying local operator...";
					progress_.printline(msg,std::cout);
					applyLocalOp(phi,psi_,tstStruct_.aOperators[i],tstStruct_.electrons,systemOrEnviron);
				} else if (stage_[i]== WFT_NOADVANCE || stage_[i]== WFT_ADVANCE) {
					size_t advance = INDEX_NOADVANCE;
					if (stage_[i] == WFT_ADVANCE) advance = INDEX_ADVANCE;
					std::ostringstream msg;
					msg<<"Applying wft...";
					progress_.printline(msg,std::cout);
					
					if (tstStruct_.aOperators.size()==1) guessPhiSectors(phi,i,systemOrEnviron);
					else phi.populateSectors(basisSE_);
					
					// OK, now that we got the partition number right, let's wft:
					waveFunctionTransformation_.setInitialVector(phi,targetVectors_[advance],
							basisS_,basisE_,basisSE_); // generalize for su(2)
					phi.collapseSectors();
					
				} else {
					throw std::runtime_error("It's 5 am, do you know what line your code is exec-ing?\n");
				}
				
				normalize(phi);
			}
			
			void calcTimeVectors(RealType Eg,const VectorWithOffsetType& phi)
			{
				std::vector<ComplexMatrixType> V(phi.sectors());
				std::vector<ComplexMatrixType> T(phi.sectors());
				
				std::vector<size_t> steps(phi.sectors());
				
				triDiag(phi,T,V,steps);
				
				std::vector<std::vector<RealType> > eigs(phi.sectors());
						
				for (size_t ii=0;ii<phi.sectors();ii++) 
					utils::diag(T[ii],eigs[ii],'V');
				
				calcTargetVectors(phi,T,V,Eg,eigs,steps);
			}
			
			void calcTargetVectors(
						const VectorWithOffsetType& phi,
						const std::vector<ComplexMatrixType>& T,
						const std::vector<ComplexMatrixType>& V,
						RealType Eg,
      						const std::vector<VectorType>& eigs,
	    					std::vector<size_t> steps)
			{
				targetVectors_[0] = phi;
				for (size_t i=1;i<times_.size();i++) {
					//targetVectors_[i] = targetVectors_[0]; <-- ONLY FOR TESTING!!
					//! Only time differences here (i.e. times_[i] not times_[i]+currentTime_) 
					calcTargetVector(targetVectors_[i],phi,T,V,Eg,eigs,times_[i],steps);
					normalize(targetVectors_[i]);
					RealType norma = std::norm(targetVectors_[i]);
					std::cerr<<"norma after normalize="<<norma<<"\n";
					if (fabs(norma-1.0)>1e-5) throw std::runtime_error("Norm is not 1\n");
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
				calcR(r,T,V,phi,Eg,eigs,t,steps);
				psimag::BLAS::GEMV('C', n2, n2, zone, &(T(0,0)), n2, &(r[0]), 1, zzero, &(tmp[0]), 1 );
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
				size_t n2)
			{
				for (size_t k=0;k<n2;k++) {
					ComplexType sum = 0.0;
					for (size_t kprime=0;kprime<n2;kprime++) {
						ComplexType tmpV = calcVTimesPhi(kprime,V,phi);
						//if (k==0) std::cerr<<"HERE:"<<tmpV<<" "<<V(0,kprime)<<"\n";
						sum += T(k,kprime)*tmpV;
					}
					//if (t!=0) throw std::runtime_error("time is not zero\n");
					RealType tmp = (eigs[k]-Eg)*t;
					ComplexType c(cos(tmp),sin(tmp));
					r[k] = c * sum;
				}
			}

			ComplexType calcVTimesPhi(size_t kprime,const ComplexMatrixType& V,const VectorWithOffsetType& phi)
			{
				ComplexType ret = 0;
				for (size_t ii=0;ii<phi.sectors();ii++) {
					size_t i = phi.sector(ii);
					size_t total = phi.effectiveSize(i);
					//size_t offset = phi.offset(i);
				
					for (size_t j=0;j<total;j++)
						ret += conj(V(j,kprime))*phi.fastAccess(i,j);
				}
				
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
				check1(V,phi2);
				return lanczosSolver.steps();
			}

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

			//! FIXME: we need to make a fast version for when we're just
			//! figuring out where the (non-zero) partition is
			void applyLocalOp(
					VectorWithOffsetType& dest,
					const VectorWithOffsetType& src,
				     	const OperatorType& A,
	  				const std::vector<size_t>& electrons,
					size_t systemOrEnviron)
			{
				if (systemOrEnviron == SHRINK_ENVIRON) applyLocalOpSystem(dest,src,A,electrons);
				else applyLocalOpEnviron(dest,src,A);
			}

			void applyLocalOpSystem(
					VectorWithOffsetType& dest,
					const VectorWithOffsetType& src,
					const OperatorType& A,
	  				const std::vector<size_t>& electrons)
			{
				TargetVectorType dest2(basisSE_.size());
				
				for (size_t i=0;i<dest2.size();i++) dest2[i] = 0;
				
				for (size_t ii=0;ii<src.sectors();ii++) {
					size_t i = src.sector(ii);
					applyLocalOpSystem(dest2,src,A,electrons,i);
				}
				dest.fromFull(dest2,basisSE_);
			}

			void applyLocalOpSystem(
					TargetVectorType& dest2,
					const VectorWithOffsetType& src,
					const OperatorType& A,
	  				const std::vector<size_t>& electrons,
					size_t i0)
			{
				size_t offset = src.offset(i0);
				size_t final = offset + src.effectiveSize(i0);
				//size_t counter=0;
				size_t ns = basisS_.size();
				size_t nx = basisS_.size()/A.data.rank();
				
				for (size_t i=offset;i<final;i++) {
					size_t x=0,y=0;
					utils::getCoordinates(x,y,basisSE_.permutation(i),ns);
					size_t x0=0,x1=0;
					utils::getCoordinates(x0,x1,basisS_.permutation(x),nx);
					int nx0 = basisS_.electrons(x)-electrons[x1];
					if (nx0<0) throw std::runtime_error("TimeStepTargetting::applyLocalOpSystem(...)\n");
					RealType sign = ((nx0%2)==0) ? 1 : A.fermionSign;
					for (int k=A.data.getRowPtr(x1);k<A.data.getRowPtr(x1+1);k++) {
						size_t x1prime = A.data.getCol(k);
						size_t xprime = basisS_.permutationInverse(x0+x1prime*nx);
						size_t j = basisSE_.permutationInverse(xprime+y*ns);
						dest2[j] += src[i]*A.data.getValue(k)*sign;
					}
				}
				
			}

			void applyLocalOpEnviron(
					VectorWithOffsetType& dest,
					const VectorWithOffsetType& src,
					const OperatorType& A)
			{
				TargetVectorType dest2(basisSE_.size());
				
				for (size_t i=0;i<dest2.size();i++) dest2[i] = 0;
				
				for (size_t ii=0;ii<src.sectors();ii++) {
					size_t i = src.sector(ii);
					applyLocalOpEnviron(dest2,src,A,i);
				}
				dest.fromFull(dest2,basisSE_);
			}

			void applyLocalOpEnviron(
					TargetVectorType& dest2,
					const VectorWithOffsetType& src,
				     	const OperatorType& A,
					size_t i0)
			{
				size_t offset = src.offset(i0);
				size_t final = offset + src.effectiveSize(i0);
				
				size_t ns = basisS_.size();
				size_t nx = A.data.rank();
				
				for (size_t i=offset;i<final;i++) {
					size_t x=0,y=0;
					utils::getCoordinates(x,y,basisSE_.permutation(i),ns);
					size_t y0=0,y1=0;
					utils::getCoordinates(y0,y1,basisE_.permutation(y),nx);
					RealType sign = basisS_.fermionicSign(x,A.fermionSign);
					for (int k=A.data.getRowPtr(y0);k<A.data.getRowPtr(y0+1);k++) {
						size_t y0prime = A.data.getCol(k);
						size_t yprime = basisE_.permutationInverse(y0prime+y1*nx);
						size_t j = basisSE_.permutationInverse(x+yprime*ns);
						dest2[j] += src[i]*A.data.getValue(k)*sign;
					}
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
				if (allStages(WFT_NOADVANCE)) {
					VectorWithOffsetType tmpVector = psi_;
					for (size_t j=0;j<tstStruct_.aOperators.size();j++) {
						applyLocalOp(phi,tmpVector,tstStruct_.aOperators[j],tstStruct_.electrons,
							systemOrEnviron);
						tmpVector = phi;
					}
					return;
				}
				applyLocalOp(phi,psi_,tstStruct_.aOperators[i],tstStruct_.electrons,
								systemOrEnviron);
			}

			void printVectors()
			{
				for (size_t i=0;i<targetVectors_.size();i++) {
					std::string label = "targetVector"+utils::ttos(i)+"_"+utils::ttos(currentTime_);
					io_.printSparseVector(targetVectors_[i],label);
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

			void test(	size_t ind,
					const VectorWithOffsetType& src,
					size_t systemOrEnviron,
				 	const std::string& label,
					size_t site)
			{
				VectorWithOffsetType dest;
				OperatorType A = tstStruct_.aOperators[ind];
				CrsMatrix<ComplexType> tmpC(model_.getOperator("c",0,0));
				CrsMatrix<ComplexType> tmpCt;
				transposeConjugate(tmpCt,tmpC);
				multiply(A.data,tmpCt,tmpC);
				A.fermionSign = 1;
				//A.data = tmpC;
				applyLocalOp(dest,src,A,tstStruct_.electrons,systemOrEnviron);
				
				
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
	};     //class TimeStepTargetting
} // namespace Dmrg
/*@}*/
#endif
