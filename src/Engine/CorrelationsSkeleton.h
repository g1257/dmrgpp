// BEGIN LICENSE BLOCK
/*
Copyright (c)  2008 , UT-Battelle, LLC
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

/*! \file CorrelationsSkeleton.h
 *
 *  Helper class for Observables
 *
 */
#ifndef CORRELATIONS_SK_H
#define CORRELATIONS_SK_H
#include "ObserverHelper.h"
#include "CrsMatrix.h"
#include "Profiling.h"
#include "ApplyOperatorLocal.h"

namespace Dmrg {
	
	//! Don't add functions to this class
	template<typename FieldType>
	struct CorrelationData {
		size_t ni;
		size_t nj;
		std::vector<PsimagLite::Matrix<FieldType> > correlationVector;
		SparseVector<FieldType> wavefunction;
	};
	
	//! Companion function:
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,CorrelationData<FieldType>& c)
	{
		os<<"ni="<<c.ni<<"\n";
		os<<"nj="<<c.nj<<"\n";
		os<<c.correlationVector;
		os<<c.wavefunction;
		return os;
	}
	
	template<typename ObserverHelperType_>
	class CorrelationsSkeleton {
		//typedef typename MatrixType::value_type FieldType;
		typedef size_t IndexType;

	public:
		typedef ObserverHelperType_ ObserverHelperType;
		//typedef typename ObserverHelperType::IoInputType IoInputType;
		typedef typename ObserverHelperType::MatrixType MatrixType;
		typedef typename ObserverHelperType::VectorType VectorType ;
		typedef typename ObserverHelperType::VectorWithOffsetType VectorWithOffsetType;
		typedef typename ObserverHelperType::BasisWithOperatorsType BasisWithOperatorsType ;
		typedef typename BasisWithOperatorsType::RealType RealType;

		typedef typename VectorType::value_type FieldType;
		typedef PsimagLite::Profiling ProfilingType;
		typedef typename BasisWithOperatorsType::OperatorType OperatorType;
		typedef ApplyOperatorLocal<BasisWithOperatorsType,VectorWithOffsetType,VectorType> ApplyOperatorType;

		enum {GROW_RIGHT,GROW_LEFT};
		enum {DIAGONAL,NON_DIAGONAL};
		enum {GS_VECTOR=ObserverHelperType::GS_VECTOR,TIME_VECTOR=ObserverHelperType::TIME_VECTOR};
		enum {LEFT_BRACKET=ObserverHelperType::LEFT_BRACKET,RIGHT_BRACKET=ObserverHelperType::RIGHT_BRACKET};

		CorrelationsSkeleton(ObserverHelperType& helper,bool verbose = false)
		: helper_(helper),verbose_(verbose) { }
		
		size_t numberOfSites() const { return helper_.basisSE().block().size(); }

		//! i can be zero here!!
		void growDirectly(MatrixType& Odest,const MatrixType& Osrc,size_t i,
				int fermionicSign,size_t ns,bool transform = true)
		{
			//ProfilingType profile("growDirectly "+utils::ttos(i)+" ns="+utils::ttos(ns));
			Odest =Osrc;
			//std::vector<int> signs;
			// from 0 --> i
			int nt=i-1;
			if (nt<0) nt=0;
			
			for (size_t s=nt;s<ns;s++) {
				// set appropriate privates which are:
				// SpermutationInverse_(s) and transform_(s)
				/*io_.rewind();
				calcSpermutation(s);
				//std::cerr<<"*****************======="<<transform_.n_row()<<"\n";
				io_.readMatrix(transform_,"#TRANSFORM_sites",s);*/
				helper_.setPointer(s);
				//std::cerr<<"%%%%%%%%%%%%%%%%%======="<<transform_.n_row()<<"\n";
				int growOption = GROW_RIGHT;
				//if (i==1 && s==0) growOption = GROW_LEFT;// <-- not needed since nt>0
				if (s==size_t(nt)) {
					growOption = GROW_LEFT;
					if (i==0) growOption = GROW_RIGHT;
				}
				/*io_.rewind();
				io_.read(electrons_,"#ELECTRONS_sites=",s);*/
				//createSigns(signs,fermionicSign);
				MatrixType Onew(helper_.columns(),helper_.columns());
				fluffUp(Onew,Odest,fermionicSign,growOption,false);
				if (!transform && s == ns-1) {
					Odest = Onew;
					continue;
				}
				helper_.transform(Odest,Onew);
			}
		}
		
		// Perfomance critical:
		void fluffUp(MatrixType& ret2,const MatrixType& O,int fermionicSign,
			     int growOption=GROW_RIGHT,bool transform = true)
		{
			size_t n = helper_.basisS().size();
			MatrixType ret(n,n);

			for (size_t e=0;e<n;e++) {
				for (size_t e2=0;e2<n;e2++) {
					ret(e,e2) = fluffUp(O,e,e2,fermionicSign,growOption);
				}
			}
			if (transform) helper_.transform(ret2,ret);
			else ret2 = ret;
		}

		// Perfomance critical:
		FieldType fluffUp(const MatrixType& O,size_t e,size_t e2,int fermionicSign,
					  int growOption=GROW_RIGHT)	
		{
			
			size_t n = O.n_row();
			size_t m = size_t(helper_.basisS().size()/n);
			RealType sign = static_cast<RealType>(1.0);
			
			// Sperm[e] = i +k*n or e= k + i*m
			// Sperm[e2] = j+k*n or e2=k+j*m
			size_t i,j,k,k2;
			if (growOption==GROW_RIGHT) {	
				if (size_t(helper_.basisS().permutation(e)/n)!=size_t(helper_.basisS().permutation(e2)/n)) return 0;
				utils::getCoordinates(i,k,helper_.basisS().permutation(e),n);
				utils::getCoordinates(j,k2,helper_.basisS().permutation(e2),n);
			} else {
				if (size_t(helper_.basisS().permutation(e)%m)!=size_t(helper_.basisS().permutation(e2)%m)) return 0;
				utils::getCoordinates(k,i,helper_.basisS().permutation(e),m);
				utils::getCoordinates(k2,j,helper_.basisS().permutation(e2),m);
				sign = helper_.fermionicSign()(k,fermionicSign); // signs[k];
			}
			FieldType ret=0;
			if (k==k2) ret=O(i,j)*sign;
			return ret;
		}
		
		void dmrgMultiply(MatrixType& result,
					const MatrixType& O1,
					const MatrixType& O2,
					int fermionicSign,
				 	size_t ns)
		{
//			size_t numberOfSites = helper_.basisSE().block().size();
//			if (ns == numberOfSites-2) {
//				dmrgMultiplyRightCorner(result,O1,O2,fermionicSign,ns);
//				return;
//			}
			size_t ni=O1.n_row();
			size_t nj=O2.n_row();

			helper_.setPointer(ns);
			size_t sprime = helper_.basisS().size(); //ni*nj;
			result.resize(sprime,sprime);

			for (size_t r=0;r<result.n_row();r++)
				for (size_t r2=0;r2<result.n_col();r2++)
					result(r,r2)=0;

			if (helper_.basisS().size()!=sprime) {
				std::cerr<<"WARNING: "<<helper_.basisS().size()<<"!="<<sprime<<"\n";
				throw std::runtime_error("problem in dmrgMultiply\n");
			}

			for (size_t r=0;r<sprime;r++) {
				size_t e,u;
				utils::getCoordinates(e,u,helper_.basisS().permutation(r),ni);
				RealType f = helper_.fermionicSign()(e,fermionicSign);
				for (size_t e2=0;e2<ni;e2++) {
					for (size_t u2=0;u2<nj;u2++) {
						size_t r2 = helper_.basisS().permutationInverse(e2 + u2*ni);
						result(r,r2) += O1(e,e2)*O2(u,u2)*f;
					}
				}
			}
//			if (result.n_row()!=helper_.rows()) {
//				std::cerr<<result.n_row()<<" "<<helper_.rows()<<"\n";
//				throw std::runtime_error("dmrgMultiply: mismatch in transform\n");
//			}
		}

		void createWithModification(MatrixType& Om,const MatrixType& O,char mod)
		{
			//ProfilingType profile("create with modification="+mod);
			if (mod == 'n' || mod == 'N') {
				Om = O;
				return;
			}
			utils::transposeConjugate(Om,O);
		}

		FieldType bracket(const MatrixType& A)
		{
			try {
				const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET);
				const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET);

				return bracket_(A,src1,src2);
			} catch (std::exception& e) {
				std::cerr<<"CAUGHT: "<<e.what();
				std::cerr<<"WARNING: CorrelationsSkeleton::bracket(...):";
				std::cerr<<" No data seen yet\n";
				return 0;
			}
		}

		FieldType bracketRightCorner(const MatrixType& A,const MatrixType& B,int fermionSign)
		{
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET);
			const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET);

			return bracketRightCorner_(A,B,fermionSign,src1,src2);
		}

		FieldType bracketRightCorner(const MatrixType& A,const MatrixType& B,
				const MatrixType& C,int fermionSign)
		{
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET);
			const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET);

			return bracketRightCorner_(A,B,C,fermionSign,src1,src2);
		}
	private:
		
		//template<typename SomeVectorType>
		RealType bracket_(const MatrixType& A,const VectorWithOffsetType& vec1,const VectorWithOffsetType& vec2)
		{
			RealType norma = std::norm(vec1);
			if (verbose_) std::cerr<<"SE.size="<<helper_.basisSE().size()<<"\n";
			
			CrsMatrix<FieldType> Acrs(A);
			FieldType sum=0;

			if (vec1.size()!=helper_.basisSE().size() || vec1.size()!=vec2.size()) throw std::runtime_error("Error\n");
			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					utils::getCoordinates(r,eta,helper_.basisSE().permutation(t),helper_.basisS().size());
					for (int k=Acrs.getRowPtr(r);k<Acrs.getRowPtr(r+1);k++) {
						size_t r2 = Acrs.getCol(k);
						size_t t2 = helper_.basisSE().permutationInverse(r2+eta*A.n_col());
						if (t2<offset || t2>=total) continue;
						sum += Acrs.getValue(k)*vec1[t]*std::conj(vec2[t2]);
					}
				}
			}
			return std::real(sum)/norma;
		}
		
		RealType bracketRightCorner_(
				const MatrixType& A,
				const MatrixType& B,
				int fermionSign,
				const VectorWithOffsetType& vec1,
				const VectorWithOffsetType& vec2)
		{

			RealType norma = std::norm(vec1);

			if (verbose_) std::cerr<<"SE.size="<<helper_.basisSE().size()<<"\n";

			CrsMatrix<FieldType> Acrs(A);
			CrsMatrix<FieldType> Bcrs(B);
			FieldType sum=0;
			size_t ni = helper_.basisS().size()/Bcrs.rank(); // = Acrs.rank()

			// some sanity checks:
			if (vec1.size()!=vec2.size() || vec1.size()!=helper_.basisSE().size())
				throw std::runtime_error("Observe::bracketRightCorner_(...): vec.size!=SE.size\n");
			if (ni!=Acrs.rank())
				throw std::runtime_error("Observe::bracketRightCorner_(...): ni!=Acrs.rank\n");

			// ok, we're ready for the main course:
			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					utils::getCoordinates(r,eta,helper_.basisSE().permutation(t),helper_.basisS().size());
					size_t r0,r1;
					utils::getCoordinates(r0,r1,helper_.basisS().permutation(r),ni);
					size_t electrons = helper_.basisSE().electrons(t);
					electrons -= helper_.basisE().electrons(eta); //helper_.electrons(eta);
					RealType sign = (electrons & 1) ? fermionSign : 1.0;
					//RealType sign = helper_.basisS().fermionicSign(r0,fermionSign);
					//sign *= helper_.basisE().fermionicSign(r1,fermionSign);

					for (int k=Acrs.getRowPtr(r0);k<Acrs.getRowPtr(r0+1);k++) {
						size_t r0prime = Acrs.getCol(k);
						for (int k2 = Bcrs.getRowPtr(eta);k2<Bcrs.getRowPtr(eta+1);k2++) {
							size_t eta2 = Bcrs.getCol(k2);
							size_t rprime = helper_.basisS().permutationInverse(r0prime+r1*ni);
							size_t t2 = helper_.basisSE().permutationInverse(rprime+eta2*helper_.basisS().size());
							if (t2<offset || t2>=total) continue;
							sum += Acrs.getValue(k)*Bcrs.getValue(k2)*vec1[t]*std::conj(vec2[t2])*sign;
						}
					}
				}
			}
			return std::real(sum)/norma;
		}

		RealType bracketRightCorner_(
				const MatrixType& A1,
				const MatrixType& A2,
				const MatrixType& B,
				int fermionSign,
				const VectorWithOffsetType& vec1,
				const VectorWithOffsetType& vec2)
		{

			RealType norma = std::norm(vec1);

			if (verbose_) std::cerr<<"SE.size="<<helper_.basisSE().size()<<"\n";

			CrsMatrix<FieldType> A1crs(A1);
			CrsMatrix<FieldType> A2crs(A2);
			CrsMatrix<FieldType> Bcrs(B);
			FieldType sum=0;
			size_t ni = helper_.basisS().size()/Bcrs.rank(); // = Acrs.rank()

			// some sanity checks:
			if (vec1.size()!=vec2.size() || vec1.size()!=helper_.basisSE().size())
				throw std::runtime_error("Observe::bracketRightCorner_(...): vec.size!=SE.size\n");
			if (ni!=A1crs.rank())
				throw std::runtime_error("Observe::bracketRightCorner_(...): ni!=A1crs.rank\n");
			if (Bcrs.rank()!=A2crs.rank())
				throw std::runtime_error("Observe::bracketRightCorner_(...): Bcrs.rank!=A2crs.rank\n");

			// ok, we're ready for the main course:
			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					utils::getCoordinates(r,eta,helper_.basisSE().permutation(t),helper_.basisS().size());
					size_t r0,r1;
					utils::getCoordinates(r0,r1,helper_.basisS().permutation(r),ni);
					RealType sign =  helper_.basisE().fermionicSign(r1,fermionSign);

					for (int k1=A1crs.getRowPtr(r0);k1<A1crs.getRowPtr(r0+1);k1++) {
						size_t r0prime = A1crs.getCol(k1);
						for (int k2=A2crs.getRowPtr(r1);k2<A2crs.getRowPtr(r1+1);k2++) {
							size_t r1prime = A2crs.getCol(k2);
								for (int k3 = Bcrs.getRowPtr(eta);k3<Bcrs.getRowPtr(eta+1);k3++) {
									size_t eta2 = Bcrs.getCol(k3);
									size_t rprime = helper_.basisS().permutationInverse(r0prime+r1prime*ni);
									size_t t2 = helper_.basisSE().permutationInverse(rprime+eta2*helper_.basisS().size());
									if (t2<offset || t2>=total) continue;
									sum += A1crs.getValue(k1)*A2crs.getValue(k2)*Bcrs.getValue(k3)*vec1[t]*std::conj(vec2[t2])*sign;
								}
						}
					}
				}
			}
			return std::real(sum)/norma;
		}

		ObserverHelperType& helper_; //<-- NB: We are not the owner
		bool verbose_;
	};  //class CorrelationsSkeleton
} // namespace Dmrg

/*@}*/
#endif
