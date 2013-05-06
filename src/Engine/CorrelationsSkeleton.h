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
#include "Matrix.h"
#include "PackIndices.h"
#include "CrsMatrix.h"
#include "Profiling.h"
#include "ApplyOperatorLocal.h"

namespace Dmrg {
	
	//! Don't add functions to this class
	template<typename FieldType>
	struct CorrelationData {
		size_t ni;
		size_t nj;
		typename PsimagLite::Vector<PsimagLite::Matrix<FieldType> >::Type correlationVector;
		PsimagLite::SparseVector<FieldType> wavefunction;
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
	
	template<typename ObserverHelperType_,typename ModelType>
	class CorrelationsSkeleton {
		typedef size_t IndexType;
		typedef PsimagLite::PackIndices PackIndicesType;

	public:
		typedef ObserverHelperType_ ObserverHelperType;
		typedef typename ObserverHelperType::MatrixType MatrixType;
		typedef typename ObserverHelperType::VectorType VectorType ;
		typedef typename ObserverHelperType::VectorWithOffsetType VectorWithOffsetType;
		typedef typename ObserverHelperType::BasisWithOperatorsType BasisWithOperatorsType ;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename ObserverHelperType::FermionSignType FermionSignType;

		typedef typename VectorType::value_type FieldType;
		typedef PsimagLite::Profiling ProfilingType;
		typedef typename BasisWithOperatorsType::OperatorType OperatorType;
		typedef ApplyOperatorLocal<BasisWithOperatorsType,VectorWithOffsetType,VectorType> ApplyOperatorType;
		typedef PsimagLite::CrsMatrix<FieldType> SparseMatrixType;

		enum {GROW_RIGHT,GROW_LEFT};
		enum {DIAGONAL,NON_DIAGONAL};
		enum {GS_VECTOR=ObserverHelperType::GS_VECTOR,
			TIME_VECTOR=ObserverHelperType::TIME_VECTOR};
		enum {LEFT_BRACKET=ObserverHelperType::LEFT_BRACKET,
			RIGHT_BRACKET=ObserverHelperType::RIGHT_BRACKET};
		static const size_t EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
		static const size_t EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;

		CorrelationsSkeleton(
				ObserverHelperType& helper,
				const ModelType& model,
				bool verbose = false)
		: helper_(helper),verbose_(verbose)
		{}

		size_t numberOfSites() const
		{
			size_t threadId = 0;
			return helper_.leftRightSuper(threadId).sites();
		}

		//! i can be zero here!!
		void growDirectly(MatrixType& Odest,const MatrixType& Osrc,size_t i,
				int fermionicSign,size_t ns,bool transform,size_t threadId)
		{
			Odest =Osrc;
			// from 0 --> i
			int nt=i-1;
			if (nt<0) nt=0;
			
			for (size_t s=nt;s<ns;s++) {
				helper_.setPointer(threadId,s);
				size_t growOption = growthDirection(s,nt,i,threadId);
				MatrixType Onew(helper_.columns(threadId),helper_.columns(threadId));
				fluffUp(Onew,Odest,fermionicSign,growOption,false,threadId);
				if (!transform && s == ns-1) {
					Odest = Onew;
					continue;
				}
				helper_.transform(Odest,Onew,threadId);
			}
		}
		
		size_t growthDirection(size_t s,int nt,size_t i,size_t threadId) const
		{
			size_t dir = helper_.direction(threadId);
			size_t growOption = (dir==EXPAND_SYSTEM) ?
					GROW_RIGHT : GROW_LEFT;

			if (s==size_t(nt)) {
				growOption = (dir==EXPAND_SYSTEM) ?
						GROW_LEFT : GROW_RIGHT;
				if (i==0) growOption = (dir==EXPAND_SYSTEM) ?
						GROW_RIGHT : GROW_LEFT;
			}
			return growOption;
		}

		// Perfomance critical:
		void fluffUp(MatrixType& ret2,
					 const MatrixType& O,
					 int fermionicSign,
					 int growOption,//=GROW_RIGHT,
					 bool transform, //= true
					 size_t threadId)
		{
			if (helper_.direction(threadId)==EXPAND_SYSTEM) {
				fluffUpSystem(ret2,O,fermionicSign,growOption,transform,threadId);
				return;
			}
			fluffUpEnviron(ret2,O,fermionicSign,growOption,transform,threadId);
		}
		
		void dmrgMultiply(MatrixType& result,
						  const MatrixType& O1,
						  const MatrixType& O2,
						  int fermionicSign,
						  size_t ns,
						  size_t threadId)
		{
			if (helper_.direction(threadId)==EXPAND_SYSTEM) {
				dmrgMultiplySystem(result,O1,O2,fermionicSign,ns,threadId);
				return;
			}
			dmrgMultiplyEnviron(result,O1,O2,fermionicSign,ns,threadId);
		}

		void createWithModification(MatrixType& Om,const MatrixType& O,char mod)
		{
			if (mod == 'n' || mod == 'N') {
				Om = O;
				return;
			}
			transposeConjugate(Om,O);
		}

		FieldType bracket(const MatrixType& A,int fermionicSign,size_t threadId)
		{
			try {
				const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET,threadId);
				const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET,threadId);

				return bracket_(A,src1,src2,fermionicSign,threadId);
			} catch (std::exception& e) {
				std::cerr<<"CAUGHT: "<<e.what();
				std::cerr<<"WARNING: CorrelationsSkeleton::bracket(...):";
				std::cerr<<" No data seen yet\n";
				return 0;
			}
		}

		FieldType bracketRightCorner(const MatrixType& A,const MatrixType& B,int fermionSign,size_t threadId)
		{
			try {
				const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET,threadId);
				const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET,threadId);
				return bracketRightCorner_(A,B,fermionSign,src1,src2,threadId);
			}  catch (std::exception& e) {
				std::cerr<<"CAUGHT: "<<e.what();
				std::cerr<<"WARNING: CorrelationsSkeleton::bracketRightCorner(...):";
				std::cerr<<" No data seen yet\n";
				return 0;
			}
		}

		FieldType bracketRightCorner(const MatrixType& A,const MatrixType& B,
				const MatrixType& C,int fermionSign,size_t threadId)
		{
			try {
				const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET,threadId);
				const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET,threadId);
				return bracketRightCorner_(A,B,C,fermionSign,src1,src2,threadId);
			} catch (std::exception& e) {
				std::cerr<<"CAUGHT: "<<e.what();
				std::cerr<<"WARNING: CorrelationsSkeleton::bracketRightCornerABC(...):";
				std::cerr<<" No data seen yet\n";                                                               
				return 0;                                                                                                                       
			}
		}

	private:
		
		void dmrgMultiplySystem(MatrixType& result,
								const MatrixType& O1,
								const MatrixType& O2,
								int fermionicSign,
								size_t ns,
								size_t threadId)
		{
			size_t ni=O1.n_row();
			size_t nj=O2.n_row();

			helper_.setPointer(threadId,ns);
			size_t sprime = helper_.leftRightSuper(threadId).left().size(); //ni*nj;
			result.resize(sprime,sprime);

			for (size_t r=0;r<result.n_row();r++)
				for (size_t r2=0;r2<result.n_col();r2++)
					result(r,r2)=0;

			if (helper_.leftRightSuper(threadId).left().size()!=sprime) {
				std::cerr<<"WARNING: "<<helper_.leftRightSuper(threadId).left().size();
				std::cerr<<"!="<<sprime<<"\n";
				throw std::runtime_error("problem in dmrgMultiply\n");
			}

			PackIndicesType pack(ni);
			for (size_t r=0;r<sprime;r++) {
				size_t e,u;
				pack.unpack(e,u,helper_.leftRightSuper(threadId).left().permutation(r));
				RealType f = helper_.fermionicSignLeft(threadId)(e,fermionicSign);
				for (size_t e2=0;e2<ni;e2++) {
					for (size_t u2=0;u2<nj;u2++) {
						size_t r2 = helper_.leftRightSuper(threadId).left().
								permutationInverse(e2 + u2*ni);
						result(r,r2) += O1(e,e2)*O2(u,u2)*f;
					}
				}
			}
		}

		void dmrgMultiplyEnviron(MatrixType& result,
								 const MatrixType& O1,
								 const MatrixType& O2,
								 int fermionicSign,
								 size_t ns,
								 size_t threadId)
		{
			size_t ni=O1.n_row();
			size_t nj=O2.n_row();

			helper_.setPointer(threadId,ns);
			size_t eprime = helper_.leftRightSuper(threadId).right().size(); //ni*nj;
			result.resize(eprime,eprime);

			for (size_t r=0;r<result.n_row();r++)
				for (size_t r2=0;r2<result.n_col();r2++)
					result(r,r2)=0;

			if (helper_.leftRightSuper(threadId).right().size()!=eprime) {
				std::cerr<<"WARNING: "<<helper_.leftRightSuper(threadId).right().size();
				std::cerr<<"!="<<eprime<<"\n";
				throw std::runtime_error("problem in dmrgMultiply\n");
			}

			PackIndicesType pack(nj);
			for (size_t r=0;r<eprime;r++) {
				size_t e,u;
				pack.unpack(e,u,helper_.leftRightSuper(threadId).right().permutation(r));
				size_t nx0 = helper_.leftRightSuper(threadId).right().
							electrons(BasisType::AFTER_TRANSFORM);
				RealType f = (nx0 & 1) ? fermionicSign : 1;
				for (size_t e2=0;e2<nj;e2++) {
					for (size_t u2=0;u2<ni;u2++) {
						size_t r2 = helper_.leftRightSuper(threadId).right().
								permutationInverse(e2 + u2*nj);
						if (r2>=eprime) throw std::runtime_error("Error\n");
						result(r,r2) += O2(e,e2)*O1(u,u2)*f;
					}
				}
			}
		}

		// Perfomance critical:
		void fluffUpSystem(
				MatrixType& ret2,
				const MatrixType& O,
				int fermionicSign,
				int growOption,
				bool transform,
			size_t threadId)
		{
			size_t n =helper_.leftRightSuper(threadId).left().size();

			MatrixType ret(n,n);

			for (size_t e=0;e<n;e++) {
				for (size_t e2=0;e2<n;e2++) {
					ret(e,e2) = fluffUpSystem_(
							O,e,e2,fermionicSign,growOption,threadId);
				}
			}
			if (transform) helper_.transform(ret2,ret,threadId);
			else ret2 = ret;
		}

		// Perfomance critical:
		void fluffUpEnviron(
				MatrixType& ret2,
				const MatrixType& O,
				int fermionicSign,
				int growOption,
				bool transform,
			size_t threadId)
		{
			size_t n =helper_.leftRightSuper(threadId).right().size();

			MatrixType ret(n,n);
			for (size_t e=0;e<n;e++) {
				for (size_t e2=0;e2<n;e2++) {
					ret(e,e2) = fluffUpEnviron_(
							O,e,e2,fermionicSign,growOption,threadId);
				}
			}
			if (transform) helper_.transform(ret2,ret,threadId);
			else ret2 = ret;
		}

		// Perfomance critical:
		FieldType fluffUpSystem_(
				const MatrixType& O,
				size_t e,size_t e2,
				int fermionicSign,
				int growOption,
			size_t threadId)
		{
			size_t n = O.n_row();
			size_t m = size_t(helper_.leftRightSuper(threadId).left().size()/n);
			RealType sign = static_cast<RealType>(1.0);

			// Sperm[e] = i +k*n or e= k + i*m
			// Sperm[e2] = j+k*n or e2=k+j*m
			size_t i,j,k,k2;
			if (growOption==GROW_RIGHT) {
				if (size_t(helper_.leftRightSuper(threadId).left().permutation(e)/n)!=
						size_t(helper_.leftRightSuper(threadId).left().permutation(e2)/n)) return 0;
				PackIndicesType pack(n);
				pack.unpack(i,k,helper_.leftRightSuper(threadId).left().permutation(e));
				pack.unpack(j,k2,helper_.leftRightSuper(threadId).left().permutation(e2));
			} else {
				if (size_t(helper_.leftRightSuper(threadId).left().permutation(e)%m)!=
						size_t(helper_.leftRightSuper(threadId).left().permutation(e2)%m)) return 0;
				PackIndicesType pack(m);
				pack.unpack(k,i,helper_.leftRightSuper(threadId).left().permutation(e));
				pack.unpack(k2,j,helper_.leftRightSuper(threadId).left().permutation(e2));
				sign = helper_.fermionicSignLeft(threadId)(k,fermionicSign);
			}
			if (k!=k2) return 0;
			return O(i,j)*sign;
		}

		// Perfomance critical:
		FieldType fluffUpEnviron_(
				const MatrixType& O,
				size_t e,size_t e2,
				int fermionicSign,
				int growOption,
			size_t threadId)
		{
			size_t n = O.n_row();
			size_t m = size_t(helper_.leftRightSuper(threadId).right().size()/n);
			RealType sign = 1;

			// Eperm[e] = i +k*n or e= k + i*m
			// Eperm[e2] = j+k*n or e2=k+j*m

			size_t i,j,k,k2;

			if (growOption==GROW_RIGHT) {
				PackIndicesType pack(n);
				pack.unpack(i,k,helper_.leftRightSuper(threadId).right().permutation(e));
				pack.unpack(j,k2,helper_.leftRightSuper(threadId).right().permutation(e2));
				size_t nx0 = helper_.leftRightSuper(threadId).left().electrons(BasisType::AFTER_TRANSFORM);
				sign = (nx0 & 1) ? fermionicSign : 1;
			} else {
				PackIndicesType pack(m);
				pack.unpack(k,i,helper_.leftRightSuper(threadId).right().permutation(e));
				pack.unpack(k2,j,helper_.leftRightSuper(threadId).right().permutation(e2));
				size_t nx0 = helper_.leftRightSuper(threadId).super().electrons(BasisType::AFTER_TRANSFORM);
				sign = (nx0 & 1) ?  fermionicSign : 1;
			}
			if (k!=k2) return 0;
			return O(i,j)*sign;
		}

		RealType bracket_(
			const MatrixType& A,
			const VectorWithOffsetType& vec1,
			const VectorWithOffsetType& vec2,
			int fermionicSign,
			size_t threadId)
		{
			if (verbose_)
				std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

			if (vec1.size()!=helper_.leftRightSuper(threadId).super().size() ||
								vec1.size()!=vec2.size())
				throw std::runtime_error(
					"CorrelationsSkeleton::bracket_(...): Error\n");

			if (helper_.direction(threadId)==EXPAND_SYSTEM) {
				return bracketSystem_(A,vec1,vec2,threadId);
			}
			return bracketEnviron_(A,vec1,vec2,fermionicSign,threadId);
		}

		RealType bracketSystem_(
						const MatrixType& A,
						const VectorWithOffsetType& vec1,
						const VectorWithOffsetType& vec2,
			size_t threadId)
		{
			SparseMatrixType Acrs(A);
			FieldType sum=0;
			PackIndicesType pack(helper_.leftRightSuper(threadId).left().size());
			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					pack.unpack(r,eta,helper_.leftRightSuper(threadId).super().
							permutation(t));
					for (int k=Acrs.getRowPtr(r);k<Acrs.getRowPtr(r+1);k++) {
						size_t r2 = Acrs.getCol(k);
						size_t t2 = helper_.leftRightSuper(threadId).super().
								permutationInverse(r2+eta*A.n_col());
						if (t2<offset || t2>=total) continue;
						sum += Acrs.getValue(k)*vec1[t]*std::conj(vec2[t2]);
					}
				}
			}
			RealType norma = std::norm(vec1);
			return std::real(sum)/norma;
		}

		RealType bracketEnviron_(
						const MatrixType& A,
						const VectorWithOffsetType& vec1,
						const VectorWithOffsetType& vec2,
						int fermionicSign,
			size_t threadId)
		{
			SparseMatrixType Acrs(A);
			FieldType sum=0;
			PackIndicesType pack(helper_.leftRightSuper(threadId).left().size());

			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					pack.unpack(r,eta,helper_.leftRightSuper(threadId).super().
							permutation(t));
					if (eta>=Acrs.row()) throw std::runtime_error("Error\n");
					size_t nx0 = helper_.leftRightSuper(threadId).left().electrons(BasisType::AFTER_TRANSFORM);
					RealType sign = (nx0 & 1) ? fermionicSign : 1;
					for (int k=Acrs.getRowPtr(eta);k<Acrs.getRowPtr(eta+1);k++) {
						size_t eta2 = Acrs.getCol(k);
						size_t t2 = helper_.leftRightSuper(threadId).super().
							permutationInverse(r+eta2*helper_.leftRightSuper(threadId).left().size());
						if (t2<offset || t2>=total) continue;
						sum += Acrs.getValue(k)*vec1[t]*std::conj(vec2[t2])*sign;
					}
				}
			}
			RealType norma = std::norm(vec1);
			return std::real(sum)/norma;
		}
		
		RealType bracketRightCorner_(
			const MatrixType& A,
			const MatrixType& B,
			int fermionSign,
			const VectorWithOffsetType& vec1,
			const VectorWithOffsetType& vec2,
			size_t threadId)
		{
			if (helper_.direction(threadId)==EXPAND_SYSTEM)
				return brRghtCrnrSystem_(A,B,fermionSign,vec1,vec2,threadId);
			return brLftCrnrEnviron_(A,B,fermionSign,vec1,vec2,threadId);
		}

		RealType brRghtCrnrSystem_(
						const MatrixType& A,
						const MatrixType& B,
						int fermionSign,
						const VectorWithOffsetType& vec1,
						const VectorWithOffsetType& vec2,
			size_t threadId)
		{
			if (verbose_) std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

			SparseMatrixType Acrs(A);
			SparseMatrixType Bcrs(B);
			FieldType sum=0;
			size_t ni = helper_.leftRightSuper(threadId).left().size()/Bcrs.row(); // = Acrs.rank()

			// some sanity checks:
			if (vec1.size()!=vec2.size() || vec1.size()!=helper_.leftRightSuper(threadId).super().size())
				throw std::runtime_error("Observe::brRghtCrnrSystem_(...): "
						"vec.size!=SE.size\n");
			if (ni!=Acrs.row())
				throw std::runtime_error("Observe::brRghtCrnrSystem_(...): "
						"ni!=Acrs.rank\n");

			// ok, we're ready for the main course:
			PackIndicesType pack1(helper_.leftRightSuper(threadId).left().size());
			PackIndicesType pack2(ni);
			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					pack1.unpack(r,eta,helper_.leftRightSuper(threadId).super().
							permutation(t));
					size_t r0,r1;
					pack2.unpack(r0,r1,helper_.leftRightSuper(threadId).left().
							permutation(r));
					size_t electrons = helper_.leftRightSuper(threadId).super().electrons(t);
					electrons -= helper_.leftRightSuper(threadId).right().electrons(eta);
					RealType sign = (electrons & 1) ? fermionSign : 1.0;

					for (int k=Acrs.getRowPtr(r0);k<Acrs.getRowPtr(r0+1);k++) {
						size_t r0prime = Acrs.getCol(k);
						for (int k2 = Bcrs.getRowPtr(eta);
								k2<Bcrs.getRowPtr(eta+1);k2++) {
							size_t eta2 = Bcrs.getCol(k2);
							size_t rprime = helper_.leftRightSuper(threadId).left().
									permutationInverse(r0prime+r1*ni);
							size_t t2 = helper_.leftRightSuper(threadId).super().permutationInverse(
									rprime+eta2*helper_.leftRightSuper(threadId).left().size());
							if (t2<offset || t2>=total) continue;
							sum += Acrs.getValue(k)*Bcrs.getValue(k2)*
									vec1[t]*std::conj(vec2[t2])*sign;
						}
					}
				}
			}
			RealType norma = std::norm(vec1);
			return std::real(sum)/norma;
		}

		RealType brLftCrnrEnviron_(
				const MatrixType& A,
				const MatrixType& B,
				int fermionSign,
				const VectorWithOffsetType& vec1,
				const VectorWithOffsetType& vec2,
			size_t threadId)
		{
			if (verbose_) std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

			SparseMatrixType Acrs(A);
			SparseMatrixType Bcrs(B);
			FieldType sum=0;
			size_t ni = Bcrs.row();

			// some sanity checks:
			if (vec1.size()!=vec2.size() || vec1.size()!=helper_.leftRightSuper(threadId).super().size())
				throw std::runtime_error("Observe::brLftCrnrEnviron_(...): "
						"vec.size!=SE.size\n");
			if (helper_.leftRightSuper(threadId).right().size()/Bcrs.row()!=Acrs.row())
				throw std::runtime_error("Observe::bracketRightCorner_(...): "
						"helper_.leftRightSuper().right().size()/Bcrs.rank()!=Acrs.rank\n");

			// ok, we're ready for the main course:
			PackIndicesType pack1(helper_.leftRightSuper(threadId).left().size());
			PackIndicesType pack2(ni);

			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					pack1.unpack(eta,r,helper_.leftRightSuper(threadId).super().
							permutation(t));
					size_t r0,r1;
					pack2.unpack(r0,r1,helper_.leftRightSuper(threadId).right().permutation(r));
					size_t electrons = helper_.leftRightSuper(threadId).left().electrons(eta);
					RealType sign = (electrons & 1) ? fermionSign : 1.0;

					for (int k=Acrs.getRowPtr(r1);k<Acrs.getRowPtr(r1+1);k++) {
						size_t r1prime = Acrs.getCol(k);
						for (int k2 = Bcrs.getRowPtr(eta);
								k2<Bcrs.getRowPtr(eta+1);k2++) {
							size_t eta2 = Bcrs.getCol(k2);
							size_t rprime = helper_.leftRightSuper(threadId).right().
									permutationInverse(r0+r1prime*ni);
							size_t t2 = helper_.leftRightSuper(threadId).super().permutationInverse(
									eta2+rprime*helper_.leftRightSuper(threadId).left().size());
							if (t2<offset || t2>=total) continue;
							sum += Acrs.getValue(k)*Bcrs.getValue(k2)*vec1[t]*
									std::conj(vec2[t2])*sign;
						}
					}
				}
			}
			RealType norma = std::norm(vec1);
			return std::real(sum)/norma;
		}

		RealType bracketRightCorner_(
			const MatrixType& A1,
			const MatrixType& A2,
			const MatrixType& B,
			int fermionSign,
			const VectorWithOffsetType& vec1,
			const VectorWithOffsetType& vec2,
			size_t threadId)
		{
			if (helper_.direction(threadId)!=EXPAND_SYSTEM) return 0;

			RealType norma = std::norm(vec1);

			if (verbose_) std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

			SparseMatrixType A1crs(A1);
			SparseMatrixType A2crs(A2);
			SparseMatrixType Bcrs(B);
			FieldType sum=0;
			size_t ni = helper_.leftRightSuper(threadId).left().size()/Bcrs.row(); // = Acrs.rank()

			// some sanity checks:
			assert(vec1.size()==vec2.size());

			if (vec1.size()==0) return 0;

			assert(vec1.size()==helper_.leftRightSuper(threadId).super().size());
			assert(ni==A1crs.row());
			assert(Bcrs.row()==A2crs.row());

			// ok, we're ready for the main course:
			PackIndicesType pack1(helper_.leftRightSuper(threadId).left().size());
			PackIndicesType pack2(ni);

			for (size_t x=0;x<vec1.sectors();x++) {
				size_t sector = vec1.sector(x);
				size_t offset = vec1.offset(sector);
				size_t total = offset + vec1.effectiveSize(sector);
				for (size_t t=offset;t<total;t++) {
					size_t eta,r;

					pack1.unpack(r,eta,
							helper_.leftRightSuper(threadId).super().permutation(t));
					size_t r0,r1;
					pack2.unpack(r0,r1,
							helper_.leftRightSuper(threadId).left().permutation(r));
					RealType sign =  helper_.leftRightSuper(threadId).right().fermionicSign(r1,fermionSign);

					for (int k1=A1crs.getRowPtr(r0);k1<A1crs.getRowPtr(r0+1);k1++) {
						size_t r0prime = A1crs.getCol(k1);
						for (int k2=A2crs.getRowPtr(r1);k2<A2crs.getRowPtr(r1+1);k2++) {
							size_t r1prime = A2crs.getCol(k2);
								for (int k3 = Bcrs.getRowPtr(eta);k3<Bcrs.getRowPtr(eta+1);k3++) {
									size_t eta2 = Bcrs.getCol(k3);
									size_t rprime = helper_.leftRightSuper(threadId).left().permutationInverse(r0prime+r1prime*ni);
									size_t t2 = helper_.leftRightSuper(threadId).super().permutationInverse(rprime+eta2*helper_.leftRightSuper(threadId).left().size());
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
