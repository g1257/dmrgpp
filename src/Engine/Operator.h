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
/** \ingroup DMRG */
/*@{*/

/*! \file Operator.h
 *
 *  A structure to represent an operator
 *  Contains the actual data, the (J,M) that indicates
 * how this operator transforms, the fermionSign which
 * indicates if this operator commutes or anticommutes
 * with operators of the same class on different sites, and
 * other properties.
 *
 */
#ifndef OPERATOR_H
#define OPERATOR_H

namespace Dmrg {
	//! This is a structure, don't add member functions here!
	struct Su2Related {
		Su2Related()
		: offset(0) // setting to zero is necessary, because
		{}		// we always print offset
				// and when running Abelian
				// it might be undefined
		size_t offset;
		PsimagLite::Vector<size_t>::Type source;
		PsimagLite::Vector<int>::Type transpose;
	};

	inline std::istream& operator>>(std::istream& is,Su2Related& x)
	{
		is>>x.offset;
		return is;
	}

	inline std::ostream& operator<<(std::ostream& os,const Su2Related& x)
	{
		os<<x.offset<<"\n";
		return os;
	}

	//! This is a structure, don't add member functions here!
	template<typename RealType_,typename SparseMatrixType_>
	struct Operator {
		typedef RealType_ RealType;
		typedef SparseMatrixType_ SparseMatrixType;
		typedef std::pair<size_t,size_t> PairType;
		typedef Su2Related Su2RelatedType;
		Operator() {}

		Operator(const SparseMatrixType& data1,
		         int fermionSign1,
		         const PairType& jm1,
		         RealType angularFactor1,
		         const Su2RelatedType& su2Related1)
		: data(data1),fermionSign(fermionSign1),jm(jm1),angularFactor(angularFactor1),su2Related(su2Related1)
		{}

		SparseMatrixType data;
		int fermionSign; // does this operator commute or anticommute with others of the same class on different sites
		PairType  jm; // angular momentum of this operator	
		RealType angularFactor;

		Su2RelatedType su2Related;
	};

	template<typename RealType_,
	         typename SparseMatrixType,
	         template<typename,typename> class SomeVectorTemplate,
	         template<typename,int> class SomeAllocatorTemplate,
	         int n>
	void fillOperator(SomeVectorTemplate<SparseMatrixType*,SomeAllocatorTemplate<SparseMatrixType*,n> >& data,
	                  SomeVectorTemplate<Operator<RealType_,SparseMatrixType>,SomeAllocatorTemplate<Operator<RealType_,SparseMatrixType>,n> >& op)
	{
		for (size_t i=0;i<data.size();i++) {
			data[i] = &(op[i].data);
		}
	}

	template<typename RealType_,
	         typename SparseMatrixType,
	         typename ConcurrencyType,
	         template<typename,typename> class SomeVectorTemplate,
	         template<typename,int> class SomeAllocatorTemplate,
	         int n>
	void gather(SomeVectorTemplate<Operator<RealType_,SparseMatrixType>,SomeAllocatorTemplate<Operator<RealType_,SparseMatrixType>,n> >& op,
	            ConcurrencyType& concurrency)
	{
		SomeVectorTemplate<SparseMatrixType*,SomeAllocatorTemplate<SparseMatrixType*,n> > data(op.size());
//		PsimagLite::Vector<int*>::Type fermionSign(op.size());
//		typename PsimagLite::Vector<typename Operator<RealType_,SparseMatrixType>::PairType*>::Type jm(op.size());
//		typename PsimagLite::Vector<RealType_*>::Type angularFactor(op.size());
//		PsimagLite::Vector<typename Operator<RealType_,SparseMatrixType>::Type::Su2RelatedType*> su2Related(op.size());

		fillOperator(data,op);
		concurrency.gather(data);
	}

	template<typename RealType_,
	         typename SparseMatrixType,
	         typename ConcurrencyType,
	         template<typename,typename> class SomeVectorTemplate,
	         template<typename,int> class SomeAllocatorTemplate,
	         int n>
	void broadcast(SomeVectorTemplate<Operator<RealType_,SparseMatrixType>,SomeAllocatorTemplate<Operator<RealType_,SparseMatrixType>,n> >& op,
	               ConcurrencyType& concurrency)
	{
		SomeVectorTemplate<SparseMatrixType*,SomeAllocatorTemplate<SparseMatrixType*,n> > data(op.size());
//		PsimagLite::Vector<int*>::Type fermionSign(op.size());
//		typename PsimagLite::Vector<typename Operator<RealType_,SparseMatrixType>::PairType*>::Type jm(op.size());
//		typename PsimagLite::Vector<RealType_*>::Type angularFactor(op.size());
//		typename PsimagLite::Vector<typename Operator<RealType_,SparseMatrixType>::Su2RelatedType*>::Type su2Related(op.size());

		fillOperator(data,op);

		concurrency.broadcast(data);
	}

	template<typename RealType,typename SparseMatrixType>
	std::istream& operator>>(std::istream& is,Operator<RealType,SparseMatrixType>& op)
	{
		is>>op.data;
		is>>op.fermionSign;
		is>>op.jm;
		is>>op.angularFactor;
		is>>op.su2Related;
		return is;
	}

	template<typename RealType,typename SparseMatrixType>
	std::ostream& operator<<(std::ostream& os,const Operator<RealType,SparseMatrixType>& op)
	{
		os<<op.data;
		os<<op.fermionSign<<"\n";
		os<<op.jm.first<<" "<<op.jm.second<<"\n";
		os<<op.angularFactor<<"\n";
		os<<op.su2Related;
		return os;
	}
} // namespace Dmrg

/*@}*/
#endif

