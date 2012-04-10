// BEGIN LICENSE BLOCK
/*
Copyright (c) 2008 , UT-Battelle, LLC
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

/*! \file Observer.h
 *
 *  A class to perform post-processing calculation of observables
 *
 */
#ifndef DMRG_OBSERVE_H
#define DMRG_OBSERVE_H
#include "ObserverHelper.h"
#include "OnePointCorrelations.h"
#include "CorrelationsSkeleton.h"
#include "TwoPointCorrelations.h"
#include "FourPointCorrelations.h"
#include "VectorWithOffsets.h" // for operator*
#include "VectorWithOffset.h" // for operator*
#include "Profiling.h"

namespace Dmrg {
	
	template<
		typename FieldType,
		typename VectorWithOffsetType,
		typename ModelType,
		typename IoInputType>
	class Observer {
		typedef SparseVector<FieldType> VectorType;
		typedef typename ModelType::ConcurrencyType ConcurrencyType;
		typedef typename ModelType::RealType RealType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;
		typedef typename ModelType::BasisWithOperatorsType
				BasisWithOperatorsType;
		typedef typename ModelType::ModelHelperType::LeftRightSuperType
				LeftRightSuperType;
		typedef ObserverHelper<IoInputType,MatrixType,VectorType,
			VectorWithOffsetType,LeftRightSuperType> ObserverHelperType;
		typedef CorrelationsSkeleton<ObserverHelperType,ModelType>
			CorrelationsSkeletonType;
		typedef OnePointCorrelations<ObserverHelperType>
			OnePointCorrelationsType;
		typedef TwoPointCorrelations<CorrelationsSkeletonType,ConcurrencyType>
			TwoPointCorrelationsType;
		typedef FourPointCorrelations<CorrelationsSkeletonType>
			FourPointCorrelationsType;
		typedef PsimagLite::Profiling ProfilingType;

		static size_t const GROW_RIGHT = CorrelationsSkeletonType::GROW_RIGHT;
		static size_t const GROW_LEFT = CorrelationsSkeletonType::GROW_LEFT;
		static size_t const DIAGONAL = CorrelationsSkeletonType::DIAGONAL;
		static size_t const NON_DIAGONAL = CorrelationsSkeletonType::NON_DIAGONAL;
		enum {GS_VECTOR=ObserverHelperType::GS_VECTOR,
			TIME_VECTOR=ObserverHelperType::TIME_VECTOR};
		enum {LEFT_BRACKET=ObserverHelperType::LEFT_BRACKET,
			RIGHT_BRACKET=ObserverHelperType::RIGHT_BRACKET};

	public:
		Observer(
				IoInputType& io,
				size_t nf,
				bool hasTimeEvolution,
				const ModelType& model,
				ConcurrencyType& concurrency,
				bool verbose=false)
		: helper_(io,nf,hasTimeEvolution,verbose),
		  concurrency_(concurrency),
		  verbose_(verbose),
		  onepoint_(helper_),
		  skeleton_(helper_,model,verbose),
		  twopoint_(helper_,skeleton_,concurrency_),
		  fourpoint_(helper_,skeleton_)
		{}

		size_t size() const { return helper_.size(); }

		RealType time() const { return helper_.time(); }

		size_t site() const { return helper_.site(); }

		size_t marker() const { return helper_.marker(); }

		void setPointer(size_t x) { helper_.setPointer(x); }

		bool endOfData() const { return helper_.endOfData(); }

		// return true if
		// we're at site 1 or n-2
		bool isAtCorner(size_t numberOfSites) const
		{
			bool es = (helper_.direction() == ProgramGlobals::EXPAND_SYSTEM);
			if (es && helper_.site() ==  numberOfSites-2) return true;
			if ((!es) && helper_.site()==1) return true;
			return false;
		}

		void setBrackets(const std::string& left,const std::string& right)
		{
			helper_.setBrackets(bracketStringToNumber(left),
					bracketStringToNumber(right));
		}

		PsimagLite::Matrix<FieldType> correlations(
				const MatrixType& O1,
				const MatrixType& O2,
				int fermionicSign,
				size_t rows,
				size_t cols)
		{
			return twopoint_(O1,O2,fermionicSign,rows,cols);
		}

		FieldType fourPoint(
				char mod1,size_t i1,const MatrixType& O1,
				char mod2,size_t i2,const MatrixType& O2,
				char mod3,size_t i3,const MatrixType& O3,
				char mod4,size_t i4,const MatrixType& O4,
				int fermionicSign)
		{
			return fourpoint_(mod1,i1,O1,mod2,i2,O2,mod3,i3,O3,mod4,i4,O4,fermionicSign);
		}

		template<typename SomeModelType>
		void fourPointDeltas(std::vector<FieldType>& fpd,
				size_t n,
				const std::vector<size_t>& gammas,
				const SomeModelType& model)
		{
			if (gammas.size()!=4) {
				std::cerr<<"Observer: fourPointDeltas(...):  wrong number of gammas ";
				std::cerr<<" expected "<<4<<" got "<<gammas.size()<<"\n";
				throw std::runtime_error("Observer::fourPointDeltas(...)\n");
			}
			for (size_t i=0;i<n;i++) {
				if (i>0 && i%2!=0) continue;
				for (size_t j=i+2;j<n;j++) {
					//if (j== n-2) continue;
					if (j%2!=0) continue;
					if (i+1>=n-1) continue;
					fpd.push_back(fourPointDelta(i,j,gammas,model));
				}
			}
		}

		template<typename ApplyOperatorType>
		FieldType onePoint(size_t site,
				   const typename ApplyOperatorType::OperatorType& A,
				   bool corner = false)
		{
			return onepoint_.template operator()<ApplyOperatorType>(site,A,corner);
		}

	private:

		size_t bracketStringToNumber(const std::string& str) const
		{
			if (str=="gs") return GS_VECTOR;
			if (str=="time") return TIME_VECTOR;
			throw std::runtime_error("Observer::bracketStringToNumber(...): must be gs or time");
		}


		template<typename SomeModelType>
		FieldType fourPointDelta(size_t i,size_t j,const std::vector<size_t>& gammas,const SomeModelType& model)
		{
			size_t hs = model.hilbertSize(0);
			size_t nx = 0;
			while(hs) {
				hs>>=1;
				nx++;
			}
			nx /= 2;
			size_t site = 0;
			const MatrixType& opC0 = model.naturalOperator("c",site,gammas[0] + 0*nx); // C_{gamma0,up}
			const MatrixType& opC1 = model.naturalOperator("c",site,gammas[1] + 1*nx); // C_{gamma1,down}
			const MatrixType& opC2 = model.naturalOperator("c",site,gammas[2] + 1*nx); // C_{gamma2,down}
			const MatrixType& opC3 = model.naturalOperator("c",site,gammas[3] + 0*nx); // C_{gamma3,up}

			//if (i+1>=n-1 || j+1>=n-1) return 0;
			return fourpoint_(
					'C',i,opC0,
			       'C',i+1,opC1,
			       'N',j,opC2,
			       'N',j+1,opC3,-1);
		}

		ObserverHelperType helper_;
		ConcurrencyType& concurrency_;
		bool verbose_;
		OnePointCorrelationsType onepoint_;
		CorrelationsSkeletonType skeleton_;
		TwoPointCorrelationsType twopoint_;
		FourPointCorrelationsType fourpoint_;
	};  //class Observer
} // namespace Dmrg

/*@}*/
#endif
