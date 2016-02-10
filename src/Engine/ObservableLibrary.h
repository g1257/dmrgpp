/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

/*! \file ObservableLibrary.h
 *
 *  A library of most used observables
 */
#ifndef OBSERVABLE_LIBRARY_H
#define OBSERVABLE_LIBRARY_H

#include "Matrix.h" // in PsimagLite
#include "PreOperatorSiteDependent.h"
#include "PreOperatorSiteIndependent.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename ObserverType>
class ObservableLibrary {
public:
	typedef typename ObserverType::ModelType ModelType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef typename ObserverType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType> ApplyOperatorType;
	typedef typename ModelType::RealType RealType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename VectorWithOffsetType::value_type FieldType;
	typedef PsimagLite::Matrix<FieldType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef PreOperatorBase<ModelType> PreOperatorBaseType;
	typedef PreOperatorSiteDependent<ModelType> PreOperatorSiteDependentType;
	typedef PreOperatorSiteIndependent<ModelType> PreOperatorSiteIndependentType;

	template<typename IoInputter>
	ObservableLibrary(
	        IoInputter& io,
	        SizeType numberOfSites,
	        bool hasTimeEvolution,
	        const ModelType& model,
	        bool verbose)
	    : numberOfSites_(numberOfSites),
	      hasTimeEvolution_(hasTimeEvolution),
	      model_(model),
	      observe_(io,numberOfSites-2,hasTimeEvolution,model,verbose),
	      szsz_(0,0)
	{
		PsimagLite::String modelName = model.params().model;
		bool hubbardLike = (modelName == "HubbardOneBand" ||
		                    modelName == "HubbardOneBandExtended");

		if (hasTimeEvolution && hubbardLike) {
			SizeType site = 0;
			// FIXME: No support for site varying operators
			fullMatrixToCrsMatrix(matrixNup_,model_.naturalOperator("nup",site,0));
			fullMatrixToCrsMatrix(matrixNdown_,model_.naturalOperator("ndown",site,0));
		}
	}

	void measure(const PsimagLite::String& label,SizeType rows,SizeType cols)
	{
		SizeType threadId = 0;

		PsimagLite::String str("WARNING: ObservableLibrary: ");
		str += "deprecated use of measure\n";
		std::cerr<<str;

		SizeType site = 0; // FIXME: No support for site varying operators
		if (label=="cc") {
			MatrixType opC = model_.naturalOperator("c",site,0); // c_{0,0} spin up
			MatrixType opCtranspose = transposeConjugate(opC);
			measureOne("OperatorC",opC,"",opCtranspose,-1,rows,cols,threadId);
			MatrixType opC2 = model_.naturalOperator("c",site,1); // c_{0,0} spin down
			MatrixType opCtranspose2 = transposeConjugate(opC2);
			measureOne("OperatorC",opC2,"",opCtranspose2,-1,rows,cols,threadId);
		} else if (label=="nn") {
			MatrixType opN = model_.naturalOperator("n",site,0);
			measureOne("OperatorN",opN,"",opN,1,rows,cols,threadId);
		} else if (label=="szsz") {
			MatrixType Sz = model_.naturalOperator("z",site,0);
			szsz_ = observe_.correlations(Sz,Sz,1,rows,cols);
			if (PsimagLite::Concurrency::root()) {
				std::cout<<"OperatorSz:\n";
				std::cout<<szsz_;
			}
		} else if (label=="s+s-") {
			// Si^+ Sj^-
			const MatrixType& sPlus = model_.naturalOperator("+",site,0);
			MatrixType sPlusT = transposeConjugate(sPlus);
			sPlusSminus_ = observe_.correlations(sPlus,sPlusT,1,rows,cols);
			if (PsimagLite::Concurrency::root()) {
				std::cout<<"OperatorSplus:\n";
				std::cout<<sPlusSminus_;
			}
		} else if (label=="s-s+") {
			// Si^- Sj^+
			const MatrixType& sMinus = model_.naturalOperator("-",site,0);
			MatrixType sMinusT = transposeConjugate(sMinus);
			sMinusSplus_= observe_.correlations(sMinus,sMinusT,1,rows,cols);
			if (PsimagLite::Concurrency::root()) {
				std::cout<<"OperatorSminus:\n";
				std::cout<<sMinusSplus_;
			}
		} else if (label=="ss") {
			if (szsz_.n_row()==0) measure("szsz",rows,cols);
			if (sPlusSminus_.n_row()==0)  measure("s+s-",rows,cols);
			if (sMinusSplus_.n_row()==0)  measure("s-s+",rows,cols);

			MatrixType spinTotal(szsz_.n_row(),szsz_.n_col());

			for (SizeType i=0;i<spinTotal.n_row();i++)
				for (SizeType j=0;j<spinTotal.n_col();j++)
					spinTotal(i,j) = static_cast<RealType>(0.5)*(sPlusSminus_(i,j) +
					                                             sMinusSplus_(i,j)) +
					        szsz_(i,j);

			if (PsimagLite::Concurrency::root()) {
				std::cout<<"SpinTotal:\n";
				std::cout<<spinTotal;
			}
		} else if (label=="dd") {

			const MatrixType& oDelta = model_.naturalOperator("d",site,0);
			MatrixType oDeltaT;
			transposeConjugate(oDeltaT,oDelta);
			measureOne("TWO-POINT DELTA-DELTA^DAGGER",oDelta,"",oDeltaT,1,
			           rows,cols,threadId);
		} else if (label=="dd4") {
			if (model_.geometry().label(0)!="ladderx") {
				PsimagLite::String str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "dd4 only available for ladderx\n";
				throw PsimagLite::RuntimeError(str.c_str());
			}
			for (SizeType g=0;g<16;g++) {
				typename PsimagLite::Vector<SizeType>::Type gammas(4,0);
				gammas[0] = (g & 1);
				gammas[1] = (g & 2)>>1;
				gammas[2] = (g & 4) >> 2;
				gammas[3] = (g & 8) >> 3;
				std::cout<<"#DD4 for the following orbitals: ";
				for (SizeType i=0;i<gammas.size();i++) std::cout<<gammas[i]<<" ";
				std::cout<<"\n";
				MatrixType fpd(numberOfSites_/2,numberOfSites_/2);
				observe_.fourPointDeltas(fpd,gammas,model_);
				for (SizeType i=0;i<fpd.n_row();i++) {
					for (SizeType j=0;j<fpd.n_col();j++) {
						std::cout<<fpd(i,j)<<" ";
					}
					std::cout<<"\n";
				}
			}
		} else if (label == "multi") {
			if (model_.params().model!="HubbardOneBand")
				throw PsimagLite::RuntimeError("multi: not for this model\n");
			MatrixType myMatrix(4,4);
			myMatrix(0,0) = myMatrix(3,3) = 1.0;
			myMatrix(1,1) = myMatrix(2,2) = -1.0;

			typename PsimagLite::Vector<FieldType>::Type result;
			observe_.multiCorrelations(result,myMatrix,rows,cols);
			for (SizeType i=0;i<result.size();i++) {
				std::cout<<i<<" "<<result[i]<<"\n";
			}
		} else {
			PsimagLite::String s = "Unknown label: " + label + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	template<typename SomeBracketType>
	void manyPoint(const SomeBracketType& bracket,
	               SizeType rows,
	               SizeType cols,
	               SizeType threadId)
	{
		if (bracket.points() == 3) {
			throw PsimagLite::RuntimeError("observe 3-point not ready yet\n");
		}

		assert(bracket.points() == 4);

		VectorMatrixType v = observe_.ladder(bracket,rows,cols,threadId);

		for (SizeType i = 0; i < v.size(); ++i) {
			std::cout<<v[0];
			std::cout<<"- - - - - - - - - - - - - \n";
		}
	}

	template<typename SomeBracketType>
	void manyPoint(const SomeBracketType& bracket,
	               SizeType rows,
	               SizeType cols)
	{
		if (bracket.points() == 3) {
			observe_.threePoint(bracket,rows,cols);
			return;
		}

		assert(bracket.points() == 4);
		observe_.fourPoint(bracket,rows,cols);
	}

	void measureTime(const PsimagLite::String& label)
	{
		PsimagLite::String str("WARNING: ObservableLibrary: ");
		str += "deprecated use of measureTime\n";
		std::cerr<<str;

		SparseMatrixType A;
		Su2RelatedType su2Related1;
		SizeType threadId = 0;

		if (label=="superDensity") {
			A.makeDiagonal(model_.hilbertSize(observe_.site(threadId)),1.0);
			std::pair<SizeType,SizeType> zeroZero(0,0);
			OperatorType opIdentity(A,1,zeroZero,1,su2Related1);
			observe_.setBrackets("time","time");
			FieldType superDensity = observe_.template
			        onePoint<ApplyOperatorType>(0,
			                                    opIdentity,
			                                    ApplyOperatorType::BORDER_NO);
			std::cout<<"SuperDensity(Weight of the timeVector)="<<
			           superDensity<<"\n";
		} else if (label=="nupNdown") {
			multiply(A,matrixNup_,matrixNdown_);
			OperatorType opA(A,1,std::pair<SizeType,SizeType>(0,0),1,su2Related1);
			PreOperatorSiteIndependentType preOperator(opA,"nupNdown",threadId);
			measureOnePoint(preOperator);
		} else if (label=="nup+ndown") {
			A = matrixNup_;
			A += matrixNdown_;
			OperatorType opA(A,1,std::pair<SizeType,SizeType>(0,0),1,su2Related1);
			PreOperatorSiteIndependentType preOperator(opA,"nup+ndown",threadId);
			measureOnePoint(preOperator);
		} else if (label=="sz") {
			A = matrixNup_;
			A += (-1)*matrixNdown_;
			OperatorType opA(A,1,std::pair<SizeType,SizeType>(0,0),1,su2Related1);
			PreOperatorSiteIndependentType preOperator(opA,"sz",threadId);
			measureOnePoint(preOperator);
		} else {
			PsimagLite::String s = "Unknown label: " + label + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void measureTheOnePoints(SizeType numberOfDofs)
	{
		PsimagLite::String str("WARNING: ObservableLibrary: ");
		str += "measureTheOnePoints is deprecated\n";
		std::cerr<<str;

		SizeType threadId = 0;
		for (SizeType dof=0;dof<numberOfDofs;dof++) {
			for (SizeType dof2=dof;dof2<numberOfDofs;dof2++) {
				PsimagLite::String str("c^\\dagger(dof=");
				str += ttos(dof) + ") c(dof=" + ttos(dof2) + ")";
				PreOperatorSiteDependentType preOperator(dof,dof2,model_,str,threadId);
				measureOnePoint(preOperator);
			}
		}
	}

	void setBrackets(const PsimagLite::String& left,const PsimagLite::String& right)
	{
		observe_.setBrackets(left,right);
	}

	bool endOfData() const { return observe_.endOfData(); }

	void measureOne(const PsimagLite::String& label1,
	                const PsimagLite::Matrix<FieldType>& op1,
	                const PsimagLite::String& label2,
	                const PsimagLite::Matrix<FieldType>& op2,
	                int fermionSign,
	                SizeType rows,
	                SizeType cols,
	                SizeType threadId)
	{
		if (hasTimeEvolution_) printSites(threadId);

		const MatrixType& v =
		        observe_.correlations(op1,op2,fermionSign,rows,cols);;
		if (PsimagLite::Concurrency::root()) {
			if (hasTimeEvolution_)
				std::cout<<"#Time="<<observe_.time(threadId)<<"\n";
			std::cout<<label1<<label2<<":\n";
			std::cout<<v;
		}
	}

	const ModelType& model() const { return model_; }

	void measureOnePoint(const PsimagLite::String& bra,
	                     const PreOperatorBaseType& preOperator,
	                     const PsimagLite::String& ket)
	{
		SizeType threadId = preOperator.threadId();
		printMarker(threadId);

		for (SizeType i0 = 0;i0<observe_.size();i0++) {
			if (!preOperator.isValid(i0+1)) continue;

			OperatorType opA = preOperator(i0+1);

			preOperator.printMatrix(opA.data,preOperator.siteDependent(),i0);

			if (i0==0) {
				std::cout<<"site <"<<bra<<"|"<<preOperator.label();
				std::cout<<"|"<<ket<<"> time\n";
			}

			observe_.setBrackets(bra,ket);
			observe_.setPointer(threadId,i0);

			onePointHookForZero(i0,opA,"gs",threadId);

			FieldType tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0,opA,ApplyOperatorType::BORDER_NO);
			std::cout<<observe_.site(threadId)<<" "<<tmp1;
			std::cout<<" "<<observe_.time(threadId)<<"\n";

			if (!observe_.isAtCorner(numberOfSites_,threadId)) continue;

			// also calculate next or prev. site:
			SizeType x = (observe_.site(threadId)==1) ? 0 : numberOfSites_-1;

			// operator might be site dependent
			if (!preOperator.isValid(x)) continue;
			OperatorType opAcorner = preOperator(x);

			// do the corner case
			observe_.setBrackets(bra,ket);
			tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0,
			                                    opAcorner,
			                                    ApplyOperatorType::BORDER_YES);
			std::cout<<x<<" "<<tmp1;
			std::cout<<" "<<observe_.time(threadId)<<"\n";
		}
	}

private:

	void measureOnePoint(const PreOperatorBaseType& preOperator)
	{
		SizeType threadId = preOperator.threadId();
		printMarker(threadId);

		for (SizeType i0 = 0;i0<observe_.size();i0++) {
			if (!preOperator.isValid(i0+1)) continue;

			OperatorType opA = preOperator(i0+1);

			preOperator.printMatrix(opA.data,preOperator.siteDependent(),i0);

			if (i0==0) {
				std::cout<<"site "<<preOperator.label()<<"(gs) ";
				if (hasTimeEvolution_)
					std::cout<<preOperator.label()<<"(timevector) time";
				std::cout<<"\n";
			}
			// for g.s. use this one:
			observe_.setBrackets("gs","gs");
			observe_.setPointer(threadId,i0);

			onePointHookForZero(i0,opA,"gs",threadId);

			FieldType tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0,opA,ApplyOperatorType::BORDER_NO);
			std::cout<<observe_.site(threadId)<<" "<<tmp1;

			if (hasTimeEvolution_) { // for time vector use this one:
				observe_.setBrackets("time","time");

				onePointHookForZero(i0,opA,"time",threadId);

				FieldType tmp2 = observe_.template
				        onePoint<ApplyOperatorType>(i0,opA,ApplyOperatorType::BORDER_NO);

				std::cout<<" "<<tmp2<<" "<<observe_.time(threadId);
			}

			std::cout<<"\n";
			// also calculate next or prev. site:
			if (observe_.isAtCorner(numberOfSites_,threadId)) {
				SizeType x = (observe_.site(threadId)==1) ? 0 : numberOfSites_-1;

				// operator might be site dependent
				if (!preOperator.isValid(x)) continue;
				OperatorType opAcorner = preOperator(x);

				// do the corner case
				// for g.s. use this one:
				observe_.setBrackets("gs","gs");
				FieldType tmp1 = observe_.template
				        onePoint<ApplyOperatorType>(i0,
				                                    opAcorner,
				                                    ApplyOperatorType::BORDER_YES);
				std::cout<<x<<" "<<tmp1;

				if (hasTimeEvolution_) {// for time vector use this one:
					observe_.setBrackets("time","time");
					FieldType tmp2 = observe_.template
					        onePoint<ApplyOperatorType>(i0,
					                                    opAcorner,
					                                    ApplyOperatorType::BORDER_YES);
					std::cout<<" "<<tmp2<<" "<<observe_.time(threadId);
				}

				std::cout<<"\n";
			}
		}
	}

	void onePointHookForZero(SizeType i0,
	                         const OperatorType& opA,
	                         const PsimagLite::String& gsOrTime,
	                         SizeType threadId)
	{
		if (hasTimeEvolution_) return;
		if (observe_.site(threadId)!=1 || observe_.isAtCorner(numberOfSites_,threadId))
			return;
		assert(observe_.site(threadId)==1);
		FieldType tmp1 = observe_.template onePointHookForZero<ApplyOperatorType>(i0,opA);
		std::cout<<0<<" "<<tmp1;
		if (hasTimeEvolution_ && gsOrTime=="time") std::cout<<"\n";
		if (!hasTimeEvolution_) std::cout<<"\n";
	}

	void printSites(SizeType threadId)
	{
		printMarker(threadId);
		std::cout<<"#Sites=";
		observe_.setPointer(threadId,0);
		if (observe_.site(threadId)==1) std::cout<<"0 ";
		if (observe_.site(threadId)==numberOfSites_-2)
			std::cout<<(numberOfSites_-1)<<" ";
		for (SizeType i=0;i<observe_.size();i++) {
			observe_.setPointer(threadId,i);
			SizeType x = observe_.site(threadId);
			std::cout<<x<<" ";
		}

		if (observe_.site(threadId)==1) std::cout<<"0";
		if (observe_.site(threadId) ==numberOfSites_-2) {
			std::cout<<(numberOfSites_-1);
		}

		std::cout<<"\n";
	}

	void printMarker(SizeType threadId) const
	{
		if (!hasTimeEvolution_) return;
		SizeType marker = observe_.marker(threadId);
		PsimagLite::String s = "INVALID MARKER";
		switch (marker) {
		case 0:
			s=" NOT ALL OPERATORS APPLIED YET";
			break;
		case 1:
			s=" ALL OPERATORS HAVE BEEN APPLIED";
		}

		std::cout<<s<<"\n";
	}

	SizeType numberOfSites_;
	bool hasTimeEvolution_;
	const ModelType& model_; // not the owner
	ObserverType observe_;
	SparseMatrixType matrixNup_,matrixNdown_;
	MatrixType szsz_,sPlusSminus_,sMinusSplus_;

}; // class ObservableLibrary

} // namespace Dmrg

/*@}*/
#endif // OBSERVABLE_LIBRARY_H

