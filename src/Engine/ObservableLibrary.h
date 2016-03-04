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
	typedef typename ObserverType::BraketType BraketType;

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
	      observe_(io,numberOfSites-2,hasTimeEvolution,model,verbose)
	{
		PsimagLite::String modelName = model.params().model;
		bool hubbardLike = (modelName == "HubbardOneBand" ||
		                    modelName == "HubbardOneBandExtended");

		if (hasTimeEvolution && hubbardLike) {
			SizeType site = 0;
			// FIXME: No support for site varying operators
			matrixNup_ = model_.naturalOperator("nup",site,0);
			matrixNdown_ = model_.naturalOperator("ndown",site,0);
		}
	}

	bool endOfData() const { return observe_.endOfData(); }

	const ModelType& model() const { return model_; }

	void interpret(const PsimagLite::String& list, SizeType rows, SizeType cols)
	{
		typename BraketType::VectorStringType vecStr;
		PsimagLite::tokenizer(list,vecStr,",");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			BraketType braket(model_, vecStr[i]);

			SizeType threadId = 0;
			if (braket.points() == 1) {
				PreOperatorSiteIndependentType preOperator(braket.op(0),
				                                           braket.opName(0),
				                                           threadId);
				measureOnePoint(braket.bra(),
				                preOperator,
				                braket.ket());
				continue;
			}

			manyPoint(0,braket,rows,cols);
		}
	}

	void measureTriage(const PsimagLite::String& label,
	                   SizeType rows,
	                   SizeType cols,
	                   SizeType orbitals,
	                   bool hasTimeEvolution)
	{
		bool ot = (label == "ot" || label== "time");

		if (hasTimeEvolution && ot) {
			measureTime("superDensity");
			measureTime("nupNdown");
			measureTime("nup+ndown");
			measureTime("sz");
			return;
		}

		if (hasTimeEvolution) setBrakets("time","time");

		const PsimagLite::String& modelName = model_.params().model;

		// Immm supports only onepoint:
		if (modelName=="Immm" && label!="onepoint") {
			PsimagLite::String str(__FILE__);
			str += " "  + ttos(__LINE__) + "\n";
			str += "Model Immm only supports onepoint\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		if (!hasTimeEvolution && label == "onepoint") {
			SizeType numberOfDofs = dofsFromModelName();
			return measureTheOnePoints(numberOfDofs);
		}

		measure(label,rows,cols,orbitals);
	}

private:

	void setBrakets(const PsimagLite::String& left,const PsimagLite::String& right)
	{
		observe_.setBrakets(left,right);
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

			observe_.setBrakets(bra,ket);
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
			observe_.setBrakets(bra,ket);
			tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0,
			                                    opAcorner,
			                                    ApplyOperatorType::BORDER_YES);
			std::cout<<x<<" "<<tmp1;
			std::cout<<" "<<observe_.time(threadId)<<"\n";
		}
	}

	void measure(const PsimagLite::String& label,SizeType rows,SizeType cols,SizeType orbitals)
	{
		PsimagLite::String modelName = model_.params().model;

		// FIXME: No support for site varying operators
		if (label=="cc") {
			BraketType braket(model_,"<gs|c?0-;c'?0-|gs>");
			manyPoint(0,braket,rows,cols); // c_{0,0} spin down
			BraketType braket2(model_,"<gs|c?1-;c'?1-|gs>");
			manyPoint(0,braket2,rows,cols); // c_{0,0} spin down
		} else if (label=="nn") {
			for (SizeType i =0; i < orbitals; ++i) {
				PsimagLite::String str = "<gs|n?" + ttos(i) + ";n?" + ttos(i) + "|gs>";
				BraketType braket(model_,str);
				manyPoint(0,braket,rows,cols);
			}

		} else if (label=="szsz") {
			resizeStorage(szsz_,rows,cols,orbitals);
			MatrixType tSzTotal;
			SizeType counter = 0;
			for (SizeType i = 0; i < orbitals; ++i) {
				for (SizeType j = i; j < orbitals; ++j) {
					PsimagLite::String str = "<gs|z?" + ttos(i) + ";z?" + ttos(j) + "|gs>";
					BraketType braket(model_,str);
					manyPoint(&szsz_[counter],braket,rows,cols);
					MatrixType tSzThis = szsz_[counter];
					RealType factor = (i != j) ? 2.0 : 1.0;
					if (counter == 0)
						tSzTotal =  factor*tSzThis;
					else
						tSzTotal +=  factor*tSzThis;
					if (PsimagLite::Concurrency::root()) {
						std::cout<<"OperatorSz orb"<<i<<"-"<<j<<":\n";
						std::cout<<tSzThis;
					}

					counter++;
				}
			}

			if (PsimagLite::Concurrency::root() && orbitals > 1) {
				std::cout<<"OperatorSz tot:\n";
				std::cout<<tSzTotal;
			}

		} else if (label=="s+s-") {
			// Si^+ Sj^-
			resizeStorage(sPlusSminus_,rows,cols,orbitals);
			MatrixType tSpTotal;
			SizeType counter = 0;
			for (SizeType i = 0; i < orbitals; ++i) {
				for (SizeType j = i; j < orbitals; ++j) {
					PsimagLite::String str = "<gs|+?" + ttos(i) + ";-?" + ttos(j) + "|gs>";
					BraketType braket(model_,str);
					manyPoint(&sPlusSminus_[counter],braket,rows,cols);
					MatrixType tSpThis = sPlusSminus_[counter];
					RealType factor = (i != j) ? 2.0 : 1.0;
					if (counter == 0)
						tSpTotal =  factor*tSpThis;
					else
						tSpTotal +=  factor*tSpThis;
					if (PsimagLite::Concurrency::root()) {
						std::cout<<"OperatorS+S- orb"<<i<<"-"<<j<<":\n";
						std::cout<<tSpThis;
					}

					counter++;
				}
			}

			if (PsimagLite::Concurrency::root() && orbitals > 1) {
				std::cout<<"OperatorS+S- tot:\n";
				std::cout<<tSpTotal;
			}

		} else if (label=="s-s+") {
			// Si^- Sj^+
			resizeStorage(sMinusSplus_,rows,cols,orbitals);
			MatrixType tSmTotal;
			SizeType counter = 0;
			for (SizeType i = 0; i < orbitals; ++i) {
				for (SizeType j = i; j < orbitals; ++j) {
					PsimagLite::String str = "<gs|-?" + ttos(i) + ";+?" + ttos(j) + "|gs>";
					BraketType braket(model_,str);
					manyPoint(&sMinusSplus_[counter],braket,rows,cols);
					MatrixType tSmThis = sMinusSplus_[counter];
					RealType factor = (i != j) ? 2.0 : 1.0;
					if (counter == 0)
						tSmTotal =  factor*tSmThis;
					else
						tSmTotal +=  factor*tSmThis;
					if (PsimagLite::Concurrency::root()) {
						std::cout<<"OperatorS-S+ orb"<<i<<"-"<<j<<":\n";
						std::cout<<tSmThis;
					}

					counter++;
				}
			}

			if (PsimagLite::Concurrency::root() && orbitals > 1) {
				std::cout<<"OperatorS+S- tot:\n";
				std::cout<<tSmTotal;
			}

		} else if (label=="ss") {
			MatrixType spinTotalTotal;
			SizeType counter = 0;
			for (SizeType x = 0; x < orbitals; ++x) {
				for (SizeType y = x; y < orbitals; ++y) {
					if (szsz_.size() == 0)
						measure("szsz",rows,cols,orbitals);
					if (sPlusSminus_.size() == 0)
						measure("s+s-",rows,cols,orbitals);
					if (sMinusSplus_.size() == 0)
						measure("s-s+",rows,cols,orbitals);

					MatrixType spinTotal(szsz_[counter].n_row(),szsz_[counter].n_col());

					RealType factorSpSm = 0.5;
					RealType factorSz = 1.0;
					for (SizeType i=0;i<spinTotal.n_row();i++)
						for (SizeType j=0;j<spinTotal.n_col();j++)
							spinTotal(i,j) = factorSpSm*(
							            sPlusSminus_[counter](i,j) + sMinusSplus_[counter](i,j)) +
							        szsz_[counter](i,j)*factorSz;

					RealType factor = (x != y) ? 2.0 : 1.0;
					if (counter == 0)
						spinTotalTotal =  factor*spinTotal;
					else
						spinTotalTotal +=  factor*spinTotal;

					if (PsimagLite::Concurrency::root()) {
						std::cout<<"SpinTotal orb"<<x<<"-"<<y<<":\n";
						std::cout<<spinTotal;
					}

					counter++;
				}
			}

			if (PsimagLite::Concurrency::root() && orbitals > 1) {
				std::cout<<"SpinTotalTotal:\n";
				std::cout<<spinTotalTotal;
			}

		} else if (label=="dd") {

			BraketType braket(model_,"<gs|d;d'|gs>");
			manyPoint(0,braket,rows,cols);

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
			SparseMatrixType myMatrixSparse(myMatrix);
			typename PsimagLite::Vector<FieldType>::Type result;
			observe_.multiCorrelations(result,myMatrixSparse,rows,cols);
			for (SizeType i=0;i<result.size();i++) {
				std::cout<<i<<" "<<result[i]<<"\n";
			}
		} else {
			PsimagLite::String s = "Unknown label: " + label + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void manyPoint(MatrixType* storage,
	               const BraketType& braket,
	               SizeType rows,
	               SizeType cols)
	{
		if (hasTimeEvolution_) {
			SizeType threadId = 0;
			printSites(threadId);
			std::cout<<"#Time="<<observe_.time(threadId)<<"\n";
		}

		std::cout<<braket.toString()<<"\n";

		if (braket.points() == 2) {
			bool needsPrinting = false;
			if (storage == 0) {
				needsPrinting = true;
				storage = new MatrixType(rows,cols);
			}

			observe_.twoPoint(*storage,braket);

			if (needsPrinting) {
				std::cout<<(*storage);
				delete storage;
				storage = 0;
			}

			return;
		}

		if (braket.points() == 3) {
			observe_.threePoint(braket,rows,cols);
			return;
		}

		assert(braket.points() == 4);
		observe_.fourPoint(braket,rows,cols);
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
			observe_.setBrakets("time","time");
			FieldType superDensity = observe_.template
			        onePoint<ApplyOperatorType>(0,
			                                    opIdentity,
			                                    ApplyOperatorType::BORDER_NO);
			std::cout<<"SuperDensity(Weight of the timeVector)="<<
			           superDensity<<"\n";
		} else if (label=="nupNdown") {
			multiply(A,matrixNup_.data,matrixNdown_.data);
			OperatorType opA(A,1,std::pair<SizeType,SizeType>(0,0),1,su2Related1);
			PreOperatorSiteIndependentType preOperator(opA,"nupNdown",threadId);
			measureOnePoint(preOperator);
		} else if (label=="nup+ndown") {
			A = matrixNup_.data;
			A += matrixNdown_.data;
			OperatorType opA(A,1,std::pair<SizeType,SizeType>(0,0),1,su2Related1);
			PreOperatorSiteIndependentType preOperator(opA,"nup+ndown",threadId);
			measureOnePoint(preOperator);
		} else if (label=="sz") {
			A = matrixNup_.data;
			A += (-1.0)*matrixNdown_.data;
			OperatorType opA(A,1,std::pair<SizeType,SizeType>(0,0),1,su2Related1);
			PreOperatorSiteIndependentType preOperator(opA,"sz",threadId);
			measureOnePoint(preOperator);
		} else {
			PsimagLite::String s = "Unknown label: " + label + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void resizeStorage(VectorMatrixType& v,SizeType rows,SizeType cols,SizeType orbitals)
	{
		if (v.size() != 0) return;
		v.resize(static_cast<SizeType>(orbitals*(orbitals+1)/2));
		for (SizeType i = 0; i < v.size(); ++i)
			v[i].resize(rows,cols);
	}

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
			observe_.setBrakets("gs","gs");
			observe_.setPointer(threadId,i0);

			onePointHookForZero(i0,opA,"gs",threadId);

			FieldType tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0,opA,ApplyOperatorType::BORDER_NO);
			std::cout<<observe_.site(threadId)<<" "<<tmp1;

			if (hasTimeEvolution_) { // for time vector use this one:
				observe_.setBrakets("time","time");

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
				observe_.setBrakets("gs","gs");
				FieldType tmp1 = observe_.template
				        onePoint<ApplyOperatorType>(i0,
				                                    opAcorner,
				                                    ApplyOperatorType::BORDER_YES);
				std::cout<<x<<" "<<tmp1;

				if (hasTimeEvolution_) {// for time vector use this one:
					observe_.setBrakets("time","time");
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

	SizeType dofsFromModelName() const
	{
		const PsimagLite::String& modelName = model_.params().model;
		SizeType site = 0; // FIXME : account for Hilbert spaces changing with site
		SizeType dofs = SizeType(log(model_.hilbertSize(site))/log(2.0));
		std::cerr<<"DOFS= "<<dofs<<" <------------------------------------\n";
		if (modelName.find("FeAsBasedSc")!=PsimagLite::String::npos) return dofs;
		if (modelName.find("FeAsBasedScExtended")!=PsimagLite::String::npos) return dofs;
		if (modelName.find("HubbardOneBand")!=PsimagLite::String::npos) return dofs;

		// max number here, site dependence taken into account elsewhere
		if (modelName.find("Immm")!=PsimagLite::String::npos) return 4;
		return 0;
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
	OperatorType matrixNup_,matrixNdown_;
	VectorMatrixType szsz_,sPlusSminus_,sMinusSplus_;

}; // class ObservableLibrary

} // namespace Dmrg

/*@}*/
#endif // OBSERVABLE_LIBRARY_H

