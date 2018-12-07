/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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
#include "Vector.h"

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
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename VectorWithOffsetType::value_type FieldType;
	typedef typename PsimagLite::Vector<FieldType>::Type VectorFieldType;
	typedef PsimagLite::Matrix<FieldType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef PreOperatorBase<ModelType> PreOperatorBaseType;
	typedef PreOperatorSiteDependent<ModelType> PreOperatorSiteDependentType;
	typedef PreOperatorSiteIndependent<ModelType> PreOperatorSiteIndependentType;
	typedef typename ObserverType::BraketType BraketType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

	template<typename IoInputter>
	ObservableLibrary(IoInputter& io,
	                  SizeType numberOfSites,
	                  bool hasTimeEvolution,
	                  const ModelType& model,
	                  SizeType start,
	                  SizeType nf,
	                  SizeType trail,
	                  bool verbose)
	    : numberOfSites_(numberOfSites),
	      hasTimeEvolution_(hasTimeEvolution),
	      model_(model),
	      observe_(io, start, nf, trail, hasTimeEvolution, model, verbose)
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
		PsimagLite::split(vecStr, list, ",");

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

	void setBrakets(const PsimagLite::String& left,
	                const PsimagLite::String& right)
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
		VectorFieldType density;

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

			SizeType tmp = density.size();
			onePointHookForZero(i0,opA,"gs",threadId,density);
			SizeType lastOne = density.size();
			if (tmp != density.size() && lastOne > 0) {
				lastOne--;
				std::cout<<"0 "<<density[lastOne];
				std::cout<<" "<<observe_.time(threadId)<<"\n";
			}

			FieldType tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0,opA,ApplyOperatorType::BORDER_NO);
			std::cout<<observe_.site(threadId)<<" "<<tmp1;
			std::cout<<" "<<observe_.time(threadId)<<"\n";
			density.push_back(tmp1);

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
			density.push_back(tmp1);
		}
	}

	MatrixType SliceOrbital(const MatrixType& m,
	                        const SizeType o1,
	                        const SizeType o2)
	{
		SizeType orbitals = 2;
		SizeType nsite = numberOfSites_/orbitals;
		MatrixType out(nsite,nsite);
		for (SizeType i = 0; i < nsite; ++i) {
			for (SizeType j = i; j < nsite; ++j) {
				SizeType k = i*orbitals + o1;
				SizeType l = j*orbitals + o2;
				out(i,j) = m(k,l);
			}
		}

		std::cout << out;
		return out;
	}

	void measure(const PsimagLite::String& label,
	             SizeType rows,
	             SizeType cols,
	             SizeType orbitals)
	{
		// FIXME: No support for site varying operators
		if (label=="cc") {
			BraketType braket(model_,"<gs|c?0-;c'?0-|gs>");
			manyPoint(0,braket,rows,cols); // c_{0,0} spin down
			BraketType braket2(model_,"<gs|c?1-;c'?1-|gs>");
			manyPoint(0,braket2,rows,cols); // c_{0,0} spin down
		} else if (label=="nn") {
			MatrixType out(rows,cols);
			int fermionicSign = 1;
			SizeType site = 1;
			for (SizeType i = 0; i < orbitals*2; ++i) {
				for (SizeType j = i; j < orbitals*2; ++j) {
					SparseMatrixType O2,O4,n1,n2;
					SparseMatrixType O1 = model_.naturalOperator("c",site,i).data; // c_i
					transposeConjugate(O2,O1); // O2 = transpose(O1)
					SparseMatrixType O3 = model_.naturalOperator("c",site,j).data; // c_j
					transposeConjugate(O4,O3); // O4 = transpose(O3)

					multiply(n1,O2,O1); // c_i^{\dagger}.c_i
					multiply(n2,O4,O3); // c_j^{\dagger}.c_j

					PsimagLite::String str = "<gs|n?" + ttos(i) + ";n?" + ttos(j) + "|gs>";
					observe_.twoPoint(out,n1,n2,fermionicSign);
					std::cout << str << std::endl;
					std::cout << out;
				}
			}

		} else if (label=="szsz") {
			resizeStorage(szsz_,rows,cols,orbitals);
			MatrixType tSzTotal;
			SizeType counter = 0;
			for (SizeType i = 0; i < orbitals; ++i) {
				for (SizeType j = i; j < orbitals; ++j) {
					PsimagLite::String str = "<gs|sz?" + ttos(i) + ";sz?" + ttos(j) + "|gs>";
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
					PsimagLite::String str = "<gs|splus?" + ttos(i) + ";sminus?" + ttos(j) + "|gs>";
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
						std::cout<<tSpTotal;
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
					PsimagLite::String str = "<gs|sminus?" + ttos(i) + ";splus?" + ttos(j) + "|gs>";
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
						std::cout<<tSmTotal;
					}
					counter++;
				}
			}

			if (PsimagLite::Concurrency::root() && orbitals > 1) {
				std::cout<<"OperatorS-S+ tot:\n";
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

		} else if (label == "pp") {
			if (model_.params().model!="TjMultiOrb" &&
			        model_.params().model!="HubbardOneBandExtendedSuper") {
				throw PsimagLite::RuntimeError("pp: not for this model\n");
			}

			typename PsimagLite::Vector<MatrixType*>::Type results;
			typename PsimagLite::Vector<PsimagLite::String>::Type names;
			ppTwopoint(results,names,rows,cols);
			ppFourpoint(results,names,rows,cols);

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
				std::cout<<"DD4 for the following orbitals: ";
				for (SizeType i=0;i<gammas.size();i++) std::cout<<gammas[i]<<" ";
				std::cout<<"\n";
				MatrixType fpd(numberOfSites_/2,numberOfSites_/2);
				observe_.fourPointDeltas(fpd,gammas,model_);
				std::cout<<fpd;
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
			for (SizeType i=0;i<result.size();i++)
				std::cout<<i<<" "<<result[i]<<"\n";

		} else if (label == "ddOrbitals") {
			if (model_.params().model!="HubbardOneBandExtendedSuper") {
				throw PsimagLite::RuntimeError("pp: not for this model\n");
			}
			typename PsimagLite::Vector<MatrixType*>::Type results;
			typename PsimagLite::Vector<PsimagLite::String>::Type names;
			ddOrbitalsTwopoint(results,names,rows,cols);
			ddOrbitalsFourpoint(results,names,rows,cols);
		} else {
			PsimagLite::String s = "Unknown label: " + label + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void ppTwopoint(typename PsimagLite::Vector<MatrixType*>::Type& result,
	                typename PsimagLite::Vector<PsimagLite::String>::Type& names,
	                SizeType rows,
	                SizeType cols)
	{
		// Two-point Pair
		MatrixType m1(rows,cols);
		MatrixType m2(rows,cols);
		ppTwo(m1,m2,0);

		if (model_.params().model=="HubbardOneBandExtendedSuper") {
			rows = rows/2;   // Actually: divided by # of orbitals
			cols = cols/2;
		}

		m1.clear(); m1.resize(rows,cols);
		m2.clear(); m2.resize(rows,cols);
		ppTwo(m1,m2,1);
		std::cout << "PairPair Correlations S^{lu}_{on}" << std::endl;
		std::cout << m1;
		std::cout << "PairPair Correlations T^{lu}_{on}" << std::endl;
		std::cout << m2;

		m1.clear(); m1.resize(rows,cols);
		m2.clear(); m2.resize(rows,cols);
		names.push_back("T_{on}");
		ppTwo(m1,m2,3);
	}

	void ppTwo(MatrixType& m,MatrixType& m2, SizeType flag)
	{
		int fermionicSign = 1; //bosons
		SizeType site = 1;
		SizeType orbitals = logBase2(model_.hilbertSize(site));
		assert(!(orbitals & 1));
		orbitals /= 2;

		if (flag==0) {
			fermionicSign = 1; //bosons
			SizeType spin0 = 0; // up
			SizeType spin1 = 1; // down

			// c_dn,0
			SparseMatrixType O1 = model_.naturalOperator("c",site,spin1).data;
			// c_up,0
			SparseMatrixType O2 = model_.naturalOperator("c",site,spin0).data;

			SparseMatrixType A,B;
			multiply(B,O1,O2); // c_dn,0 . c_up,0.
			transposeConjugate(A,B);
			observe_.twoPoint(m,A,B,fermionicSign);
			//std::cout << m;
			std::cout << "PairPair Correlations S^{l}_{on}" << std::endl;
			SliceOrbital(m,0,0);
			std::cout << "PairPair Correlations S^{u}_{on}" << std::endl;
			SliceOrbital(m,1,1);

		} else if (flag==1) {
			SizeType orb1 = 0;
			SizeType orb2 = 1;
			SizeType orb3 = 1;
			SizeType orb4 = 0;
			int sign = -1;
			PsimagLite::String onsiteOrNot = "two";
			// notice - orb3 and orb4 order had to be fliped to preserve
			// i1 > i2 > i3 > i4 thin site ordering, this should add multiplication by (-1.0);
			//std::cout << "PairPair Correlations S^{lu}_{on}" << std::endl;
			ppFour(m,m2,orb1,orb2,orb4,orb3,onsiteOrNot,sign);

		} else if (flag==2) {
			SizeType orb1 = 0;
			SizeType orb2 = 1;
			SizeType orb3 =	1;
			SizeType orb4 = 0;
			int sign = 1;
			PsimagLite::String onsiteOrNot = "two";
			// notice - orb3 and orb4 order had to be fliped to preserve
			// i1 > i2 > i3 > i4 ordering, this adds multiplication by (-1.0);
			std::cout << "PairPair Correlations T^{lu}_{on}" << std::endl;
			ppFour(m,m2,orb1,orb2,orb4,orb3,onsiteOrNot,sign);
		} else if (flag==3) {
			SizeType orb1 = 0;
			SizeType orb2 = 1;
			SizeType orb3 = 0;
			SizeType orb4 = 1;
			ppupupdndn(m, m2, orb1, orb2, orb3, orb4);

			std::cout << "PairPair Correlations T^{ab-upup}_{on}" << std::endl;
			std::cout << m;
			std::cout << "PairPair Correlations T^{ab-dndn}_{on}" << std::endl;
			std::cout << m2;
		} else {
			PsimagLite::String s = "Unknown flag: " + ttos(flag);
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void ppupupdndn(MatrixType& m,
	                MatrixType& m2,
	                SizeType orb1,
	                SizeType orb2,
	                SizeType orb3,
	                SizeType orb4) const
	{
		PsimagLite::String string = "two";
		if (string!="four" && string!="two") {
			throw PsimagLite::RuntimeError("ppFour: only string = 'two' or 'four' is allowed \n");
		}
		SizeType rows = m.n_row();
		SizeType cols = m.n_col();
		assert(rows == cols);
		typename PsimagLite::Vector<PairSizeType>::Type pairs;
		int sign = 1;
		VectorSizeType gammas(1,1+sign);
		SizeType orbitals = 2;
		SizeType bigSize = rows*orbitals*rows*orbitals*2;
		m.resize(bigSize, bigSize, static_cast<typename MatrixType::value_type>(0.0));

		SizeType offset = (string=="four") ? orbitals : 1;
		SizeType jmax = (string=="four") ? cols-1 : cols;
		for (SizeType i = 0; i < rows; ++i) {
			SizeType thini1 = i*orbitals + orb1;
			SizeType thini2 = (string=="four") ? (i+1)*orbitals + orb2 : i*orbitals+orb2;
			for (SizeType j = i + offset; j < jmax; ++j) {
				SizeType thinj1 = j*orbitals + orb3;
				SizeType thinj2 = (string=="four") ? (j+1)*orbitals + orb4 : j*orbitals + orb4;
				for (SizeType spin0 = 0; spin0 < 2; ++spin0) {
					SizeType spin1=spin0;
					pairs.push_back(PairSizeType(thini1+thini2*rows*orbitals+rows*orbitals*rows*orbitals*spin0,
					                             thinj1+thinj2*rows*orbitals+rows*orbitals*rows*orbitals*spin1));
				}
			}
		}

		typedef typename ObserverType::Parallel4PointDsType Parallel4PointDsType;
		typedef PsimagLite::Parallelizer<Parallel4PointDsType> ParallelizerType;
		ParallelizerType threaded4PointDs(PsimagLite::Concurrency::codeSectionParams);

		Parallel4PointDsType helper4PointDs(m,
		                                    observe_.fourpoint(),
		                                    model_,
		                                    gammas,
		                                    pairs,
		                                    Parallel4PointDsType::MODE_THINupdn);

		threaded4PointDs.loopCreate(helper4PointDs);

		MatrixType mup(rows,cols);
		MatrixType mdown(rows,cols);
		for (SizeType i = 0; i < rows; ++i) {
			SizeType thini1 = i*orbitals + orb1;
			SizeType thini2 = (string=="four") ? (i+1)*orbitals + orb2 : i*orbitals+orb2;
			for (SizeType j = i + offset; j < jmax; ++j) {
				SizeType thinj1 = j*orbitals + orb3;
				SizeType thinj2 = (string=="four") ? (j+1)*orbitals + orb4 : j*orbitals + orb4;

				RealType vsign = -1.0;
				SizeType spin0=0;
				SizeType spin1=0;
				mup(i,j) += vsign*m(thini1+thini2*rows*orbitals+spin0*rows*orbitals*rows*orbitals,
				                    thinj1+thinj2*rows*orbitals+spin1*rows*orbitals*rows*orbitals);

				spin0=1;
				spin1=1;
				mdown(i,j) += vsign*m(thini1+thini2*rows*orbitals+spin0*rows*orbitals*rows*orbitals,
				                      thinj1+thinj2*rows*orbitals+spin1*rows*orbitals*rows*orbitals);
			}
		}

		m = mup;
		m2 = mdown;
		//std::cout << mSinglet;
		//std::cout << mTriplet;

	}

	void ddOrbitalsTwopoint(typename PsimagLite::Vector<MatrixType*>::Type& result,
	                        typename PsimagLite::Vector<PsimagLite::String>::Type& names,
	                        SizeType rows,
	                        SizeType cols)
	{
		// Two-point Pair
		MatrixType* m1 = new MatrixType(rows,cols);
		names.push_back("S^{l}_{on}");
		std::cout << "PairPair Correlations S^{l}_{on}" << std::endl;
		ddOrbitalsTwo(*m1,0);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("S^{u}_{on}");
		std::cout << "PairPair Correlations S^{u}_{on}" << std::endl;
		ddOrbitalsTwo(*m1,1);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("S^{lu}_{on}");
		std::cout << "PairPair Correlations S^{lu}_{on}" << std::endl;
		ddOrbitalsTwo(*m1,2);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{lu}_{on}" << std::endl;
		ddOrbitalsTwo(*m1,3);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{up-up}_{on}" << std::endl;
		ddOrbitalsTwo(*m1,4);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{dn-dn}_{on}" << std::endl;
		ddOrbitalsTwo(*m1,5);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{up*up+dn*dn}_{on}" << std::endl;
		ddOrbitalsTwo(*m1,6);
		result.push_back(m1);
	}

	void ddOrbitalsFourpoint(typename PsimagLite::Vector<MatrixType*>::Type& result,
	                         typename PsimagLite::Vector<PsimagLite::String>::Type& names,
	                         SizeType rows,
	                         SizeType cols) const
	{
		// Singlet four-points
		MatrixType* m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations S^{l}_{nn}" << std::endl;
		names.push_back("S^{l}_{nn}");
		ddOrbitalsFour(*m1,0,0,0,0,-1);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations S^{u}_{nn}" << std::endl;
		names.push_back("S^{u}_{nn}");
		ddOrbitalsFour(*m1,1,1,1,1,-1);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations S^{lu}_{nn}" << std::endl;
		names.push_back("S^{lu}_{nn}");
		ddOrbitalsFour(*m1,0,1,1,0,-1);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations S^{l-u}_{nn}" << std::endl;
		names.push_back("S^{l-u}_{nn}");
		ddOrbitalsFour(*m1,0,0,1,1,-1);
		result.push_back(m1);

		// Triplet four-points
		m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations T^{l}_{nn}" << std::endl;
		names.push_back("T^{l}_{nn}");
		ddOrbitalsFour(*m1,0,0,0,0,1);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations T^{u}_{nn}" << std::endl;
		names.push_back("T{^{u}_{nn}");
		ddOrbitalsFour(*m1,1,1,1,1,1);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations T^{lu}_{nn}" << std::endl;
		names.push_back("T^{lu}_{nn}");
		ddOrbitalsFour(*m1,0,1,1,0,1);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		std::cout << "PairPair Correlations T^{l-u}_{nn}" << std::endl;
		names.push_back("T^{l-u}_{nn}");
		ddOrbitalsFour(*m1,0,0,1,1,1);
		result.push_back(m1);
	}

	void ddOrbitalsTwo(MatrixType& m, SizeType flag)
	{
		int fermionicSign = 1;
		SizeType site = 1;
		SizeType orbitals = logBase2(model_.hilbertSize(site));
		assert(!(orbitals & 1));
		orbitals /= 2;

		if (flag==0) {
			SizeType orb1 = 0;  // lower orbital
			SizeType spin0 = 0; // up
			SizeType spin1 = 1; // down
			// c_dn,0
			SparseMatrixType O1 = model_.naturalOperator("c",site,orb1+spin1*orbitals).data;
			// c_up,0
			SparseMatrixType O2 = model_.naturalOperator("c",site,orb1+spin0*orbitals).data;

			SparseMatrixType A,B;
			multiply(B,O1,O2); // c_dn,0 . c_up,0.
			transposeConjugate(A,B);
			observe_.twoPoint(m,A,B,fermionicSign);

			std::cout << m;
		} else if (flag==1) {
			SizeType orb1 = 1;  // lower orbital
			SizeType spin0 = 0; // up
			SizeType spin1 = 1; // down
			// c_dn,0
			SparseMatrixType O1 = model_.naturalOperator("c",site,orb1+spin1*orbitals).data;
			// c_up,0
			SparseMatrixType O2 = model_.naturalOperator("c",site,orb1+spin0*orbitals).data;
			SparseMatrixType A,B;
			multiply(B,O1,O2); // c_dn,0 . c_up,0.
			transposeConjugate(A,B);
			observe_.twoPoint(m,A,B,fermionicSign);
			std::cout << m;
		} else if (flag==2) {
			SizeType spin0 = 0; // up
			SizeType spin1 = 1; // down
			SparseMatrixType O1 = model_.naturalOperator("c",site,1+spin1*orbitals).data; // c_dn,1
			SparseMatrixType O2 = model_.naturalOperator("c",site,0+spin0*orbitals).data; // c_up,0
			SparseMatrixType O3 = model_.naturalOperator("c",site,1+spin0*orbitals).data; // c_up,1
			SparseMatrixType O4 = model_.naturalOperator("c",site,0+spin1*orbitals).data; // c_dn,0

			SparseMatrixType A,B,tmp1,tmp2;
			multiply(tmp1,O1,O2); // c_dn,1 . c_up,0
			multiply(tmp2,O3,O4); // c_up,1 . c_dn,0

			FieldType mult1, mult2;
			mult1 = 1.0; mult2 = -1.0;
			operatorPlus(B,tmp1,mult1,tmp2,mult2); // B = 1.0*tmp1 + (-1.0)*tmp2 = Singlet
			transposeConjugate(A,B); // A = transpose(B)
			observe_.twoPoint(m,A,B,fermionicSign);
			std::cout << m;
		} else if (flag==3) {
			SizeType spin0 = 0; // up
			SizeType spin1 = 1; // down
			SparseMatrixType O1 = model_.naturalOperator("c",site,1+spin1*orbitals).data; // c_dn,1
			SparseMatrixType O2 = model_.naturalOperator("c",site,0+spin0*orbitals).data; // c_up,0
			SparseMatrixType O3 = model_.naturalOperator("c",site,1+spin0*orbitals).data; // c_up,1
			SparseMatrixType O4 = model_.naturalOperator("c",site,0+spin1*orbitals).data; // c_dn,0

			SparseMatrixType A,B,tmp1,tmp2;
			multiply(tmp1,O1,O2); // c_dn,1 . c_up,0
			multiply(tmp2,O3,O4); // c_up,1 . c_dn,0

			FieldType mult1, mult2;
			mult1 = 1.0; mult2 = 1.0;
			operatorPlus(B,tmp1,mult1,tmp2,mult2); // B = 1.0*tmp1 + (1.0)*tmp2 = Triplet
			transposeConjugate(A,B); // A = transpose(B)
			observe_.twoPoint(m,A,B,fermionicSign);
			std::cout << m;
		} else if (flag==4) {
			SizeType orb0 = 0;  // lower orbital
			SizeType orb1 = 1;  // upper orbital
			SizeType spin0 = 0; // up
			// c_up,0
			SparseMatrixType O1 = model_.naturalOperator("c",site,orb0+spin0*orbitals).data;
			// c_up,1
			SparseMatrixType O2 = model_.naturalOperator("c",site,orb1+spin0*orbitals).data;
			SparseMatrixType A,B;
			multiply(B,O1,O2);      // c_up,0 . c_up,1
			transposeConjugate(A,B);
			observe_.twoPoint(m,A,B,fermionicSign);
			std::cout << m;
		} else if (flag==5) {
			SizeType orb0 = 0;  // lower orbital
			SizeType orb1 = 1;  // upper orbital
			SizeType spin1 = 1; // dn
			// c_dn,0
			SparseMatrixType O1 = model_.naturalOperator("c",site,orb0+spin1*orbitals).data;
			// c_dn,1
			SparseMatrixType O2 = model_.naturalOperator("c",site,orb1+spin1*orbitals).data;
			SparseMatrixType A,B;
			multiply(B,O1,O2);      // c_dn,0 . c_dn,1
			transposeConjugate(A,B);
			observe_.twoPoint(m,A,B,fermionicSign);
			std::cout << m;
		} else if (flag==6) {
			SizeType spin0 = 0; // up
			SizeType spin1 = 1; // down
			SparseMatrixType O1 = model_.naturalOperator("c",site,1+spin1*orbitals).data; // c_up,0
			SparseMatrixType O2 = model_.naturalOperator("c",site,0+spin0*orbitals).data; // c_up,1
			SparseMatrixType O3 = model_.naturalOperator("c",site,1+spin0*orbitals).data; // c_dn,0
			SparseMatrixType O4 = model_.naturalOperator("c",site,0+spin1*orbitals).data; // c_dn,1

			SparseMatrixType A,B,tmp1,tmp2;
			multiply(tmp1,O1,O2); // c_dn,1 . c_up,0
			multiply(tmp2,O3,O4); // c_up,1 . c_dn,0

			FieldType mult1, mult2;
			mult1 = 1.0; mult2 = 1.0;
			// B = 1.0*tmp1 + (1.0)*tmp2 = Triplet = up*up + dn*dn
			operatorPlus(B,tmp1,mult1,tmp2,mult2);
			transposeConjugate(A,B); // A = transpose(B)
			observe_.twoPoint(m,A,B,fermionicSign);
			std::cout << m;
		} else {
			PsimagLite::String s = "Unknown flag: " + ttos(flag);
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void ddOrbitalsFour(MatrixType& m,
	                    SizeType orb1,
	                    SizeType orb2,
	                    SizeType orb3,
	                    SizeType orb4,
	                    int sign) const
	{
		SizeType rows = m.n_row();
		SizeType cols = m.n_row();

		for (SizeType i = 0; i < rows; ++i) {
			for (SizeType j = i + 2; j < cols; ++j) {
				m(i,j) = ddOrbitalsFour2(i,j,orb1,orb2,orb3,orb4,sign);
			}
		}
		std::cout << m;
	}

	FieldType ddOrbitalsFour2(SizeType i,
	                          SizeType j,
	                          SizeType orb1,
	                          SizeType orb2,
	                          SizeType orb3,
	                          SizeType orb4,
	                          int sign) const
	{
		SizeType i1 = i;
		SizeType i2 = i + 1;
		SizeType j1 = j;
		SizeType j2 = j + 1;

		SizeType threadId = 0;
		FieldType sum = 0.0;
		SizeType site = 0;
		SizeType orbitals = logBase2(model_.hilbertSize(site));
		assert(!(orbitals & 1));
		orbitals /= 2;

		for (SizeType spin0 = 0; spin0 < 2; ++spin0) {
			// c(i1,orb1,spin0)
			PsimagLite::String str("<gs|c[" + ttos(site) + "]?" +
			                       ttos(orb1+spin0*orbitals) + ";");

			// c(i2,orb2,1-spin0)
			str += "c[" + ttos(site) + "]?" + ttos(orb2+(1-spin0)*orbitals) + ";";

			for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
				// c(i2,orb2,spin1)
				PsimagLite::String str2("c[" + ttos(site) + "]?" +
				                        ttos(orb3+spin1*orbitals) + "';");

				// c(i3,orb1,1-spin1)
				str2 += "c[" + ttos(site) + "]?" + ttos(orb4+(1-spin1)*orbitals) + "'|gs>";

				BraketType braket(model_, str + str2);
				SizeType val = spin0 + spin1 + 1;
				int signTerm = (val & 1) ? sign : 1;
				sum +=  signTerm*observe_.fourpoint()(i1,i2,j1,j2,braket,threadId);
			}
		}

		return sum;
	}

	void ppFourpoint(typename PsimagLite::Vector<MatrixType*>::Type& result,
	                 typename PsimagLite::Vector<PsimagLite::String>::Type& names,
	                 SizeType rows,
	                 SizeType cols) const
	{
		/*	Alberto Nocera
		 *
//		// Singlet four-points
//		MatrixType* m1 = new MatrixType(rows,cols);
//		std::cout << "PairPair Correlations SrSr_{nn}" << std::endl;
//		names.push_back("S^{l}_{nn}");
//		ppFour(*m1,0,0,0,0,0,-1);
//		result.push_back(m1);

//		m1 = new MatrixType(rows,cols);
//		std::cout << "PairPair Correlations SrSl_{nn}" << std::endl;
//		names.push_back("S^{u}_{nn}");
//		ppFour(*m1,0,0,0,0,1,-1);
//		result.push_back(m1);
		*
		*/
		if (model_.params().model=="HubbardOneBandExtendedSuper") {
			rows = rows/2;   // Actually: divided by # of orbitals
			cols = cols/2;
		}

		PsimagLite::String onsiteOrNot = "four";
		// Singlet four-points
		MatrixType m1(rows,cols);
		MatrixType m2(rows,cols);
//		std::cout << "PairPair Correlations S^{l}_{nn}" << std::endl;
//		names.push_back("S^{l}_{nn}");
//		ppFour(m1,m2,0,0,0,0,onsiteOrNot,-1);


//		m1.clear();
//		m1.resize(rows,cols);
//		std::cout << "PairPair Correlations S^{u}_{nn}" << std::endl;
//		names.push_back("S^{u}_{nn}");
//		ppFour(m1,m2,1,1,1,1,onsiteOrNot,-1);


		m1.clear();
		m2.clear();
		m1.resize(rows,cols);
		m2.resize(rows,cols);
		ppFour(m1,m2,0,1,0,1,onsiteOrNot,-1);
		std::cout << "PairPair Correlations S^{lu}_{nn}" << std::endl;
		std::cout << m1;
		std::cout << "PairPair Correlations T^{lu}_{nn}" << std::endl;
		std::cout << m2;

		m1.clear();
		m2.clear();
		m1.resize(rows,cols);
		m2.resize(rows,cols);
		ppFour(m1,m2,0,1,1,0,onsiteOrNot,-1);
		std::cout << "PairPair Correlations S_RL^{lu}_{nn}" << std::endl;
		std::cout << m1;
		std::cout << "PairPair Correlations T_RL^{lu}_{nn}" << std::endl;
		std::cout << m2;


//		m1 = new MatrixType(rows,cols);
//		std::cout << "PairPair Correlations S^{l-u}_{nn}" << std::endl;
//		names.push_back("S^{l-u}_{nn}");
//		ppFour(*m1,0,0,1,1,onsiteOrNot,-1);
//		result.push_back(m1);

		// Triplet four-points
//		m1 = new MatrixType(rows,cols);
//		std::cout << "PairPair Correlations T^{l}_{nn}" << std::endl;
//		names.push_back("T^{l}_{nn}");
//		ppFour(*m1,0,0,0,0,onsiteOrNot,1);
//		result.push_back(m1);

//		m1 = new MatrixType(rows,cols);
//		std::cout << "PairPair Correlations T^{u}_{nn}" << std::endl;
//		names.push_back("T^{u}_{nn}");
//		ppFour(*m1,1,1,1,1,onsiteOrNot,1);
//		result.push_back(m1);

//		m1 = new MatrixType(rows,cols);
//		std::cout << "PairPair Correlations T^{l-u}_{nn}" << std::endl;
//		names.push_back("T^{l-u}_{nn}");
//		ppFour(*m1,0,0,1,1,onsiteOrNot,1);
//		result.push_back(m1);
	}

	void ppFour(MatrixType& m, MatrixType& m2,
	            SizeType orb1,
	            SizeType orb2,
	            SizeType orb3,
	            SizeType orb4,
	            const PsimagLite::String& string,
	            int sign) const
	{
		if (string!="four" && string!="two") {
			throw PsimagLite::RuntimeError("ppFour: only string = 'two' or 'four' is allowed \n");
		}
		SizeType rows = m.n_row();
		SizeType cols = m.n_col();
		assert(rows == cols);
		typename PsimagLite::Vector<PairSizeType>::Type pairs;
		VectorSizeType gammas(1,1+sign);
		SizeType orbitals = 2;
		SizeType bigSize = rows*orbitals*rows*orbitals*2;
		m.resize(bigSize, bigSize, static_cast<typename MatrixType::value_type>(0.0));

		SizeType offset = (string=="four") ? orbitals : 1;
		SizeType jmax = (string=="four") ? cols-1 : cols;
		for (SizeType i = 0; i < rows; ++i) {
			SizeType thini1 = i*orbitals + orb1;
			SizeType thini2 = (string=="four") ? (i+1)*orbitals + orb2 : i*orbitals+orb2;
			for (SizeType j = i + offset; j < jmax; ++j) {
				SizeType thinj1 = j*orbitals + orb3;
				SizeType thinj2 = (string=="four") ? (j+1)*orbitals + orb4 : j*orbitals + orb4;
				for (SizeType spin0 = 0; spin0 < 2; ++spin0) {
					for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
						pairs.push_back(PairSizeType(thini1+thini2*rows*orbitals+rows*orbitals*
						                             rows*orbitals*spin0,
						                             thinj1+thinj2*rows*orbitals+
						                             rows*orbitals*rows*orbitals*spin1));
					}
				}
			}
		}

		typedef typename ObserverType::Parallel4PointDsType Parallel4PointDsType;
		typedef PsimagLite::Parallelizer<Parallel4PointDsType> ParallelizerType;
		ParallelizerType threaded4PointDs(PsimagLite::Concurrency::codeSectionParams);

		Parallel4PointDsType helper4PointDs(m,
		                                    observe_.fourpoint(),
		                                    model_,
		                                    gammas,
		                                    pairs,
		                                    Parallel4PointDsType::MODE_THIN);

		threaded4PointDs.loopCreate(helper4PointDs);

		MatrixType mTriplet(rows,cols);
		MatrixType mSinglet(rows,cols);
		for (SizeType i = 0; i < rows; ++i) {
			SizeType thini1 = i*orbitals + orb1;
			SizeType thini2 = (string=="four") ? (i+1)*orbitals + orb2 : i*orbitals+orb2;
			for (SizeType j = i + offset; j < jmax; ++j) {
				SizeType thinj1 = j*orbitals + orb3;
				SizeType thinj2 = (string=="four") ? (j+1)*orbitals + orb4 : j*orbitals + orb4;
				for (SizeType spin0 = 0; spin0 < 2; ++spin0) {
					for (SizeType spin1 = 0; spin1 < 2; ++spin1) {

						RealType TripletSign = (spin0==spin1) ? -1.0:-1.0;
						mTriplet(i,j) += TripletSign*m(thini1+thini2*rows*orbitals+
						                               spin0*rows*orbitals*rows*orbitals,
						                               thinj1+thinj2*rows*orbitals+
						                               spin1*rows*orbitals*rows*orbitals);

						RealType SingletSign = (spin0==spin1) ? -1.0:1.0;
						mSinglet(i,j) += SingletSign*m(thini1+thini2*rows*orbitals+
						                               spin0*rows*orbitals*rows*orbitals,
						                               thinj1+thinj2*rows*orbitals+
						                               spin1*rows*orbitals*rows*orbitals);
					}
				}
			}
		}

		m = mSinglet;
		m2 = mTriplet;
		//std::cout << mSinglet;
		//std::cout << mTriplet;

	}

	FieldType ppFour2(SizeType i1,
	                  SizeType i2,
	                  SizeType j1,
	                  SizeType j2,
	                  SizeType orb1,
	                  SizeType orb2,
	                  SizeType orb3,
	                  SizeType orb4,
	                  int sign) const
	{

		SizeType thini1 = i1*2 + orb1;
		SizeType thini2 = i2*2 + orb2;
		SizeType thinj1 = j1*2 + orb3;
		SizeType thinj2 = j2*2 + orb4;

		int fermionicSign = -1;
		SizeType threadId = 0;
		FieldType sum = 0.0;
		SizeType site = 0;
		//		SizeType orbitals = 1;

		for (SizeType spin0 = 0; spin0 < 2; ++spin0) {
			// c(i1,orb1,spin0)
			SparseMatrixType O1 = model_.naturalOperator("c",site,spin0).data;
			// c(i2,orb2,1-spin0)
			SparseMatrixType O2 = model_.naturalOperator("c",site,1-spin0).data;
			for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
				// c(i2,orb2,spin1)
				SparseMatrixType O3 = model_.naturalOperator("c",site,spin1).data;
				// c(i3,orb1,1-spin1)
				SparseMatrixType O4 = model_.naturalOperator("c",site,1-spin1).data;
				SizeType val = spin0 + spin1 + 1;
				int signTerm = (val & 1) ? sign : 1;
				sum +=  signTerm*observe_.fourpoint()('N',
				                                      thini1,
				                                      O1,
				                                      'N',
				                                      thini2,
				                                      O2,
				                                      'C',
				                                      thinj1,
				                                      O3,
				                                      'C',
				                                      thinj2,
				                                      O4,
				                                      fermionicSign,
				                                      threadId);
			}
		}

		return sum;
	}

	void manyPoint(MatrixType* storage,
	               const BraketType& braket,
	               SizeType rows,
	               SizeType cols)
	{
		if (hasTimeEvolution_) {
			SizeType threadId = 0;
			printSites(threadId);
			std::cout<<"Time="<<observe_.time(threadId)<<"\n";
		}

		std::cout<<braket.toString()<<"\n";
		observe_.setBrakets(braket.bra(), braket.ket());

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

		if (braket.points() == 3)
			return observe_.threePoint(braket,rows,cols);


		if (braket.points() == 4)
			return observe_.fourPoint(braket,rows,cols);

		observe_.anyPoint(braket);
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
			const FieldType f1 = (-1.0);
			A += f1*matrixNdown_.data;
			OperatorType opA(A,1,std::pair<SizeType,SizeType>(0,0),1,su2Related1);
			PreOperatorSiteIndependentType preOperator(opA,"sz",threadId);
			measureOnePoint(preOperator);
		} else {
			PsimagLite::String s = "Unknown label: " + label + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void resizeStorage(VectorMatrixType& v,
	                   SizeType rows,
	                   SizeType cols,
	                   SizeType orbitals)
	{
		if (v.size() != 0) return;
		v.resize(static_cast<SizeType>(orbitals*(orbitals+1)/2));
		for (SizeType i = 0; i < v.size(); ++i)
			v[i].resize(rows,cols);
	}

	void measureOnePoint(const PreOperatorBaseType& preOperator)
	{
		const PsimagLite::String& modelName = model_.params().model;
		SizeType threadId = preOperator.threadId();
		VectorFieldType density;

		for (SizeType i0 = 0;i0<observe_.size();i0++) {
			if (!preOperator.isValid(i0+1)) continue;

			OperatorType opA = preOperator(i0+1);
			FieldType tmp1;

			preOperator.printMatrix(opA.data,preOperator.siteDependent(),i0);

			if (i0==0) {
				std::cout<<"site "<<preOperator.label()<<"(gs) ";
				if (hasTimeEvolution_)
					std::cout<<preOperator.label()<<"(timevector) time";
				//std::cout<<"\n";
			}
			// for g.s. use this one:
			observe_.setBrakets("gs","gs");
			observe_.setPointer(threadId,i0);

			onePointHookForZero(i0,opA,"gs",threadId, density);

			tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0,opA,ApplyOperatorType::BORDER_NO);
			density.push_back(tmp1);

			if (hasTimeEvolution_) { // for time vector use this one:
				observe_.setBrakets("time","time");

				onePointHookForZero(i0,opA,"time",threadId,density);

				FieldType tmp2 = observe_.template
				        onePoint<ApplyOperatorType>(i0,opA,ApplyOperatorType::BORDER_NO);

				std::cout<<" "<<tmp2<<" "<<observe_.time(threadId);
			}

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
				density.push_back(tmp1);

				if (hasTimeEvolution_) {// for time vector use this one:
					observe_.setBrakets("time","time");
					FieldType tmp2 = observe_.template
					        onePoint<ApplyOperatorType>(i0,
					                                    opAcorner,
					                                    ApplyOperatorType::BORDER_YES);
					std::cout<<" "<<tmp2<<" "<<observe_.time(threadId);
				}

				//std::cout<<"\n";
			}
		}

		if (modelName=="HubbardOneBandExtendedSuper") {
			SizeType nsite=observe_.size()/2+1;
			SizeType orbitals = 2;
			for (SizeType i0 = 0;i0<nsite;i0++)
				std::cout << i0 << " " << density[i0*orbitals+0]
				          << " " << density[i0*orbitals+1] << std::endl;
		} else {
			for (SizeType i0 = 0;i0<observe_.size();i0++)
				std::cout << i0 << " " << density[i0] << std::endl;
		}
	}

	void onePointHookForZero(SizeType i0,
	                         const OperatorType& opA,
	                         const PsimagLite::String& gsOrTime,
	                         SizeType threadId,
	                         VectorFieldType& density)
	{
		if (hasTimeEvolution_) return;
		if (observe_.site(threadId)!=1 || observe_.isAtCorner(numberOfSites_,threadId))
			return;
		assert(observe_.site(threadId)==1);
		FieldType tmp1 = observe_.template onePointHookForZero<ApplyOperatorType>(i0,opA);
		density.push_back(tmp1);
		if (hasTimeEvolution_ && gsOrTime=="time") std::cout<<"\n";
		if (!hasTimeEvolution_) std::cout<<"\n";
	}

	SizeType logBase2(SizeType x) const
	{
		SizeType counter = 0;
		while (x > 0) {
			x >>= 1;
			counter++;
		}

		return (counter == 0) ? counter : counter - 1;
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
		std::cout<<"Sites=";
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

