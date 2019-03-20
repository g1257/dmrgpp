/*
Copyright (c) 2009-2014-2019, UT-Battelle, LLC
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
#include "Concurrency.h"
#include "Vector.h"
#include "ProgramGlobals.h"
#include "ApplyOperatorLocal.h"

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
	typedef typename ObserverType::BraketType BraketType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename ObserverType::PointerForSerializerType PointerForSerializerType;

	template<typename IoInputter>
	ObservableLibrary(IoInputter& io,
	                  SizeType numberOfSites,
	                  const ModelType& model,
	                  SizeType start,
	                  SizeType nf,
	                  SizeType trail)
	    : numberOfSites_(numberOfSites),
	      model_(model),
	      observe_(io, start, nf, trail, model.params())
	{}

	bool endOfData() const { return observe_.helper().endOfData(); }

	const ModelType& model() const { return model_; }

	void interpret(const PsimagLite::String& list, SizeType rows, SizeType cols)
	{
		typename BraketType::VectorStringType vecStr;
		PsimagLite::split(vecStr, list, ",");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			BraketType braket(model_, vecStr[i]);

			if (braket.points() == 1) {
				measureOnePoint(braket.bra(),
				                braket.op(0),
				                braket.opName(0),
				                braket.ket());
				continue;
			}

			manyPoint(0,braket,rows,cols);
		}
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
					observe_.twoPoint(out,
					                  n1,
					                  n2,
					                  ProgramGlobals::FermionOrBosonEnum::BOSON,
					                  "gs",
					                  "gs");
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
			ppTwopoint(results, names, rows, cols, "gs", "gs");
			ppFourpoint(results, names, rows, cols);

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
			observe_.multiCorrelations(result, myMatrixSparse, rows, cols, "gs", "gs");
			for (SizeType i=0;i<result.size();i++)
				std::cout<<i<<" "<<result[i]<<"\n";

		} else if (label == "ddOrbitals") {
			if (model_.params().model!="HubbardOneBandExtendedSuper") {
				throw PsimagLite::RuntimeError("pp: not for this model\n");
			}
			typename PsimagLite::Vector<MatrixType*>::Type results;
			typename PsimagLite::Vector<PsimagLite::String>::Type names;
			ddOrbitalsTwopoint(results, names, rows, cols, "gs", "gs");
			ddOrbitalsFourpoint(results,names,rows,cols);
		} else {
			err("Unknown label: " + label + "\n");
		}
	}

private:

	void measureOnePoint(const PsimagLite::String& bra,
	                     const OperatorType& opA,
	                     PsimagLite::String label,
	                     const PsimagLite::String& ket)
	{
		for (SizeType i0 = 0; i0 < observe_.helper().size(); ++i0) {

			if (i0==0) {
				std::cout<<"Using Matrix A:\n";
				std::cout<<opA.data.toDense();
				std::cout<<"site <"<<bra<<"|"<<label;
				std::cout<<"|"<<ket<<"> time\n";
			}

			// observe_.setBrakets(bra,ket);
			PointerForSerializerType ptr(i0);
			cornerLeftOrRight(1, ptr, bra, opA, ket);

			FieldType tmp1 = observe_.template
			        onePoint<ApplyOperatorType>(i0, opA, ApplyOperatorType::BORDER_NO, bra, ket);
			std::cout<<observe_.helper().site(ptr)<<" "<<tmp1;
			std::cout<<" "<<observe_.helper().time(ptr)<<"\n";

			cornerLeftOrRight(numberOfSites_ - 2, ptr, bra, opA, ket);
		}
	}

	void cornerLeftOrRight(SizeType site,
	                       const PointerForSerializerType& ptr,
	                       const PsimagLite::String& bra,
	                       const OperatorType& opA,
	                       const PsimagLite::String& ket)
	{
		if (observe_.helper().site(ptr) != site) return;

		bool atCorner = observe_.isAtCorner(numberOfSites_, ptr);

		if (site == 1 && !atCorner) {
			const typename ApplyOperatorType::BorderEnum borderNo =
			        ApplyOperatorType::BorderEnum::BORDER_NO;
			FieldType tmp1 = observe_.template onePointHookForZero<ApplyOperatorType>(ptr.get(),
			                                                                          opA,
			                                                                          borderNo,
			                                                                          bra,
			                                                                          ket);
			std::cout<<"0 "<<tmp1<<" "<<observe_.helper().time(ptr)<<"\n";
			return;
		}

		if (!atCorner)
			return;

		// also calculate next or prev. site:
		SizeType x = (observe_.helper().site(ptr) == 1) ? 0 : numberOfSites_ - 1;

		// operator might be site dependent
		OperatorType opAcorner = opA;

		// do the corner case
		// observe_.setBrakets(bra,ket);
		FieldType tmp1 = observe_.template
		        onePoint<ApplyOperatorType>(ptr.get(),
		                                    opAcorner,
		                                    ApplyOperatorType::BORDER_YES,
		                                    bra,
		                                    ket);
		std::cout<<x<<" "<<tmp1;
		std::cout<<" "<<observe_.helper().time(ptr)<<"\n";
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

	void ppTwopoint(typename PsimagLite::Vector<MatrixType*>::Type& result,
	                typename PsimagLite::Vector<PsimagLite::String>::Type& names,
	                SizeType rows,
	                SizeType cols,
	                PsimagLite::String bra,
	                PsimagLite::String ket)
	{
		// Two-point Pair
		MatrixType m1(rows,cols);
		MatrixType m2(rows,cols);
		ppTwo(m1, m2, 0, bra, ket);

		if (model_.params().model=="HubbardOneBandExtendedSuper") {
			rows = rows/2;   // Actually: divided by # of orbitals
			cols = cols/2;
		}

		m1.clear(); m1.resize(rows,cols);
		m2.clear(); m2.resize(rows,cols);
		ppTwo(m1, m2, 1, bra, ket);
		std::cout << "PairPair Correlations S^{lu}_{on}" << std::endl;
		std::cout << m1;
		std::cout << "PairPair Correlations T^{lu}_{on}" << std::endl;
		std::cout << m2;

		m1.clear(); m1.resize(rows,cols);
		m2.clear(); m2.resize(rows,cols);
		names.push_back("T_{on}");
		ppTwo(m1, m2, 3, bra, ket);
	}

	void ppTwo(MatrixType& m,
	           MatrixType& m2,
	           SizeType flag,
	           PsimagLite::String bra,
	           PsimagLite::String ket)
	{
		SizeType site = 1;
		SizeType orbitals = logBase2(model_.hilbertSize(site));
		assert(!(orbitals & 1));
		orbitals /= 2;

		if (flag==0) {
			SizeType spin0 = 0; // up
			SizeType spin1 = 1; // down

			// c_dn,0
			SparseMatrixType O1 = model_.naturalOperator("c",site,spin1).data;
			// c_up,0
			SparseMatrixType O2 = model_.naturalOperator("c",site,spin0).data;

			SparseMatrixType A,B;
			multiply(B,O1,O2); // c_dn,0 . c_up,0.
			transposeConjugate(A,B);
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
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
			err("Unknown flag: " + ttos(flag));
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
					pairs.push_back(PairSizeType(thini1 + thini2*rows*orbitals +
					                             rows*orbitals*rows*orbitals*spin0,
					                             thinj1 + thinj2*rows*orbitals +
					                             rows*orbitals*rows*orbitals*spin1));
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
	                        SizeType cols,
	                        PsimagLite::String bra,
			                PsimagLite::String ket)
	{
		// Two-point Pair
		MatrixType* m1 = new MatrixType(rows,cols);
		names.push_back("S^{l}_{on}");
		std::cout << "PairPair Correlations S^{l}_{on}" << std::endl;
		ddOrbitalsTwo(*m1, 0, bra, ket);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("S^{u}_{on}");
		std::cout << "PairPair Correlations S^{u}_{on}" << std::endl;
		ddOrbitalsTwo(*m1, 1, bra, ket);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("S^{lu}_{on}");
		std::cout << "PairPair Correlations S^{lu}_{on}" << std::endl;
		ddOrbitalsTwo(*m1, 2, bra, ket);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{lu}_{on}" << std::endl;
		ddOrbitalsTwo(*m1, 3, bra, ket);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{up-up}_{on}" << std::endl;
		ddOrbitalsTwo(*m1, 4, bra, ket);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{dn-dn}_{on}" << std::endl;
		ddOrbitalsTwo(*m1, 5, bra, ket);
		result.push_back(m1);

		m1 = new MatrixType(rows,cols);
		names.push_back("T^{lu}_{on}");
		std::cout << "PairPair Correlations T^{up*up+dn*dn}_{on}" << std::endl;
		ddOrbitalsTwo(*m1, 6, bra, ket);
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

	void ddOrbitalsTwo(MatrixType& m,
	                   const SizeType flag,
	                   PsimagLite::String bra,
	                   PsimagLite::String ket)
	{
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
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);

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
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
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
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
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
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
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
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
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
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
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
			observe_.twoPoint(m, A, B, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
			std::cout << m;
		} else {
			err("Unknown flag: " + ttos(flag));
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
				sum +=  signTerm*observe_.fourpoint()(i1,i2,j1,j2,braket);
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

	void manyPoint(MatrixType* storage,
	               const BraketType& braket,
	               SizeType rows,
	               SizeType cols)
	{
		if (hasTimeEvolution_) {
			PointerForSerializerType ptr(0);
			printSites();
			std::cout<<"Time="<<observe_.helper().time(ptr)<<"\n";
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

		if (braket.points() == 3)
			return observe_.threePoint(braket,rows,cols);


		if (braket.points() == 4)
			return observe_.fourPoint(braket,rows,cols);

		observe_.anyPoint(braket);
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

	void printSites() const
	{
		std::cout<<"Sites=";
		PointerForSerializerType ptr(0);
		if (observe_.helper().site(ptr) == 1) std::cout<<"0 ";
		if (observe_.helper().site(ptr) == numberOfSites_ - 2)
			std::cout<<(numberOfSites_ - 1)<<" ";
		for (SizeType i = 0; i < observe_.helper().size(); ++i) {
			ptr.setPointer(i);
			SizeType x = observe_.helper().site(ptr);
			std::cout<<x<<" ";
		}

		if (observe_.helper().site(ptr) == 1) std::cout<<"0";
		if (observe_.helper().site(ptr) == numberOfSites_ - 2)
			std::cout<<(numberOfSites_ - 1);

		std::cout<<"\n";
	}

	SizeType numberOfSites_;
	bool hasTimeEvolution_;
	const ModelType& model_; // not the owner
	ObserverType observe_;
	VectorMatrixType szsz_,sPlusSminus_,sMinusSplus_;

}; // class ObservableLibrary

} // namespace Dmrg

/*@}*/
#endif // OBSERVABLE_LIBRARY_H

