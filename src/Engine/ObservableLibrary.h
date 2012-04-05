// BEGIN LICENSE BLOCK
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file ObservableLibrary.h
 *
 *  A library of most used observables
 */
#ifndef OBSERVABLE_LIBRARY_H
#define OBSERVABLE_LIBRARY_H

#include "Matrix.h" // in PsimagLite

namespace Dmrg {
	
	template<typename ObserverType,typename TargettingType>
	class ObservableLibrary {
	public:
		typedef typename TargettingType::ModelType ModelType;
		typedef typename ModelType::OperatorType OperatorType;
		typedef typename OperatorType::Su2RelatedType Su2RelatedType;
		typedef typename TargettingType::ApplyOperatorType ApplyOperatorType;
		typedef typename ModelType::ConcurrencyType ConcurrencyType;
		typedef typename ModelType::RealType RealType;

		typedef typename OperatorType::SparseMatrixType SparseMatrixType;
		typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
		typedef typename VectorWithOffsetType::value_type FieldType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;

		template<typename IoInputter>
		ObservableLibrary(
				IoInputter& io,
				size_t numberOfSites,
				bool hasTimeEvolution,
				const ModelType& model,
				ConcurrencyType& concurrency,
				bool verbose)
		: numberOfSites_(numberOfSites),
		  hasTimeEvolution_(hasTimeEvolution),
		  model_(model),
		  concurrency_(concurrency),
		  observe_(io,numberOfSites-2,hasTimeEvolution,model,
				  concurrency,verbose),
				 szsz_(0,0)
		{
			if (hasTimeEvolution) {
				size_t site = 0;
				// FIXME: No support for site varying operators
				fullMatrixToCrsMatrix(matrixNup_,model_.naturalOperator("nup",site,0));
				fullMatrixToCrsMatrix(matrixNdown_,model_.naturalOperator("ndown",site,0));
			}
		}

		void measure(const std::string& label,size_t rows,size_t cols)
		{
			// Note that I can't print sites when there no time evolution
			// since the DmrgSerializer doens't have sites yet
			// as opposed to the TimeSerializer
			if (hasTimeEvolution_) printSites();
			size_t site = 0; // FIXME: No support for site varying operators
			if (label=="cc") {
				MatrixType opC = model_.naturalOperator("c",site,0); // c_{0,0} spin up
				MatrixType opCtranspose = transposeConjugate(opC);
				measureOne("OperatorC",opC,opCtranspose,-1,rows,cols);
			} else if (label=="nn") {
				MatrixType opN = model_.naturalOperator("n",site,0);
				measureOne("OperatorN",opN,opN,1,rows,cols);
			} else if (label=="szsz") {
				MatrixType Sz = model_.naturalOperator("z",site,0);
				szsz_ = observe_.correlations(Sz,Sz,1,rows,cols);
				if (concurrency_.root()) {
					std::cout<<"OperatorSz:\n";
					std::cout<<szsz_;
				}
			} else if (label=="s+s-") {
				// Si^+ Sj^-
				const MatrixType& sPlus = model_.naturalOperator("+",site,0);
				MatrixType sPlusT = transposeConjugate(sPlus);
				sPlusSminus_ = observe_.correlations(sPlus,sPlusT,1,rows,cols);
				if (concurrency_.root()) {
						std::cout<<"OperatorSplus:\n";
						std::cout<<sPlusSminus_;
					}
			} else if (label=="s-s+") {
				// Si^- Sj^+
				const MatrixType& sMinus = model_.naturalOperator("-",site,0);
				MatrixType sMinusT = transposeConjugate(sMinus);
				sMinusSplus_= observe_.correlations(sMinus,sMinusT,1,rows,cols);
				if (concurrency_.root()) {
					std::cout<<"OperatorSminus:\n";
					std::cout<<sMinusSplus_;
				}
			} else if (label=="ss") {
				if (szsz_.n_row()==0) measure("szsz",rows,cols);
				if (sPlusSminus_.n_row()==0)  measure("s+s-",rows,cols);
				if (sMinusSplus_.n_row()==0)  measure("s-s+",rows,cols);

				MatrixType spinTotal(szsz_.n_row(),szsz_.n_col());

				for (size_t i=0;i<spinTotal.n_row();i++)
					for (size_t j=0;j<spinTotal.n_col();j++)
						spinTotal(i,j) = 0.5*(sPlusSminus_(i,j) +
								sMinusSplus_(i,j)) + szsz_(i,j);

				if (concurrency_.root()) {
					std::cout<<"SpinTotal:\n";
					std::cout<<spinTotal;
				}
			} else if (label=="dd") {

				const MatrixType& oDelta = model_.naturalOperator("d",site,0);
				MatrixType oDeltaT;
				transposeConjugate(oDeltaT,oDelta);
				measureOne("TWO-POINT DELTA-DELTA^DAGGER",oDelta,oDeltaT,1,
						rows,cols);
			} else if (label=="dd4") {
				for (size_t g=0;g<16;g++) {
					std::vector<FieldType> fpd;
					std::vector<size_t> gammas(4,0);
					gammas[0] = (g & 1);
					gammas[1] = (g & 2)>>1;
					gammas[2] = (g & 4) >> 2;
					gammas[3] = (g & 8) >> 3;
					observe_.fourPointDeltas(fpd,numberOfSites_,
							gammas,model_);
					for (size_t step=0;step<fpd.size();step++) {
						// step --> (i,j) FIXME
						std::cout<<step<<" "<<g<<" "<<fpd[step]<<"\n";
					}
				}
			} else {
				std::string s = "Unknown label: " + label + "\n";
				throw std::runtime_error(s.c_str());
			}
		}

		void measureTime(const std::string& label)
		{
			SparseMatrixType A;
			Su2RelatedType su2Related1;

			if (label=="superDensity") {
				A.makeDiagonal(model_.hilbertSize(observe_.site()),1.0);
				std::pair<size_t,size_t> zeroZero(0,0);
				OperatorType opIdentity(A,1,zeroZero,1,su2Related1);
				observe_.setBrackets("time","time");
				FieldType superDensity = observe_.template
						onePoint<ApplyOperatorType>(0,opIdentity);
				std::cout<<"SuperDensity(Weight of the timeVector)="<<
						superDensity<<"\n";
			} else if (label=="nupNdown") {
				multiply(A,matrixNup_,matrixNdown_);
				measureOnePoint(A,"nupNdown");
			} else if (label=="nup+ndown") {
				A = matrixNup_;
				A += matrixNdown_;
				measureOnePoint(A,"nup+ndown");
			} else if (label=="sz") {
				A = matrixNup_;
				A += (-1)*matrixNdown_;
				measureOnePoint(A,"sz");
			} else {
				std::string s = "Unknown label: " + label + "\n";
					throw std::runtime_error(s.c_str());
			}
		}

		void setBrackets(const std::string& left,const std::string& right)
		{
			observe_.setBrackets(left,right);
		}

		bool endOfData() const { return observe_.endOfData(); }

		void measureTheOnePoints(size_t numberOfDofs)
		{
			size_t site = 0; // FIXME : account for Hilbert spaces changing with site
			for (size_t dof=0;dof<numberOfDofs;dof++) {
				for (size_t dof2=dof;dof2<numberOfDofs;dof2++) {
					MatrixType opCup = model_.naturalOperator("c",site,dof);
					MatrixType opCdown = model_.naturalOperator("c",site,dof2);
					MatrixType opCupTranspose;
					transposeConjugate(opCupTranspose,opCup);
					MatrixType Afull = opCupTranspose * opCdown;
					SparseMatrixType A(Afull);
					std::string str("c^\\dagger(dof=");
					str += ttos(dof) + ") c(dof=" + ttos(dof2) + ")";
					measureOnePoint(A,str);
				}
			}

		}

	private:

		void measureOne(const std::string& label,
			const PsimagLite::Matrix<FieldType>& op1,
			const PsimagLite::Matrix<FieldType>& op2,
			int fermionSign,
			size_t rows,
			size_t cols)
		{
			const MatrixType& v =
				observe_.correlations(op1,op2,fermionSign,rows,cols);;
			if (concurrency_.root()) {
				if (hasTimeEvolution_)
					std::cout<<"#Time="<<observe_.time()<<"\n";
				std::cout<<label<<":\n";
				std::cout<<v;
			}
		}

		void measureOnePoint(const SparseMatrixType& A,const std::string& label)
		{
			Su2RelatedType su2Related1;
			printMarker();
			std::cout<<"#Using Matrix A:\n";
			for (size_t i=0;i<A.row();i++) {
				for (size_t j=0;j<A.col();j++)
					std::cout<<"#A("<<i<<","<<j<<")="<<A(i,j)<<" ";
				std::cout<<"\n";
			}
			OperatorType opA(A,1,std::pair<size_t,size_t>(0,0),1,su2Related1);
			std::cout<<"site "<<label<<"(gs) ";
			if (hasTimeEvolution_) std::cout<<label<<"(timevector) time";
			std::cout<<"\n";
			for (size_t i0 = 0;i0<observe_.size();i0++) {
				// for g.s. use this one:
				observe_.setBrackets("gs","gs");
				FieldType tmp1 = observe_.template
						onePoint<ApplyOperatorType>(i0,opA);
				std::cout<<observe_.site()<<" "<<tmp1;

				if (hasTimeEvolution_) { // for time vector use this one:
					observe_.setBrackets("time","time");
					FieldType tmp2 = observe_.template
							 onePoint<ApplyOperatorType>(i0,opA);

					std::cout<<" "<<tmp2<<" "<<observe_.time();
				}
				std::cout<<"\n";
				// also calculate next or prev. site:
				if (observe_.isAtCorner(numberOfSites_)) {
					size_t x = (observe_.site()==1) ? 0 : numberOfSites_-1;
					// do the corner case
					// for g.s. use this one:
					observe_.setBrackets("gs","gs");
					bool doCorner = true;
					FieldType tmp1 = observe_.template
							onePoint<ApplyOperatorType>(i0,opA,doCorner);
					std::cout<<x<<" "<<tmp1;

					if (hasTimeEvolution_) {// for time vector use this one:
						observe_.setBrackets("time","time");
						FieldType tmp2 = observe_.template
								 onePoint<ApplyOperatorType>(i0,opA,doCorner);
						std::cout<<" "<<tmp2<<" "<<observe_.time();
					}
					std::cout<<"\n";
				}
			}
		}

		void printSites()
		{
			printMarker();
			std::cout<<"#Sites=";
			observe_.setPointer(0);
			if (observe_.site()==1) std::cout<<"0 ";
			if (observe_.site()==numberOfSites_-2)
				std::cout<<(numberOfSites_-1)<<" ";
			for (size_t i=0;i<observe_.size();i++) {
				observe_.setPointer(i);
				size_t x = observe_.site();
				std::cout<<x<<" ";
			}
			if (observe_.site()==1) std::cout<<"0";
			if (observe_.site() ==numberOfSites_-2) {
				std::cout<<(numberOfSites_-1);
			}

			std::cout<<"\n";
		}

		void printMarker()
		{
			if (!hasTimeEvolution_) return;
			size_t marker = observe_.marker();
			std::string s = "INVALID MARKER";
			switch (marker) {
			case 0:
				s=" NOT ALL OPERATORS APPLIED YET";
				break;
			case 1:
				s=" ALL OPERATORS HAVE BEEN APPLIED";
			}

			std::cout<<s<<"\n";
		}

		size_t numberOfSites_;
		bool hasTimeEvolution_;
		const ModelType& model_; // not the owner
		ConcurrencyType& concurrency_; // not the owner
		ObserverType observe_;
		SparseMatrixType matrixNup_,matrixNdown_;
		MatrixType szsz_,sPlusSminus_,sMinusSplus_;

	}; // class ObservableLibrary


} // namespace Dmrg 

/*@}*/
#endif // OBSERVABLE_LIBRARY_H
