
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
/** \ingroup PsimagLite */
/*@{*/

/*! \file ContinuedFraction.h
 *
 * A continued fraction as explained in, e.g.,
 * E. Dagotto, Rev. Mod. Phys., 66, 763, (2004).
 */

#ifndef CONTINUED_FRACTION_H
#define CONTINUED_FRACTION_H
#include <iostream>
#include "Complex.h"
#include "TypeToString.h"
#include "ProgressIndicator.h"
#include "Random48.h"
#include "PlotParams.h"
#include "ParametersForSolver.h"
#include "IoSimple.h"
#include "FreqEnum.h"

namespace PsimagLite {
template<typename TridiagonalMatrixType_>
class ContinuedFraction  {
public:

	typedef TridiagonalMatrixType_ TridiagonalMatrixType;
	typedef typename TridiagonalMatrixType::value_type MatrixElementType;
	typedef typename Real<MatrixElementType>::Type RealType;
	typedef typename std::complex<RealType> ComplexType;
	typedef typename TridiagonalMatrixType::value_type FieldType;
	typedef Matrix<FieldType> MatrixType;
	typedef Matrix<RealType> MatrixRealType;
	typedef typename Vector<std::pair<RealType,ComplexType> >::Type PlotDataType;
	typedef PlotParams<RealType> PlotParamsType;
	typedef ParametersForSolver<RealType> ParametersType;

	ContinuedFraction(const TridiagonalMatrixType& ab,
	                  const MatrixType& reortho,
	                  const ParametersType& params)
	    : progress_("ContinuedFraction"),
	      freqEnum_(FREQ_REAL),
	      ab_(ab),
	      reortho_(reortho),
	      Eg_(params.Eg),
	      weight_(params.weight),
	      isign_(params.isign)
	{
		diagonalize();
	}

	ContinuedFraction(FreqEnum freqEnum = FREQ_REAL) : progress_("ContinuedFraction"),
	    freqEnum_(freqEnum),ab_(),reortho_(),Eg_(0),weight_(0),isign_(1) { }

	ContinuedFraction(IoSimple::In& io)
	    : progress_("ContinuedFraction"), freqEnum_(FREQ_REAL),ab_(io)
	{
		String f;
		try {
			io.readline(f,"#FreqEnum=");
		} catch(std::exception& e) {
			std::cerr<<"ContinuedFraction: FreqEnum assumed REAL\n";
			f = "Real";
			io.rewind();
		}

		if (f == "Matsubara") freqEnum_ = FREQ_MATSUBARA;
		try {
			io.readMatrix(reortho_,"#ReorthogonalizationMatrix");
		} catch(std::exception& e) {
			assert(reortho_.n_row() == 0);
			std::cerr<<"ContinuedFraction: #ReorthogonalizationMatrix (nrow=";
			std::cerr<<reortho_.n_row()<<") DISABLED\n";
			io.rewind();
		}

		io.readline(weight_,"#CFWeight=");
		io.readline(Eg_,"#CFEnergy=");
		io.readline(isign_,"#CFIsign=");
		io.read(eigs_,"#CFEigs");
		io.read(intensity_,"#CFIntensities");
		diagonalize();
	}

	template<typename IoOutputType>
	void save(IoOutputType& io) const
	{
		io.setPrecision(12);
		ab_.save(io);

		String f = (freqEnum_ == FREQ_MATSUBARA) ? "Matsubara" : "Real";
		io.print("#FreqEnum=",f);
		io.printMatrix(reortho_,"#ReorthogonalizationMatrix");

		io.print("#CFWeight=",weight_);

		io.print("#CFEnergy=",Eg_);

		io.print("#CFIsign=" ,isign_);

		io.printVector(eigs_,"#CFEigs");
		io.printVector(intensity_,"#CFIntensities");
	}

	void set(const TridiagonalMatrixType& ab,
	         const MatrixRealType& reortho,
	         const RealType& Eg,
	         RealType weight,
	         int isign)
	{
		ab_ = ab;
		reortho_ = reortho;
		Eg_ = Eg;
		weight_ = weight;
		isign_ = isign;

		diagonalize();
	}

	void plot(PlotDataType& result,const PlotParamsType& params) const
	{
		if (freqEnum_ == FREQ_MATSUBARA || params.numberOfMatsubaras > 0) {
			plotMatsubara(result,params);
			return;
		}

		if (freqEnum_ == FREQ_REAL || params.numberOfMatsubaras == 0) {
			plotReal(result,params);
		}
	}

	void plotReal(PlotDataType& result,const PlotParamsType& params) const
	{
		SizeType counter = 0;
		SizeType n = SizeType((params.omega2 - params.omega1)/params.deltaOmega);
		if (result.size()==0) result.resize(n);
		for (RealType omega=params.omega1;omega<params.omega2;omega+=params.deltaOmega) {
			ComplexType z(omega,params.delta);
			ComplexType res = iOfOmega(z,Eg_,isign_);
			std::pair<RealType,ComplexType> p(omega,res);
			result[counter++] = p;
			if (counter>=result.size()) break;
		}
	}

	void plotMatsubara(PlotDataType& result,const PlotParamsType& params) const
	{
		SizeType counter = 0;
		SizeType n = params.numberOfMatsubaras;
		if (result.size()==0) result.resize(n);
		for (SizeType omegaIndex = 0; omegaIndex < params.numberOfMatsubaras; ++omegaIndex) {
			ComplexType z(params.delta, matsubara(omegaIndex,params));
			ComplexType res = iOfOmega(z,Eg_,isign_);
			std::pair<RealType,ComplexType> p(std::imag(z),res);
			result[counter++] = p;
			if (counter>=result.size()) break;
		}
	}
	//! Cases:
	//! (1) < phi0|A (z+(E0-e_k))^{-1}|A^\dagger|phi0> and
	//! (2) < phi0|A^\dagger (z-(E0-e_k))^{-1}|A|phi0>
	//! (There are actually 4 cases for the off-diagonal gf because
	//! A has two cases:
	//! (1) A = c_i + c_j and
	//! (2) A = c_i - c_j
	ComplexType iOfOmega(const ComplexType& z,RealType offset,int isign) const

	{
		if (weight_==0) return ComplexType(0,0);

		ComplexType sum = 0;
		for (SizeType l=0;l<intensity_.size();l++)
			sum +=intensity_[l]/(z-isign*(offset - eigs_[l]));

		return sum*weight_;
	}

	SizeType size() const { return ab_.size(); }

	FreqEnum freqType() const  { return freqEnum_; }

private:

	void diagonalize()
	{
		if (weight_==0) return;
		MatrixType T;
		ab_.buildDenseMatrix(T);

		if (reortho_.n_row()>0) {
			MatrixType tmp = T * reortho_;
			T = multiplyTransposeConjugate(reortho_,tmp);
		}

		eigs_.resize(T.n_row());
		diag(T,eigs_,'V');
		intensity_.resize(T.n_row());
		for (SizeType i=0;i<T.n_row();i++) {
			intensity_[i]= T(0,i)*T(0,i);
		}
	}

	RealType matsubara(int ind,const PlotParamsType& params) const
	{
		int halfNs = static_cast<int>(params.omega1);
		RealType factor = 2.0*M_PI/params.beta;
		int ind2 = ind - halfNs;
		if (ind2 >= 0) return factor*(ind2 + 1);
		return factor*ind2;
	}

	ProgressIndicator progress_;
	FreqEnum freqEnum_;
	TridiagonalMatrixType ab_;
	MatrixRealType reortho_;
	RealType Eg_;
	RealType weight_;
	int isign_;
	typename Vector<RealType>::Type eigs_;
	typename Vector<RealType>::Type intensity_;
}; // class ContinuedFraction
} // namespace PsimagLite
/*@}*/
#endif  //CONTINUED_FRACTION_H

