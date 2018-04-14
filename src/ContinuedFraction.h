
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
#include <typeinfo>

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
	                  const ParametersType& params)
	    : progress_("ContinuedFraction"),
	      freqEnum_(FREQ_REAL),
	      ab_(ab),
	      Eg_(params.Eg),
	      weight_(params.weight),
	      isign_(params.isign)
	{
		diagonalize();
	}

	ContinuedFraction(FreqEnum freqEnum = FREQ_REAL) : progress_("ContinuedFraction"),
	    freqEnum_(freqEnum),ab_(),Eg_(0),weight_(0),isign_(1) { }

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

		io.readline(weight_,"#CFWeight=");
		io.readline(Eg_,"#CFEnergy=");
		io.readline(isign_,"#CFIsign=");
		io.read(eigs_,"#CFEigs");
		io.read(intensity_,"#CFIntensities");
		diagonalize();
	}

	template<typename SomeIoOutputType>
	void save(SomeIoOutputType&) const
	{
		String name(typeid(SomeIoOutputType).name());
		std::cerr<<"WARNING: cannot save ContinuedFraction";
		std::cerr<<"to output type "<<name<<"\n";
	}

	void save(IoSimple::Out& io) const
	{
		io.setPrecision(12);
		ab_.save(io);

		String f = (freqEnum_ == FREQ_MATSUBARA) ? "Matsubara" : "Real";
		io.write(f, "#FreqEnum=");

		io.write(weight_, "#CFWeight=");

		io.write(Eg_, "#CFEnergy=");

		io.write(isign_, "#CFIsign=");

		io.write(eigs_,"#CFEigs");
		io.write(intensity_,"#CFIntensities");
	}

	void set(const TridiagonalMatrixType& ab,
	         const RealType& Eg,
	         RealType weight,
	         int isign)
	{
		ab_ = ab;
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
			std::pair<RealType,ComplexType> p(PsimagLite::imag(z),res);
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
		eigs_.resize(T.rows());
		diag(T,eigs_,'V');
		intensity_.resize(T.rows());
		for (SizeType i=0;i<T.rows();i++) {
			intensity_[i]= T(0,i)*T(0,i);
		}
	}

	RealType matsubara(int ind,const PlotParamsType& params) const
	{
		return (2.0*ind + 1.0)*M_PI/params.beta;
	}

	ProgressIndicator progress_;
	FreqEnum freqEnum_;
	TridiagonalMatrixType ab_;
	RealType Eg_;
	RealType weight_;
	int isign_;
	typename Vector<RealType>::Type eigs_;
	typename Vector<RealType>::Type intensity_;
}; // class ContinuedFraction
} // namespace PsimagLite
/*@}*/
#endif  //CONTINUED_FRACTION_H

