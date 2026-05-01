
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
#include "Complex.h"
#include "FreqEnum.h"
#include "Io/IoSimple.h"
#include "Matsubaras.h"
#include "ParametersForSolver.h"
#include "ProgressIndicator.h"
#include "Random48.h"
#include "RealFrequencyRange.hh"
#include "TridiagonalMatrix.h"
#include "TypeToString.h"
#include <iostream>
#include <typeinfo>

namespace PsimagLite {
template <typename RealType> class ContinuedFraction {
public:

	using TridiagonalMatrixType = TridiagonalMatrix<RealType>;
	using ComplexType           = typename std::complex<RealType>;
	using MatrixType            = Matrix<RealType>;
	using MatrixRealType        = Matrix<RealType>;
	using PlotDataType          = std::vector<std::pair<RealType, ComplexType>>;
	using ParametersType        = ParametersForSolver<RealType>;

	ContinuedFraction(const TridiagonalMatrixType& ab, const ParametersType& params)
	    : progress_("ContinuedFraction")
	    , freqEnum_(FreqEnum::REAL)
	    , ab_(ab)
	    , Eg_(params.Eg)
	    , weight_(params.weight)
	    , isign_(params.isign)
	{
		diagonalize();
	}

	ContinuedFraction(FreqEnum freqEnum = FreqEnum::REAL)
	    : progress_("ContinuedFraction")
	    , freqEnum_(freqEnum)
	    , ab_()
	    , Eg_(0)
	    , weight_(0)
	    , isign_(1)
	{ }

	ContinuedFraction(IoSimple::In& io)
	    : progress_("ContinuedFraction")
	    , freqEnum_(FreqEnum::REAL)
	    , ab_(io)
	{
		String f;
		try {
			io.readline(f, "#FreqEnum=");
		} catch (std::exception& e) {
			std::cerr << "ContinuedFraction: FreqEnum assumed REAL\n";
			f = "Real";
			io.rewind();
		}

		if (f == "Matsubara")
			freqEnum_ = FreqEnum::MATSUBARA;

		io.readline(weight_, "#CFWeight=");
		io.readline(Eg_, "#CFEnergy=");
		io.readline(isign_, "#CFIsign=");
		io.read(eigs_, "#CFEigs");
		io.read(intensity_, "#CFIntensities");
		diagonalize();
	}

	template <typename SomeIoOutputType> void write(SomeIoOutputType&, String) const
	{
		String name(typeid(SomeIoOutputType).name());
		std::cerr << "WARNING: cannot save ContinuedFraction";
		std::cerr << "to output type " << name << "\n";
	}

	void write(IoSimple::Out& io, String) const
	{
		io.setPrecision(12);
		ab_.write(io);

		String f = (freqEnum_ == FreqEnum::MATSUBARA) ? "Matsubara" : "Real";
		io.write(" ", "#FreqEnum=" + f);

		io.write(weight_, "#CFWeight");

		io.write(Eg_, "#CFEnergy");

		io.write(isign_, "#CFIsign");

		io.write(eigs_, "#CFEigs");
		io.write(intensity_, "#CFIntensities");
	}

	void set(const TridiagonalMatrixType& ab, const RealType& Eg, RealType weight, int isign)
	{
		ab_     = ab;
		Eg_     = Eg;
		weight_ = weight;
		isign_  = isign;

		diagonalize();
	}

	void plot(PlotDataType& result, const RealFrequencyRange<RealType>& params) const
	{
		SizeType total = params.total();
		if (total == 0) {
			return;
		}

		if (freqEnum_ != FreqEnum::REAL) {
			throw std::runtime_error("ContinuedFraction::plot() matsubaras/realfreq. "
			                         "mismatch: Real expected\n");
		}

		result.resize(total);
		for (SizeType i = 0; i < total; ++i) {
			RealType                         omega = params.omega(i);
			ComplexType                      z(omega, params.delta());
			ComplexType                      res = iOfOmega(z, Eg_, isign_);
			std::pair<RealType, ComplexType> p(omega, res);
			result[i] = p;
		}
	}

	void plot(PlotDataType& result, const Matsubaras<RealType>& matsubaras) const
	{
		SizeType n = matsubaras.total();
		if (n == 0) {
			return;
		}

		if (freqEnum_ != FreqEnum::MATSUBARA) {
			throw std::runtime_error("ContinuedFraction::plot() matsubaras/realfreq. "
			                         "mismatch: Matsubaras expected\n");
		}

		result.resize(n);
		for (SizeType omegaIndex = 0; omegaIndex < n; ++omegaIndex) {
			ComplexType z(matsubaras.delta(), matsubaras.omega(omegaIndex));
			ComplexType res = iOfOmega(z, Eg_, isign_);
			std::pair<RealType, ComplexType> p(PsimagLite::imag(z), res);
			result[omegaIndex] = p;
		}
	}

	//! Cases:
	//! (1) < phi0|A (z+(E0-e_k))^{-1}|A^\dagger|phi0> and
	//! (2) < phi0|A^\dagger (z-(E0-e_k))^{-1}|A|phi0>
	//! (There are actually 4 cases for the off-diagonal gf because
	//! A has two cases:
	//! (1) A = c_i + c_j and
	//! (2) A = c_i - c_j
	ComplexType iOfOmega(const ComplexType& z, RealType offset, int isign) const

	{
		if (PsimagLite::real(weight_) == 0 && PsimagLite::imag(weight_) == 0)
			return ComplexType(0, 0);

		ComplexType sum = 0;
		for (SizeType l = 0; l < intensity_.size(); l++)
			sum += intensity_[l] / (z - isign * (offset - eigs_[l]));

		return sum * weight_;
	}

	SizeType size() const { return ab_.size(); }

	FreqEnum freqType() const { return freqEnum_; }

private:

	void diagonalize()
	{
		if (PsimagLite::real(weight_) == 0 && PsimagLite::imag(weight_) == 0)
			return;

		MatrixType T;
		ab_.buildDenseMatrix(T);
		eigs_.resize(T.rows());
		diag(T, eigs_, 'V');
		intensity_.resize(T.rows());
		for (SizeType i = 0; i < T.rows(); i++) {
			intensity_[i] = T(0, i) * T(0, i);
		}
	}

	ProgressIndicator               progress_;
	FreqEnum                        freqEnum_;
	TridiagonalMatrixType           ab_;
	RealType                        Eg_;
	std::complex<RealType>          weight_;
	int                             isign_;
	typename Vector<RealType>::Type eigs_;
	typename Vector<RealType>::Type intensity_;
}; // class ContinuedFraction
} // namespace PsimagLite
/*@}*/
#endif // CONTINUED_FRACTION_H
