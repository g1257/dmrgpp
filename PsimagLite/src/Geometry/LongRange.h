/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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
/** \ingroup PsimagLite */
/*@{*/

/*! \file LongRange.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef PSI_GEOM_LONG_RANGE_H
#define PSI_GEOM_LONG_RANGE_H
#include "GeometryBase.h"

namespace PsimagLite {

template <typename ComplexOrRealType, typename InputType>
class LongRange : public GeometryBase<ComplexOrRealType, InputType> {

	using MatrixType = Matrix<ComplexOrRealType>;
	using RealType   = typename Real<ComplexOrRealType>::Type;

public:

	LongRange()
	    : linSize_(0)
	    , dofs_(0)
	    , maxConnections_(0)
	{ }

	LongRange(SizeType linSize, std::string goptions, InputType& io)
	    : linSize_(linSize)
	    , maxConnections_(0)
	{
		bool                                   hasEntangler = false;
		typename Real<ComplexOrRealType>::Type entangler    = 0;

		try {
			io.readline(entangler, "GeometryEntangler=");
			hasEntangler = true;
		} catch (std::exception&) { }

		io.readline(dofs_, "DegreesOfFreedom=");

		if (hasEntangler) {
			const SizeType n = dofs_ * linSize;
			matrix_.resize(n, n);
			setEntangler(entangler);
		} else {
			io.read(matrix_, "Connectors");
		}

		if (goptions != "ConstantValues" and goptions != "compact" and goptions != "none") {
			throw RuntimeError(std::string("GeometryOptions must be either ")
			                   + "ConstantValues or compact or none, not " + goptions
			                   + "\n");
		}

		if (goptions == "compact") {
			reinterpretMatrix();
		} else {
			if (dofs_ != matrix_.rows() / linSize) {
				throw RuntimeError("Wrong Connectors matrix size\n");
			}
		}

		checkConnectors(matrix_, linSize_);

		try {
			io.readline(maxConnections_, "GeometryMaxConnections=");
		} catch (std::exception& e) {
			if (!hasEntangler) {
				std::cerr << "Please add GeometryMaxConnections=0 or "
				             "some other number\n";
				throw e;
			}
		}
	}

	virtual void set(MatrixType& m, SizeType orbitals) const
	{
		m = matrix_;
		if (orbitals != dofs_)
			throw RuntimeError("General geometry connectors: wrong size\n");
	}

	virtual SizeType maxConnections() const
	{
		return (maxConnections_ == 0) ? linSize_ * linSize_ * 0.25 : maxConnections_;
	}

	virtual SizeType dirs() const { return 1; }

	SizeType handle(SizeType i, SizeType j) const { return (i < j) ? i : j; }

	SizeType getVectorSize(SizeType dirId) const
	{
		assert(dirId == 0);
		throw RuntimeError("LongRange::getVectorSize(): unimplemented\n");
	}

	bool connected(SizeType i1, SizeType i2) const { return true; }

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType, SizeType) const { return 0; }

	bool fringe(SizeType, SizeType, SizeType) const { return true; }

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax, SizeType emin, SizeType siteNew2) const
	{
		assert(siteNew2 >= emin);
		SizeType tmp = siteNew2 - emin + smax + 1;
		assert(tmp < linSize_);
		return tmp;
	}

	String label() const { return "LongRange"; }

	SizeType findReflection(SizeType site) const { return linSize_ - site - 1; }

	SizeType length(SizeType i) const
	{
		assert(i == 0);
		return linSize_;
	}

	SizeType translate(SizeType site, SizeType dir, SizeType amount) const
	{
		assert(dir == 0);

		site += amount;
		while (site >= linSize_)
			site -= linSize_;
		return site;
	}

	template <class Archive> void write(Archive&, const unsigned int)
	{
		throw RuntimeError("LongRange::write(): unimplemented\n");
	}

private:

	void setEntangler(ComplexOrRealType value)
	{
		const SizeType n = matrix_.rows();
		for (SizeType i = 0; i < n; ++i)
			for (SizeType j = i + 1; j < n; ++j)
				matrix_(i, j) = value;
	}

	static void checkConnectors(const MatrixType& matrix, SizeType linSize)
	{
		PsimagLite::String str;
		if (matrix.rows() != matrix.cols())
			str = "LongRange: Connectors matrix isn't square\n";

		if (matrix.rows() % linSize != 0)
			str += "LongRange: Connectors matrix isn't divisible "
			       "by number of sites\n";

		if (hasDiagonal(matrix)) {
			str += "LongRange: Connectors matrix has non-zero "
			       "diagonal value(s)\n";
		}

		if (hasLowerTriangle(matrix)) {
			str += "LongRange: Connectors matrix has non-zero(es) "
			       "in lower triangle\n";
		}

		if (str.empty())
			return;
		std::cerr << "WARNING: " << str;
	}

	static bool hasDiagonal(const MatrixType& matrix)
	{
		SizeType n = matrix.rows();
		if (n != matrix.cols())
			return true;
		for (SizeType i = 0; i < n; ++i) {
			if (PsimagLite::norm(matrix(i, i)) != 0)
				return true;
		}

		return false;
	}

	static bool hasLowerTriangle(const MatrixType& matrix)
	{
		for (SizeType i = 0; i < matrix.rows(); ++i) {
			for (SizeType j = 0; j < i; ++j) {
				if (PsimagLite::norm(matrix(i, j)) != 0)
					return true;
			}
		}

		return false;
	}

	// If orbitals == 1 then order is
	// site0 site1 value
	//
	// If orbitals > 1 then order is
	// site0 orb0 site1 orb1 value
	void reinterpretMatrix()
	{
		MatrixType values = matrix_;
		assert(linSize_ > 0);
		assert(dofs_ > 0);
		SizeType n = linSize_ * dofs_;
		matrix_.clear();
		matrix_.resize(n, n);
		if (dofs_ > 1) {
			assert(values.cols() == 5);
			for (SizeType i = 0; i < values.rows(); ++i) {
				SizeType counter = 0;
				SizeType site0   = complexToInteger(values(i, counter++));
				SizeType orb0    = complexToInteger(values(i, counter++));
				SizeType site1   = complexToInteger(values(i, counter++));
				SizeType orb1    = complexToInteger(values(i, counter++));
				matrix_(orb0 + site0 * dofs_, orb1 + site1 * dofs_)
				    = values(i, counter++);
			}
		} else {
			assert(values.cols() == 3);
			for (SizeType i = 0; i < values.rows(); ++i) {
				SizeType counter      = 0;
				SizeType site0        = complexToInteger(values(i, counter++));
				SizeType site1        = complexToInteger(values(i, counter++));
				matrix_(site0, site1) = values(i, counter++);
			}
		}
	}

	static SizeType complexToInteger(const ComplexOrRealType& value)
	{
		if (std::imag(value) != 0) {
			throw RuntimeError(std::string("Expected an integer in Connectors matrix ")
			                   + +" with compact option, not a complex number\n");
		}

		RealType val = std::real(value);
		if (!isInt64(val)) {
			throw RuntimeError(
			    std::string("Expected an integer in Connectors matrix ")
			    + +" with compact option, not a floating point number\n");
		}

		SizeType valInt = static_cast<SizeType>(val);
		return valInt;
	}

	// Credit for the idea: https://sillycross.github.io/index.html
	static bool isInt64(double d)
	{
		if (-9223372036854775808.0 <= d && d < 9223372036854775808.0) {
			return d == static_cast<double>(static_cast<int64_t>(d));
		} else {
			return false;
		}
	}

	SizeType   linSize_;
	SizeType   dofs_;
	SizeType   maxConnections_;
	MatrixType matrix_;
}; // class LongRange
} // namespace PsimagLite

/*@}*/
#endif // LADDER_H
