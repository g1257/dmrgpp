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

template<typename ComplexOrRealType, typename InputType>
class LongRange : public GeometryBase<ComplexOrRealType, InputType> {

	typedef Matrix<ComplexOrRealType> MatrixType;

public:

	LongRange() : linSize_(0), orbitals_(0), maxConnections_(0) {}

	LongRange(SizeType linSize,InputType& io)
	    : linSize_(linSize), maxConnections_(0)
	{
		io.read(matrix_, "Connectors");
		assert(matrix_.rows()%linSize == 0);
		orbitals_ = static_cast<SizeType>(matrix_.rows()/linSize);
		try {
			io.readline(maxConnections_,"GeometryMaxConnections=");
		} catch (std::exception&) {}
	}

	virtual void set(MatrixType& m, SizeType& orbitals) const
	{
		m = matrix_;
		orbitals = orbitals_;
	}

	virtual SizeType maxConnections() const
	{
		return (maxConnections_ == 0) ? linSize_*linSize_*0.25 : maxConnections_;
	}

	virtual SizeType dirs() const { return 1; }

	SizeType handle(SizeType i,SizeType j) const
	{
		return (i<j) ? i : j;
	}

	SizeType getVectorSize(SizeType dirId) const
	{
		assert(dirId == 0);
		throw RuntimeError("LongRange::getVectorSize(): unimplemented\n");
	}

	bool connected(SizeType i1,SizeType i2) const
	{
		return true;
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType,SizeType) const
	{
		return 0;
	}

	bool fringe(SizeType,SizeType,SizeType) const
	{
		return true;
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax,SizeType emin,SizeType siteNew2) const
	{
		assert(siteNew2 >= emin);
		SizeType tmp = siteNew2 - emin + smax+1;
		assert(tmp < linSize_);
		return tmp;
	}

	String label() const
	{
		return "LongRange";
	}

	SizeType findReflection(SizeType site) const
	{
		return linSize_ - site -1;
	}

	SizeType length(SizeType i) const
	{
		assert(i==0);
		return linSize_;
	}

	SizeType translate(SizeType site,SizeType dir,SizeType amount) const
	{
		assert(dir==0);

		site+=amount;
		while (site>=linSize_) site -= linSize_;
		return site;
	}

	template<class Archive>
	void serialize(Archive&, const unsigned int)
	{
		throw RuntimeError("LongRange::serialize(): unimplemented\n");
	}

	SizeType memResolv(MemResolv&,
	                   SizeType,
	                   String) const
	{
		throw RuntimeError("LongRange::memResolv(): unimplemented\n");
	}

private:

	SizeType linSize_;
	SizeType orbitals_;
	SizeType maxConnections_;
	MatrixType matrix_;
}; // class LongRange
} // namespace PsimagLite

/*@}*/
#endif // LADDER_H

