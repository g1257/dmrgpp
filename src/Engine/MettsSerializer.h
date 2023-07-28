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

/*! \file MettsSerializer.h
 *
 *  Serialize time data
 */
#ifndef METTS_SERIALIZER_H
#define METTS_SERIALIZER_H

#include "Io/IoSelector.h"
#include "TypeToString.h"

namespace Dmrg
{

template <typename VectorType>
class MettsSerializer
{

	typedef typename VectorType::value_type VectorElementType;
	typedef typename PsimagLite::Real<VectorElementType>::Type RealType;

public:

	// Unfortunately we need a default ctor
	// to build an array of these
	MettsSerializer() { }

	MettsSerializer(RealType currentBeta,
	    SizeType site,
	    const typename PsimagLite::Vector<VectorType>::Type& targetVectors)
	    : currentBeta_(currentBeta)
	    , site_(site)
	    , targetVectors_(targetVectors)
	{
	}

	MettsSerializer(typename PsimagLite::IoSelector::In& io)
	{
		RealType x = 0;
		PsimagLite::String s = "BETA";

		io.read(x, s);

		if (x < 0)
			throw PsimagLite::RuntimeError("MettsSerializer:: time cannot be negative\n");

		currentBeta_ = x;

		s = "TargetCentralSite";
		int xi = 0;
		io.read(xi, s);
		if (xi < 0)
			throw PsimagLite::RuntimeError("MettsSerializer:: site cannot be negative\n");

		site_ = xi;

		s = "TNUMBEROFVECTORS";
		io.read(xi, s);
		if (xi <= 0)
			throw PsimagLite::RuntimeError("MettsSerializer:: n. of vectors must be positive\n");

		targetVectors_.resize(xi);
		for (SizeType i = 0; i < targetVectors_.size(); i++) {
			s = "targetVector" + ttos(i);
			targetVectors_[i].read(io, s);
		}
	}

	SizeType size(SizeType i = 0) const
	{
		return targetVectors_[i].size();
	}

	RealType beta() const { return currentBeta_; }

	SizeType site() const
	{
		return site_;
	}

	const VectorType& vector(SizeType i = 0) const
	{
		return targetVectors_[i];
	}

	template <typename IoOutputter>
	void write(IoOutputter& io,
	    typename PsimagLite::EnableIf<
		PsimagLite::IsOutputLike<IoOutputter>::True,
		int>::Type
	    = 0) const
	{
		PsimagLite::String s = "BETA=" + ttos(currentBeta_);
		io.printline(s);
		s = "TargetCentralSite=" + ttos(site_);
		io.printline(s);
		s = "TNUMBEROFVECTORS=" + ttos(targetVectors_.size());
		io.printline(s);
		for (SizeType i = 0; i < targetVectors_.size(); i++) {
			PsimagLite::String label = "targetVector" + ttos(i) + "_" + ttos(currentBeta_);
			targetVectors_[i].write(io, label);
		}
	}

private:

	RealType currentBeta_;
	SizeType site_;
	typename PsimagLite::Vector<VectorType>::Type targetVectors_;
}; // class MettsSerializer
} // namespace Dmrg

/*@}*/
#endif // METTS_SERIALIZER_H
