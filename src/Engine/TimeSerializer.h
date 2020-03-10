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
/** \ingroup DMRG */
/*@{*/

/*! \file TimeSerializer.h
 *
 *  Serialize time data
 */
#ifndef TIME_SERIAL_H
#define TIME_SERIAL_H

#include "Io/IoSelector.h"
#include "TypeToString.h"
#include "Vector.h"
#include "StageEnum.h"

namespace Dmrg {

template<typename VectorType>
class TimeSerializer {

	typedef typename VectorType::value_type VectorElementType;
	typedef typename PsimagLite::Real<VectorElementType>::Type RealType;

public:

	typedef typename PsimagLite::Vector<StageEnum>::Type VectorStageEnumType;

	TimeSerializer(SizeType currentTimeStep,
	               RealType currentTime,
	               SizeType site,
	               const typename PsimagLite::Vector<VectorType>::Type& targetVectors,
	               const VectorStageEnumType& stages,
	               PsimagLite::String name)
	    : currentTimeStep_(currentTimeStep),
	      currentTime_(currentTime),
	      site_(site),
	      targetVectors_(targetVectors),
	      stages_(stages),
	      name_(name)
	{}

	TimeSerializer(typename PsimagLite::IoSelector::In& io,
	               PsimagLite::String prefix)
	{
		prefix += "/TimeSerializer/";

		PsimagLite::String s = prefix + "CurrentTimeStep";
		io.read(currentTimeStep_, s);

		s = prefix + "Time";
		io.read(currentTime_, s);

		s = prefix + "TargetCentralSite";
		int xi = 0;
		io.read(xi, s);
		if (xi < 0)
			err("TimeSerializer:: site cannot be negative\n");
		site_ = xi;

		s = prefix + "TNUMBEROFVECTORS";
		io.read(xi, s);
		if (xi <= 0)
			err("TimeSerializer:: n. of vectors must be positive\n");
		targetVectors_.resize(xi);
		for (SizeType i = 0; i < targetVectors_.size(); ++i) {
			s = prefix + "targetVector"+ttos(i);
			targetVectors_[i].read(io, s);
		}

		s = prefix + "Stages";
		io.read(stages_, s);

		try {
			io.read(name_, prefix + "Name");
		} catch (...) {
			// reading an old file, set name to LEGACY
			name_ = "LEGACY";
		}
	}

	void write(PsimagLite::IoSelector::Out& io, PsimagLite::String prefix) const
	{
		prefix += "/TimeSerializer";
		io.createGroup(prefix);
		prefix += "/";

		io.write(currentTime_, prefix + "Time");
		io.write(currentTimeStep_, prefix + "CurrentTimeStep");
		io.write(site_, prefix + "TargetCentralSite");
		io.write(targetVectors_.size(), prefix + "TNUMBEROFVECTORS");

		for (SizeType i=0;i<targetVectors_.size();i++) {
			PsimagLite::String label = "targetVector" + ttos(i);
			targetVectors_[i].write(io, prefix + label);
		}

		io.write(stages_, prefix + "Stages");
		io.write(name_, prefix + "Name");
	}

	SizeType numberOfVectors() const
	{
		return  targetVectors_.size();
	}

	SizeType currentTimeStep() const { return currentTimeStep_; }

	RealType time() const { return currentTime_; }

	SizeType site() const { return  site_; }

	PsimagLite::String name() const { return name_; }

	const VectorType& vector(SizeType i) const
	{
		if (i < targetVectors_.size())
			return targetVectors_[i];
		throw PsimagLite::RuntimeError("Not so many time vectors\n");
	}

	const VectorStageEnumType& stages() const { return stages_; }

private:

	SizeType currentTimeStep_;
	RealType currentTime_;
	SizeType site_;
	typename PsimagLite::Vector<VectorType>::Type targetVectors_;
	VectorStageEnumType stages_;
	PsimagLite::String name_;
}; // class TimeSerializer
} // namespace Dmrg 

/*@}*/
#endif
