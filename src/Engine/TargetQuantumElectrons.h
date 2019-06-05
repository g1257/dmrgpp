/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

/*! \file TargetQuantumElectrons.h
 *
 *
 */
#ifndef TargetQuantumElectrons_H
#define TargetQuantumElectrons_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {
//! Hubbard Model Parameters
template<typename RealType, typename QnType>
class TargetQuantumElectrons {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename QnType::PairSizeType PairSizeType;
	typedef typename QnType::VectorQnType VectorQnType;

	template<typename IoInputType>
	TargetQuantumElectrons(IoInputType& io)
	    : totalNumberOfSites_(0),
	      isSu2_(false)
	{
		io.readline(totalNumberOfSites_, "TotalNumberOfSites=");
		int tmp = 0;
		try {
			io.readline(tmp,"UseSu2Symmetry=");
		} catch (std::exception&) {}

		isSu2_ = (tmp > 0);

		bool hasNqns = false;
		SizeType nqns = 0;
		try {
			io.readline(nqns, "NumberOfTargetQns=");
			hasNqns = true;
		} catch (std::exception&) {}

		if (!hasNqns) {
			readOneTarget(io, "");
			return;
		}

		for (SizeType i = 0; i < nqns; ++i)
			readOneTarget(io, ttos(i));
	}

	SizeType sizeOfOther() const
	{
		const SizeType n = vqn_.size();
		if (n == 0) return 0;
		const SizeType answer = vqn_[0].other.size();
		for (SizeType i = 1; i < n; ++i) {
			if (vqn_[i].other.size() == answer) continue;
			err("sizeOfOther must be the same for all target qns\n");
		}

		return answer;
	}

	const QnType& qn(SizeType ind) const
	{
		assert(ind < vqn_.size());
		return vqn_[ind];
	}

	void updateQuantumSector(VectorQnType& quantumSector,
	                         SizeType sites,
	                         ProgramGlobals::DirectionEnum direction,
	                         SizeType step,
	                         const VectorQnType& adjustQuantumNumbers) const
	{
		const SizeType maxSites = totalNumberOfSites_;

		if (direction == ProgramGlobals::DirectionEnum::INFINITE &&
		        sites < maxSites &&
		        adjustQuantumNumbers.size() > step) {
			if (quantumSector.size() != 1)
				err("adjustQuantumNumbers only with single target\n");
			quantumSector[0] = adjustQuantumNumbers[step];
			return;
		} else {
			quantumSector = vqn_;
		}

		const SizeType n = quantumSector.size();

		for (SizeType i = 0; i < n; ++i)
			quantumSector[i].scale(sites,
			                       totalNumberOfSites_,
			                       direction,
			                       isSu2_);
	}

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/TargetQuantumElectrons";
		io.createGroup(label);
		io.write(label + "/TotalNumberOfSites", totalNumberOfSites_);
		io.write(label + "/isSu2", isSu2_);
		vqn_.write(label + "/qn", io);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream &os,
	                                const TargetQuantumElectrons& p)
	{
		if (p.vqn_.size() == 0) return os;
		const QnType& qn = p.vqn_[0];
		os<<"TargetElectronsTotal="<<p.totalElectrons<<"\n";
		os<<"TargetOther="<<p.other<<"\n";
		if (p.isSu2)
			os<<"TargetSpinTimesTwo="<<p.twiceJ<<"\n";
		if (p.vqn_.size() > 1)
			os<<"FIXME TODO: More than one qn found\n";
		return os;
	}

private:

	TargetQuantumElectrons(const TargetQuantumElectrons&);

	TargetQuantumElectrons& operator=(const TargetQuantumElectrons&);

	template<typename IoInputType>
	void readOneTarget(IoInputType& io,
	                   const PsimagLite::String label)
	{
		QnType qn(QnType::zero());
		VectorSizeType qnOther;
		const bool allowUpDown = true;

		PsimagLite::String msg("TargetQuantumElectrons: ");
		bool hasTwiceJ = false;
		try {
			io.readline(qn.jmPair.first, "TargetSpinTimesTwo" + label + "=");
			hasTwiceJ = true;
		} catch (std::exception&) {}

		SizeType ready = 0;
		if (allowUpDown) {
			SizeType electronsUp = 0;
			SizeType electronsDown = 0;
			try {
				io.readline(electronsUp,"TargetElectronsUp" + label + "=");
				io.readline(electronsDown,"TargetElectronsDown" + label + "=");
				SizeType tmp = electronsUp + electronsDown;
				qn.oddElectrons = (tmp & 1);
				qnOther.push_back(tmp);
				qnOther.push_back(electronsUp);
				ready = 2;
			} catch (std::exception&) {}
		}

		try {
			SizeType tmp = 0;
			io.readline(tmp, "TargetElectronsTotal" + label + "=");
			qn.oddElectrons = (tmp & 1);
			qnOther.push_back(tmp);
			ready++;
		} catch (std::exception&) {}

		try {
			SizeType szPlusConst = 0;
			io.readline(szPlusConst,"TargetSzPlusConst" + label + "=");
			qnOther.push_back(szPlusConst);
		} catch (std::exception&) {}

		if (ready == 3) {
			msg += "Provide either up/down or total/sz but not both.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		bool flag = false;
		try {
			int dummy = 0;
			io.readline(dummy, "TargetExtra" + label + "=");
			flag = true;
		} catch (std::exception&) {}

		if (flag) err("Instead of TargetExtra" + label + "= please use a vector\n");

		try {
			VectorSizeType extra;
			io.read(extra,"TargetExtra" + label);
			for (SizeType i = 0; i < extra.size(); ++i)
				qnOther.push_back(extra[i]);
		} catch (std::exception&) {}

		qn.other.fromStdVector(qnOther);

		if (isSu2_ && !hasTwiceJ) {
			msg += "Please provide TargetSpinTimesTwo when running with SU(2).\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (isSu2_)
			qn.oddElectrons = (totalNumberOfSites_ & 1);

		vqn_.push_back(qn);
	}

	SizeType totalNumberOfSites_;
	bool isSu2_;
	VectorQnType vqn_;
};
} // namespace Dmrg

/*@}*/
#endif

