/*
Copyright (c) 2009-2013-2018, UT-Battelle, LLC
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

/*! \file WaveFunctionTransfSu2.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_BASE_H
#define WFT_BASE_H
#include "ProgramGlobals.h"
#include "PackIndices.h"

namespace Dmrg {

template<typename DmrgWaveStructType_,typename VectorWithOffsetType_>
class WaveFunctionTransfBase {

public:

	typedef DmrgWaveStructType_ DmrgWaveStructType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef PsimagLite::PackIndices PackIndicesType;

	struct WftOptions {
		typedef typename DmrgWaveStructType::SparseElementType ComplexOrRealType;
		typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

		enum AccelEnum {ACCEL_NONE, ACCEL_PATCHES, ACCEL_BLOCKS, ACCEL_SVD};

		WftOptions(ProgramGlobals::DirectionEnum dir1,
		           PsimagLite::String options,
		           bool f,
		           SizeType c,
		           RealType d)
		    : dir(dir1),
		      twoSiteDmrg(options.find("twositedmrg") != PsimagLite::String::npos),
		      accel((twoSiteDmrg) ? ACCEL_BLOCKS : ACCEL_PATCHES),
		      kronLoadBalance(options.find("KronLoadBalance") != PsimagLite::String::npos),
		      firstCall(f),
		      counter(c),
		      denseSparseThreshold(d)
		{
			if (options.find("wftAccelSvd") != PsimagLite::String::npos)
				accel = ACCEL_SVD;

			if (options.find("wftNoAccel") != PsimagLite::String::npos)
				accel = ACCEL_NONE;

			if (accel == ACCEL_SVD && twoSiteDmrg)
				err("wftAccelSvd not yet supported with twositedmrg\n");
		}

		void read(PsimagLite::IoSelector::In& io, PsimagLite::String label)
		{
			io.read(dir, label + "/dir");
			io.read(twoSiteDmrg, label + "/twoSiteDmrg");
			io.read(accel, label + "/accel");
			io.read(kronLoadBalance, label + "/kronLoadBalance");
			io.read(firstCall, label + "/firstCall");
			io.read(counter, label + "/counter");
			io.read(denseSparseThreshold, label + "/denseSparseThreshold");
		}

		void write(PsimagLite::IoSelector::Out& io, PsimagLite::String label) const
		{
			io.createGroup(label);
			io.write(dir, label + "/dir");
			io.write(twoSiteDmrg, label + "/twoSiteDmrg");
			io.write(accel, label + "/accel");
			io.write(kronLoadBalance, label + "/kronLoadBalance");
			io.write(firstCall, label + "/firstCall");
			io.write(counter, label + "/counter");
			io.write(denseSparseThreshold, label + "/denseSparseThreshold");
		}

		ProgramGlobals::DirectionEnum dir;
		bool twoSiteDmrg;
		AccelEnum accel;
		bool kronLoadBalance;
		bool firstCall;
		SizeType counter;
		RealType denseSparseThreshold;

	private:

		void accelMustBeNone(SizeType x) const
		{
			if (accel != ACCEL_NONE) {
				err("WFTOptions: More than one accel mode specified. Specify none or 1\n");
				return;
			}

			if (x == 0 || twoSiteDmrg) return;

			err("WFTOptions: onesitedmrg only with ACCEL_NONE or ACCEL_PATCHES\n");
		}
	};

	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	virtual void transformVector(VectorWithOffsetType& psiDest,
	                             const VectorWithOffsetType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             const VectorSizeType& nk) const = 0;

	virtual ~WaveFunctionTransfBase() {}

protected:



}; // class WaveFunctionTransfBase
} // namespace Dmrg

/*@}*/
#endif // WFT_BASE_H

