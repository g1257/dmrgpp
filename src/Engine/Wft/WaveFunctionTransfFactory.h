/*
Copyright (c) 2009-2012-2018, UT-Battelle, LLC
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

/*! \file WaveFunctionTransfFactory.h
 *
 *  This class implements the wave function transformation factory,
 *  see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_FACTORY_H
#define WFT_FACTORY_H
#include "Utils.h"
#include "ProgressIndicator.h"
#include "WaveFunctionTransfLocal.h"
#include "WaveFunctionTransfSu2.h"
#include "WaveStructCombined.h"
#include "Io/IoSelector.h"
#include "Random48.h"

namespace Dmrg {
template<typename LeftRightSuperType,typename VectorWithOffsetType_>
class WaveFunctionTransfFactory {

	typedef PsimagLite::IoSelector IoType;

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef WaveStructCombined<LeftRightSuperType> WaveStructCombinedType;
	typedef typename WaveStructCombinedType::VectorVectorRealType VectorVectorRealType;
	typedef typename WaveStructCombinedType::VectorMatrixType VectorMatrixType;
	typedef typename WaveStructCombinedType::VectorQnType VectorQnType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef WaveFunctionTransfBase<WaveStructCombinedType,VectorWithOffsetType>
	WaveFunctionTransfBaseType;
	typedef WaveFunctionTransfLocal<WaveStructCombinedType,VectorWithOffsetType>
	WaveFunctionTransfLocalType;
	typedef WaveFunctionTransfSu2<WaveStructCombinedType,VectorWithOffsetType>
	WaveFunctionTransfSu2Type;
	typedef typename WaveFunctionTransfBaseType::WftOptionsType WftOptionsType;
	typedef typename WaveStructCombinedType::WaveStructSvdType WaveStructSvdType;

	template<typename SomeParametersType>
	WaveFunctionTransfFactory(SomeParametersType& params)
	    : isEnabled_(!(params.options.find("nowft")!=PsimagLite::String::npos)),
	      wftOptions_(ProgramGlobals::DirectionEnum::INFINITE,
	                  params.options,
	                  true,
	                  true,
	                  params.denseSparseThreshold,
	                  params.options.find("blasNotThreadSafe") == PsimagLite::String::npos),
	      progress_("WaveFunctionTransf"),
	      filenameIn_(params.checkpoint.filename()),
	      filenameOut_(params.filename),
	      wftImpl_(0),
	      rng_(3433117),
	      noLoad_(false),
	      save_(params.options.find("noSaveWft") == PsimagLite::String::npos)
	{
		if (!isEnabled_) return;

		bool b = (params.options.find("restart")!=PsimagLite::String::npos ||
		        params.autoRestart);

		if (b) {
			if (params.options.find("noloadwft")!=PsimagLite::String::npos)
				noLoad_=true;
			else
				read();
		} else {
			if (params.options.find("noloadwft")!=PsimagLite::String::npos) {
				PsimagLite::String str("Error: noloadwft needs restart or checkpoint\n");
				throw PsimagLite::RuntimeError(str.c_str());
			}
		}

		if (BasisType::useSu2Symmetry())
			wftImpl_=new WaveFunctionTransfSu2Type(waveStructCombined_,
			                                       wftOptions_);
		else
			wftImpl_=new WaveFunctionTransfLocalType(waveStructCombined_,
			                                         wftOptions_);
	}

	~WaveFunctionTransfFactory()
	{
		if (!isEnabled_) return;
		IoType::Out ioOut(filenameOut_, IoType::ACC_RDW);
		write(ioOut);
		delete wftImpl_;
	}

	void setStage(ProgramGlobals::DirectionEnum stage)
	{
		if (stage == wftOptions_.dir) return;
		wftOptions_.dir = stage;
		wftOptions_.bounce = true;
	}

	void triggerOn()
	{
		bool allow = false;
		switch (wftOptions_.dir) {
		case ProgramGlobals::DirectionEnum::INFINITE:
			allow = false;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_SYSTEM:
		case ProgramGlobals::DirectionEnum::EXPAND_ENVIRON:
			allow = true;
		}

		// FIXME: Must check the below change when using SU(2)!!
		//if (m<0) allow = false; // isEnabled_=false;

		if (noLoad_) allow = false;

		if (!isEnabled_ || !allow) return;

		waveStructCombined_.beforeWft(wftOptions_.dir,
		                              wftOptions_.twoSiteDmrg,
		                              wftOptions_.bounce);
		PsimagLite::OstringStream msg;
		msg<<"Window open, ready to transform vectors";
		progress_.printline(msg,std::cout);
	}

	// FIXME: change name to transformVector
	void setInitialVector(VectorWithOffsetType& dest,
	                      const VectorWithOffsetType& src,
	                      const LeftRightSuperType& lrs,
	                      const VectorSizeType& nk) const
	{
		bool allow=false;
		switch (wftOptions_.dir) {
		case ProgramGlobals::DirectionEnum::INFINITE:
			allow=false;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_SYSTEM:
			allow=true;

		case ProgramGlobals::DirectionEnum::EXPAND_ENVIRON:
			allow=true;
		}

		// FIXME: Must check the below change when using SU(2)!!
		//if (m<0) allow = false; // isEnabled_=false;

		if (noLoad_) allow = false;

		if (isEnabled_ && allow) {
#ifndef NDEBUG
			RealType eps = 1e-12;
			RealType x = norm(src);
			bool b = (x<eps);
			if (b) std::cerr<<"norm="<<x<<"\n";
			assert(!b);
#endif
			createVector(dest,src,lrs,nk);
		} else {
			createRandomVector(dest);
		}
	}

	void triggerOff(const LeftRightSuperType& lrs)
	{
		bool allow=false;
		switch (wftOptions_.dir) {
		case ProgramGlobals::DirectionEnum::INFINITE:
			allow=false;
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_SYSTEM:
			allow=true;

		case ProgramGlobals::DirectionEnum::EXPAND_ENVIRON:
			allow=true;
		}

		// FIXME: Must check the below change when using SU(2)!!
		//if (m<0) allow = false; // isEnabled_=false;

		if (noLoad_) allow = false;

		if (!isEnabled_ || !allow) return;
		afterWft(lrs);
		PsimagLite::OstringStream msg;
		msg<<"Window closed, no more transformations, please";
		progress_.printline(msg,std::cout);
	}

	void createRandomVector(VectorWithOffsetType& y) const
	{
		for (SizeType jj=0;jj<y.sectors();jj++) {
			SizeType j = y.sector(jj);
			createRandomVector(y, j);
		}

		if (!isEnabled_) return; // don't make noise unless enabled
		PsimagLite::OstringStream msg;
		msg<<"Yes, I'm awake, but there's nothing heavy to do now";
		progress_.printline(msg,std::cout);
	}

	void createRandomVector(VectorWithOffsetType& y, SizeType i0) const
	{
		SizeType total = y.effectiveSize(i0);
		typename VectorWithOffsetType::value_type tmp;
		RealType atmp=0;
		for (SizeType i=0;i<total;i++) {
			myRandomT(tmp);
			y.fastAccess(i0,i)=tmp;
			atmp += PsimagLite::real(tmp*PsimagLite::conj(tmp));
		}

		assert(fabs(atmp)>1e-10);
		atmp = 1.0 / sqrt (atmp);
		for (SizeType i=0;i<total;i++) y.fastAccess(i0,i) *= atmp;

	}

	void push(const BlockDiagonalMatrixType& transform,
	          ProgramGlobals::DirectionEnum direction,
	          const LeftRightSuperType& lrs,
	          const VectorMatrixType& vts,
	          const VectorVectorRealType& s,
	          const VectorQnType& qns)
	{
		if (!isEnabled_) return;

		waveStructCombined_.push(transform, direction, vts, s, qns, wftOptions_.dir);
		waveStructCombined_.setLrs(lrs);
		PsimagLite::OstringStream msg;
		msg<<"OK, pushing option="<<ProgramGlobals::toString(direction);
		msg<<" and stage="<<ProgramGlobals::toString(wftOptions_.dir);
		progress_.printline(msg,std::cout);

		if (noLoad_) {
			SizeType center = computeCenter(lrs,direction);
			updateNoLoad(lrs,center);
		}
	}

	const BlockDiagonalMatrixType& getTransform(ProgramGlobals::SysOrEnvEnum dir) const
	{
		return waveStructCombined_.getTransform(dir);
	}

	SizeType size(ProgramGlobals::SysOrEnvEnum sysOrEnv) const
	{
		return waveStructCombined_.size(sysOrEnv);
	}

	bool isEnabled() const { return isEnabled_; }

	const WftOptionsType options() const { return wftOptions_; }

	void write(PsimagLite::IoSelector::Out& ioMain)
	{
		if (!isEnabled_) return;
		if (!save_) return;

		PsimagLite::String label = "Wft";
		writePartial(ioMain, label);
		waveStructCombined_.write(ioMain, label + "/WaveStructCombined");
	}

	void write(PsimagLite::IoSelector::Out& ioMain) const
	{
		if (!isEnabled_) return;
		if (!save_) return;

		PsimagLite::String label = "Wft";
		writePartial(ioMain, label);
		waveStructCombined_.write(ioMain, label + "/WaveStructCombined");
	}

private:

	void writePartial(PsimagLite::IoSelector::Out& ioMain,
	                  PsimagLite::String prefix) const
	{
		assert(isEnabled_);
		assert(save_);

		ioMain.createGroup(prefix);
		ioMain.write(isEnabled_, prefix + "/isEnabled");
		wftOptions_.write(ioMain, prefix + "/WftOptions");
	}

	void read()
	{
		if (!isEnabled_)
			throw PsimagLite::RuntimeError("WFT::read(...) called but wft is disabled\n");

		PsimagLite::IoSelector::In ioMain(filenameIn_);
		PsimagLite::String label = "Wft";
		ioMain.read(isEnabled_, label + "/isEnabled");
		wftOptions_.read(ioMain, label + "/WftOptions");
		waveStructCombined_.read(ioMain, label + "/WaveStructCombined");
		ioMain.close();
	}

	void myRandomT(std::complex<RealType> &value) const
	{
		value = std::complex<RealType>(rng_() - 0.5, rng_() - 0.5);
	}

	void myRandomT(RealType &value) const
	{
		value = rng_() - 0.5;
	}

	void afterWft(const LeftRightSuperType& lrs)
	{
		waveStructCombined_.setLrs(lrs);
		waveStructCombined_.afterWft(wftOptions_.dir);
		wftOptions_.firstCall = false;
		wftOptions_.bounce = false;
	}

	void createVector(VectorWithOffsetType& psiDest,
	                  const VectorWithOffsetType& psiSrc,
	                  const LeftRightSuperType& lrs,
	                  const VectorSizeType& nk) const
	{

		RealType norm1 = norm(psiSrc);
		if (norm1 < 1e-5)
			std::cerr<<"WFT Factory: norm1 = " << norm1 << " < 1e-5\n";

		wftImpl_->transformVector(psiDest, psiSrc, lrs, nk);

		RealType norm2 = norm(psiDest);
		PsimagLite::OstringStream msg;
		msg<<"Transformation completed ";
		if (fabs(norm1-norm2)>1e-5) {
			msg<<"WARNING: orig. norm= "<<norm1<<" resulting norm= "<<norm2;
		}

		if (norm2 < 1e-5)
			std::cerr << "WFT Factory: norm2 = " << norm2 << " < 1e-5\n";

		progress_.printline(msg,std::cout);
	}

	SizeType computeCenter(const LeftRightSuperType& lrs,
	                       ProgramGlobals::DirectionEnum direction) const
	{
		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			SizeType total = lrs.left().block().size();
			assert(total>0);
			total--;
			return lrs.left().block()[total];
		}

		return lrs.right().block()[0];
	}

	void updateNoLoad(const LeftRightSuperType& lrs,SizeType center)
	{
		sitesSeen_.push_back(center);
		SizeType numberOfSites = lrs.super().block().size();
		if (checkSites(numberOfSites)) {
			noLoad_=false;
			PsimagLite::OstringStream msg;
			msg<<" now available";
			progress_.printline(msg,std::cout);
		}
	}

	bool checkSites(SizeType numberOfSites) const
	{
		assert(numberOfSites>0);
		for (SizeType i=1;i<numberOfSites-1;i++) {
			bool seen = (std::find(sitesSeen_.begin(),
			                       sitesSeen_.end(),
			                       i) != sitesSeen_.end());
			if (!seen) return false;
		}

		return true;
	}

	bool isEnabled_;
	WftOptionsType wftOptions_;
	PsimagLite::ProgressIndicator progress_;
	PsimagLite::String filenameIn_;
	PsimagLite::String filenameOut_;
	WaveStructCombinedType waveStructCombined_;
	WaveFunctionTransfBaseType* wftImpl_;
	PsimagLite::Random48<RealType> rng_;
	bool noLoad_;
	const bool save_;
	VectorSizeType sitesSeen_;
}; // class WaveFunctionTransformation
} // namespace Dmrg

/*@}*/
#endif

