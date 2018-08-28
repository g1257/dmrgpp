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

/*! \file WaveStructSvd.h
 *
 *  DOC NEEDED FIXME (This file should go in Wft/ directory perhaps)
 */
#ifndef WAVE_STRUCT_SVD_H
#define WAVE_STRUCT_SVD_H
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg {

template<typename LeftRightSuperType_>
struct WaveStructSvd {

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::QnType QnType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename BasisWithOperatorsType::VectorQnType VectorQnType;

	WaveStructSvd() {}

	WaveStructSvd(const BlockDiagonalMatrixType& u,
	                   const VectorMatrixType& vts,
	                   const VectorVectorRealType& s,
	                   const VectorQnType& qns)
	    : u_(u), vts_(vts), s_(s), qns_(qns)
	{}

	WaveStructSvd& operator=(const WaveStructSvd& other)
	{
		u_ = other.u_;
		vts_ = other.vts_;
		s_ = other.s_;
		qns_ = other.qns_;
		return *this;
	}

	const BlockDiagonalMatrixType& u() const { return u_; }

	const VectorMatrixType& vts() const { return vts_; }

	const VectorVectorRealType& s() const { return s_; }

	const VectorQnType& qns() const { return qns_; }

	void read(PsimagLite::IoNg::In& io, PsimagLite::String prefix)
	{
		io.read(*this, prefix);
	}

	void read(PsimagLite::String prefix, PsimagLite::IoNgSerializer& io)
	{
		u_.read(prefix + "/u", io);
		io.read(vts_, prefix + "/vts");
		io.read(s_, prefix + "/s");
		QnType::readVector(qns_, prefix + "/qns", io);
	}

	void write(PsimagLite::IoNg::Out& io, PsimagLite::String prefix) const
	{
		io.createGroup(prefix);
		io.write(u_, prefix + "/u");
		io.write(vts_, prefix + "/vts");
		io.write(s_, prefix + "/s");
		io.write(qns_, prefix + "/qns");
	}

	void write(PsimagLite::String prefix, PsimagLite::IoNgSerializer& io) const
	{
		io.createGroup(prefix);
		u_.write(prefix + "/u", io);
		io.write(prefix + "/vts", vts_);
		io.write(prefix + "/s", s_);
		io.write(prefix + "/qns", qns_);
	}

	void clear()
	{
		u_.clear();
		vts_.clear();
		s_.clear();
		qns_.clear();
	}

private:

	BlockDiagonalMatrixType u_;
	VectorMatrixType vts_;
    VectorVectorRealType s_;
    VectorQnType qns_;
}; // struct WaveStructSvd

} // namespace Dmrg 

/*@}*/
#endif // WAVE_STRUCT_SVD_H
