/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.0.0]
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

/*! \file GeometryDirection.h
 *
 *  FIXME DOC. TBW
 */
#ifndef GEOMETRY_DIR_H
#define GEOMETRY_DIR_H
#include <cassert>
#include "BoostSerializationHeaders.h"

namespace PsimagLite {

template<typename ComplexOrRealType,typename GeometryBaseType>
class GeometryDirection {

	typedef Matrix<ComplexOrRealType> MatrixType;
	typedef typename Real<ComplexOrRealType>::Type RealType;

public:

	struct Auxiliary {

		Auxiliary(bool c, SizeType d, SizeType e)
		    : constantValues(c), dirId(d), edof(e)
		{}

		bool constantValues;
		SizeType dirId;
		SizeType edof;
	}; // struct Auxiliary

	template<typename IoInputter>
	GeometryDirection(IoInputter& io,
	                  const Auxiliary& aux,
	                  const GeometryBaseType* geometryFactory)
	    : aux_(aux),
	      geometryBase_(geometryFactory)
	{
		SizeType n = (aux.edof & 2) ? 0 : getVectorSize();
		dataType_ = aux.edof;
		orbitals_ = 1;

		if (aux.edof & 2) {
			geometryBase_->set(rawHoppings_, orbitals_);
			return;
		}

		String connectors = "Connectors";
		String savedPrefix = io.prefix();
		io.prefix() += "dir" + ttos(aux.dirId) + ":";

		if (aux.edof & 1) {
			if (n == 0) n = 1;
			for (SizeType i=0;i<n;i++) {
				MatrixType m;
				String extraString =  (n > 1 && io.version() > 2) ? ttos(i) : "";
				io.read(m, connectors + extraString);
				dataMatrices_.push_back(m);
				if (orbitals_ < m.rows()) orbitals_ = m.rows();
				if (orbitals_ < m.cols()) orbitals_ = m.cols();
			}
		} else {
			io.read(dataNumbers_, connectors);
			if (dataNumbers_.size()!=n) {
				String s(__FILE__);
				s += " " + ttos(dataNumbers_.size()) + " != " + ttos(n) + "\n";
				throw RuntimeError(s.c_str());
			}
		}

		io.prefix() = savedPrefix;
	}

	template<class Archive>
	void serialize(Archive&, const unsigned int)
	{}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType,
	                   String) const
	{
		return 0;
	}

	ComplexOrRealType operator()(SizeType i,
	                             SizeType edof1,
	                             SizeType j,
	                             SizeType edof2) const
	{
		SizeType ii = edof1 + i*orbitals_;
		SizeType jj = edof2 + j*orbitals_;
		if (dataType_ & 2) {
			assert((dataType_&1) || (edof1 == 0 && edof2 == 0));
			return rawHoppings_(ii,jj);
		}

		SizeType h = (constantValues()) ? 0 : geometryBase_->handle(i,j);

		bool isMatrix = (dataType_&1);
		if (!isMatrix) {
			assert(dataNumbers_.size()>h);
			return dataNumbers_[h];
		}

		if (ii < jj) return PsimagLite::conj(operator()(j,edof2,i,edof1));

		assert(dataMatrices_.size()>h);

		bool b = (dataMatrices_[h].rows()>edof1 &&
		          dataMatrices_[h].cols()>edof2);

		assert(b ||  (dataMatrices_[h].rows()>edof2 &&
		              dataMatrices_[h].cols()>edof1));

		ComplexOrRealType tmp = (b) ?
		            dataMatrices_[h](edof1,edof2) : dataMatrices_[h](edof2,edof1);
		int signChange = geometryBase_->signChange(i,j);
		return tmp * static_cast<RealType>(signChange);
	}

	SizeType orbitals() const { return orbitals_; }

	bool constantValues() const
	{
		return aux_.constantValues;
	}

	friend std::ostream& operator<<(std::ostream& os, const Auxiliary& a)
	{
		os<<"constantValues="<<a.constantValues<<"\n";
		os<<"dirId="<<a.dirId<<"\n";
		os<<"edof="<<a.edof<<"\n";

		return os;
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const GeometryDirection& gd)
	{
		os<<"#GeometryDirectionAuxiliary\n";
		os<<gd.aux_;
		bool isMatrix = (gd.dataType_&1);
		if (!isMatrix) {
			os<<"#GeometryNumbersSize="<<gd.dataNumbers_.size()<<"\n";
			os<<"#GeometryNumbers=";
			for (SizeType i=0;i<gd.dataNumbers_.size();i++) {
				os<<gd.dataNumbers_[i]<<" ";
			}
			os<<"\n";
		} else {
			os<<"#GeometryMatrixSize="<<gd.dataMatrices_.size()<<"\n";
			for (SizeType i=0;i<gd.dataMatrices_.size();i++)
				os<<gd.dataMatrices_[i];
		}
		return os;
	}

private:

	SizeType getVectorSize()
	{
		return (aux_.constantValues) ? 1 : geometryBase_->getVectorSize(aux_.dirId);
	}

	Auxiliary aux_;
	SizeType dataType_;
	SizeType orbitals_;
	const GeometryBaseType* geometryBase_;
	typename Vector<ComplexOrRealType>::Type dataNumbers_;
	typename Vector<MatrixType>::Type dataMatrices_;
	MatrixType rawHoppings_;
}; // class GeometryDirection
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_DIR_H

