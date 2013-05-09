// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup PsimagLite */
/*@{*/

/*! \file GeometryDirection.h
 *
 *  FIXME DOC. TBW
 */
#ifndef GEOMETRY_DIR_H
#define GEOMETRY_DIR_H
#include <cassert>
#include "String.h"

namespace PsimagLite {
	
	template<typename RealType,typename GeometryFactoryType>
	class GeometryDirection {

		typedef Matrix<RealType> MatrixType;

		public:

			enum {NUMBERS,MATRICES};

			template<typename IoInputter>
			GeometryDirection(IoInputter& io,size_t dirId,size_t edof,
					  const String& options,
					  const GeometryFactoryType& geometryFactory)
			: dirId_(dirId),geometryFactory_(&geometryFactory)
			{
				size_t n = getVectorSize(options);
				dataType_ = edof;
				if (edof==NUMBERS) {
					io.read(dataNumbers_,"Connectors");
					if (dataNumbers_.size()!=n) {
						String s(__FILE__);
						s += " " + ttos(dataNumbers_.size()) + " != " + ttos(n) + "\n";
					 	throw std::runtime_error(s.c_str());
					}
				} else {
					for (size_t i=0;i<n;i++) {
						MatrixType m;
						io.readMatrix(m,"Connectors");
// 						if (m.n_row()!=edof || m.n_col()!=edof)
// 							throw std::runtime_error("GeometryDirection\n");
						dataMatrices_.push_back(m);
					}
				}
			}

// 			const RealType& operator()(size_t i,size_t j) const
// 			{
// 				if (dataType_==NUMBERS) return dataNumbers_[0];
// 				return dataMatrices_[0](i,j);
// 			}

			RealType operator()(size_t i,size_t edof1,size_t j,size_t edof2) const
			{
				size_t h = (constantValues()) ? 0 : geometryFactory_->handle(i,j);
				
				if (dataType_==NUMBERS) {
					assert(dataNumbers_.size()>h);
					return dataNumbers_[h];
				}
				assert(dataMatrices_.size()>h);
				
				bool b = (dataMatrices_[h].n_row()>edof1 &&
				          dataMatrices_[h].n_col()>edof2);
				
				assert(b ||  (dataMatrices_[h].n_row()>edof2 &&
				              dataMatrices_[h].n_col()>edof1));

				RealType tmp = (b) ? dataMatrices_[h](edof1,edof2) : dataMatrices_[h](edof2,edof1);
				int signChange = geometryFactory_->signChange(i,j);
				return tmp * signChange;
			}

			size_t nRow() const
			{
				return (dataType_==NUMBERS) ? 1 : dataMatrices_[0].n_row();
			}
			
			size_t nCol() const
			{
				return (dataType_==NUMBERS) ? 1 : dataMatrices_[0].n_col();
			}
			
			size_t size() const
			{
				return (dataType_==NUMBERS) ? dataNumbers_.size() : dataMatrices_.size();
			}

			bool constantValues() const
			{
				return (size()==1) ? true : false;
			}

			template<typename RealType_,typename GeometryFactoryType_>
			friend std::ostream& operator<<(std::ostream& os,const GeometryDirection<RealType_,GeometryFactoryType_>& gd);

		private:

			size_t getVectorSize(const String& s)
			{
				if (s.find("ConstantValues")!=String::npos)
					return 1;

				return geometryFactory_->getVectorSize(dirId_);
			}
			
			size_t dirId_;
			const GeometryFactoryType* geometryFactory_;
			size_t dataType_;
			typename Vector<RealType>::Type dataNumbers_;
			typename Vector<MatrixType>::Type dataMatrices_;
	}; // class GeometryDirection

	template<typename RealType,typename GeometryFactoryType>
	std::ostream& operator<<(std::ostream& os,const GeometryDirection<RealType,GeometryFactoryType>& gd)
	{
		os<<"#GeometrydirId="<<gd.dirId_<<"\n";
		if (gd.dataType_==GeometryDirection<RealType,GeometryFactoryType>::NUMBERS) {
			os<<"#GeometryNumbersSize="<<gd.dataNumbers_.size()<<"\n";
			os<<"#GeometryNumbers=";
			for (size_t i=0;i<gd.dataNumbers_.size();i++) {
				os<<gd.dataNumbers_[i]<<" ";
			}
			os<<"\n";
		} else {
			os<<"#GeometryMatrixSize="<<gd.dataMatrices_.size()<<"\n";
			for (size_t i=0;i<gd.dataMatrices_.size();i++) 
				os<<gd.dataMatrices_[i];
		}
		return os;
	}
} // namespace PsimagLite 

/*@}*/
#endif // GEOMETRY_DIR_H
