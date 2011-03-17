// BEGIN LICENSE BLOCK
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file GeometryDirection.h
 *
 *  FIXME DOC. TBW
 */
#ifndef GEOMETRY_DIR_H
#define GEOMETRY_DIR_H



namespace Dmrg {
	
	template<typename RealType,typename GeometryFactoryType>
	class GeometryDirection {
			typedef PsimagLite::Matrix<RealType> MatrixType;
			enum {NUMBERS,MATRICES};

		public:

			template<typename IoInputter>
			GeometryDirection(IoInputter& io,size_t dirId,size_t edof,
					  const std::string& options,
					  const GeometryFactoryType& geometryFactory)
			: dirId_(dirId),geometryFactory_(&geometryFactory)
			{
				size_t n = getVectorSize(options);
				//std::cerr<<"vectorsize="<<n<<"\n";
				if (edof==1) {
					io.read(dataNumbers_,"Connectors");
					dataType_ = NUMBERS;
					if (dataNumbers_.size()!=n)
						throw std::runtime_error("GeometryDirection Numbers\n");
				} else {
					for (size_t i=0;i<n;i++) {
						MatrixType m;
						io.readMatrix(m,"Connectors");
						if (m.n_row()!=edof || m.n_col()!=edof)
							throw std::runtime_error("GeometryDirection\n");
						dataMatrices_.push_back(m);
					}
					dataType_ = MATRICES;
				}
			}

			const RealType& operator()(size_t i,size_t j) const
			{
				if (dataType_==NUMBERS) return dataNumbers_[0];
				return dataMatrices_[0](i,j);
			}

			const RealType& operator()(size_t i,size_t edof1,size_t j,size_t edof2) const
			{
				
				size_t h = geometryFactory_->handle(i,j);
				if (dataType_==NUMBERS) return dataNumbers_[h];
				return dataMatrices_[h](edof1,edof2);
			}

			size_t size() const
			{
				if (dataType_==NUMBERS) return dataNumbers_.size();
				return dataMatrices_.size();
			}

			bool constantValues() const
			{
				if (size()==1) return true;
				return false;
			}


			template<typename RealType_,typename GeometryFactoryType_>
			friend std::ostream& operator<<(std::ostream& os,const GeometryDirection<RealType_,GeometryFactoryType_>& gd);

		private:

			size_t getVectorSize(const std::string& s)
			{
				if (s.find("ConstantValues")!=std::string::npos)
					return 1;

				return geometryFactory_->getVectorSize(dirId_);
			}
			
			size_t dirId_;
			const GeometryFactoryType* geometryFactory_;
			size_t dataType_;
			std::vector<RealType> dataNumbers_;
			std::vector<MatrixType> dataMatrices_;
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
} // namespace Dmrg 

/*@}*/
#endif // GEOMETRY_DIR_H
