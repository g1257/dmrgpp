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

/*! \file GeometryTerm.h
 *
 * Each geometry term represents a Hamiltonian connection term of the form
 * X_{ij} A_i B_j, where A and B are operators and X_{ij} are numbers.
 * Note that on-site Hamiltonian terms aren't included in the geometry since
 * they're trivially handled by the DMRG algorithm.
 */
#ifndef GEOMETRY_TERM_H
#define GEOMETRY_TERM_H

#include "Utils.h"
#include "GeometryDirection.h"
#include "GeometryBase.h"

namespace Dmrg {
	
	template<typename RealType>
	class GeometryTerm {
			typedef GeometryDirection<RealType,GeometryBase> GeometryDirectionType;
		public:
			
			template<typename IoInputter>
			GeometryTerm(IoInputter& io,size_t termId,size_t linSize,bool debug=false) :
				linSize_(linSize)
			{
				int x;
				io.readline(x,"DegreesOfFreedom=");
				if (x<=0) throw std::runtime_error("DegreesOfFreedom<=0 is an error\n");
				//std::cerr<<"DegreesOfFreedom "<<x<<"\n";
				edof_ = x;
				std::string s;
				io.readline(s,"GeometryKind=");
				//std::cerr<<"GeometryKind "<<s<<"\n";

				std::string gOptions;
				io.readline(gOptions,"GeometryOptions=");
				//std::cerr<<"GeometryOptions "<<gOptions<<"\n";

				geometryBase_.init(io,s,linSize);

				for (size_t i=0;i<geometryBase_.dirs();i++) {
					directions_.push_back(GeometryDirectionType(io,i,edof_,gOptions,geometryBase_));
				}
				
				cachedValues_.resize(linSize*linSize*edof_*edof_);

				for (size_t i1=0;i1<linSize;i1++)
					for (size_t i2=0;i2<linSize;i2++)
						for (size_t edof1=0;edof1<edof_;edof1++)
							for (size_t edof2=0;edof2<edof_;edof2++)
								cachedValues_[pack(i1,edof1,i2,edof2)]=
									calcValue(i1,edof1,i2,edof2);
				if (debug) {
					std::cerr<<"Cached values:\n";
					std::cerr<<cachedValues_;
					std::cerr<<"-----------\n";
				}
			}
			
			const RealType& operator()(size_t i1,size_t edof1,size_t i2,size_t edof2) const
			{
				size_t p = pack(i1,edof1,i2,edof2);
				return cachedValues_[p];
			}
			
			//assumes 1<smax+1 < emin
			const RealType& operator()(size_t smax,size_t emin,
				size_t i1,size_t edof1,size_t i2,size_t edof2) const
			{

				bool bothFringe = (geometryBase_.fringe(i1,smax,emin) && geometryBase_.fringe(i2,smax,emin));
				size_t siteNew1 = i1;
				size_t siteNew2 = i2;
				size_t edofNew1 = edof1;
				size_t edofNew2 = edof2;
				if (bothFringe) {
					if (i2<i1) {
						siteNew1 = i2;
						siteNew2 = i1;
						edofNew1 = edof2;
						edofNew2 = edof1;
					}
					siteNew2 = geometryBase_.getSubstituteSite(smax,emin,siteNew2);
				}
				
				size_t p = pack(siteNew1,edofNew1,siteNew2,edofNew2);
				return cachedValues_[p];
			}
			
			bool connected(size_t smax,size_t emin,size_t i1,size_t i2) const
			{
				if (i1==i2) return false;

				bool bothFringe = (geometryBase_.fringe(i1,smax,emin) && geometryBase_.fringe(i2,smax,emin));

				if (!bothFringe) return geometryBase_.connected(i1,i2);
				//std::cerr<<"fringe= "<<i1<<" "<<i2<<"\n";
				return true;
			}

			bool connected(size_t i1,size_t i2) const
			{
				return geometryBase_.connected(i1,i2);
			}

			std::string label() const
			{
				return geometryBase_.label();
			}

			template<typename RealType_>	
			friend std::ostream& operator<<(std::ostream& os,const GeometryTerm<RealType_>& gt);
	
		private:	
			
			RealType calcValue(size_t i1,size_t edof1,size_t i2,size_t edof2) const
			{
				if (!geometryBase_.connected(i1,i2)) return 0.0;

				size_t dir = geometryBase_.calcDir(i1,i2);
				if (directions_[dir].constantValues()) {
					return directions_[dir](edof1,edof2);
				}
				
				return directions_[dir](i1,edof1,i2,edof2);
			}
			
			size_t pack(size_t i1,size_t edof1,size_t i2,size_t edof2) const
			{
				return edof1+i1*edof_+(edof2+i2*edof_)*linSize_*edof_;
			}

			size_t linSize_;
			size_t edof_;
			GeometryBase geometryBase_;
			std::vector<GeometryDirectionType> directions_;
			std::vector<RealType> cachedValues_;
	}; // class GeometryTerm

	template<typename RealType>
	std::ostream& operator<<(std::ostream& os,const GeometryTerm<RealType>& gt)
	{
		os<<"#GeometryDirections="<<gt.directions_.size()<<"\n";
		for (size_t i=0;i<gt.directions_.size();i++) os<<gt.directions_[i];
		return os;
	}
} // namespace Dmrg

/*@}*/
#endif // GEOMETRY_TERM_H
