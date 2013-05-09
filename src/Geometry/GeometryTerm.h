// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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

/*! \file GeometryTerm.h
 *
 * Each geometry term represents a Hamiltonian connection term of the form
 * X_{ij} A_i B_j, where A and B are operators and X_{ij} are numbers.
 * Note that on-site Hamiltonian terms aren't included in the geometry since
 * they're trivially handled by the DMRG algorithm.
 */
#ifndef GEOMETRY_TERM_H
#define GEOMETRY_TERM_H

#include "GeometryDirection.h"
#include "GeometryFactory.h"
#include <cassert>
#include "String.h"

namespace PsimagLite {
	
	template<typename RealType>
	class GeometryTerm {
			typedef GeometryDirection<RealType,GeometryFactory> GeometryDirectionType;
		public:
			
			typedef typename GeometryFactory::AdditionalDataType AdditionalDataType;
			
			template<typename IoInputter>
			GeometryTerm(IoInputter& io,size_t termId,size_t linSize,bool debug=false) :
				linSize_(linSize),maxEdof_(0)
			{
				int x;
				io.readline(x,"DegreesOfFreedom=");
				if (x<=0) throw std::runtime_error("DegreesOfFreedom<=0 is an error\n");
				//std::cerr<<"DegreesOfFreedom "<<x<<"\n";
				size_t edof = (x==1) ? GeometryDirectionType::NUMBERS : GeometryDirectionType::MATRICES;
				String s;
				io.readline(s,"GeometryKind=");
				//std::cerr<<"GeometryKind "<<s<<"\n";

				String gOptions;
				io.readline(gOptions,"GeometryOptions=");
				//std::cerr<<"GeometryOptions "<<gOptions<<"\n";

				geometryFactory_.init(io,s,linSize);

				for (size_t i=0;i<geometryFactory_.dirs();i++) {
					directions_.push_back(GeometryDirectionType(io,i,edof,gOptions,geometryFactory_));
				}

				findMaxEdof();
				cacheValues();

				if (debug) {
					std::cerr<<"Cached values:\n";
					std::cerr<<cachedValues_;
					std::cerr<<"-----------\n";
				}
			}
			
			const RealType& operator()(size_t i1,size_t edof1,size_t i2,size_t edof2) const
			{
				int k1 = geometryFactory_.index(i1,edof1,maxEdof_);
				int k2 = geometryFactory_.index(i2,edof2,maxEdof_);
				assert(k1>=0 && k2>=0);
				return cachedValues_(k1,k2);
			}

			//assumes 1<smax+1 < emin
			const RealType& operator()(size_t smax,size_t emin,size_t i1,size_t edof1,size_t i2,size_t edof2) const
			{
				bool bothFringe = (geometryFactory_.fringe(i1,smax,emin) && geometryFactory_.fringe(i2,smax,emin));
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
					siteNew2 = geometryFactory_.getSubstituteSite(smax,emin,siteNew2);
				}
				return operator()(siteNew1,edofNew1,siteNew2,edofNew2);
			}
			
			bool connected(size_t smax,size_t emin,size_t i1,size_t i2) const
			{
				if (i1==i2) return false;

				bool bothFringe = (geometryFactory_.fringe(i1,smax,emin) && geometryFactory_.fringe(i2,smax,emin));

				if (!bothFringe) return geometryFactory_.connected(i1,i2);
				//std::cerr<<"fringe= "<<i1<<" "<<i2<<"\n";
				return true;
			}

			bool connected(size_t i1,size_t i2) const
			{
				return geometryFactory_.connected(i1,i2);
			}

			String label() const
			{
				return geometryFactory_.label();
			}

			size_t maxConnections() const
			{
				return geometryFactory_.maxConnections();
			}

			void fillAdditionalData(AdditionalDataType& additionalData,size_t ind,size_t jnd) const
			{
				geometryFactory_.fillAdditionalData(additionalData,ind,jnd);
			}

			size_t findReflection(size_t site) const
			{
				return geometryFactory_.findReflection(site);
			}

			size_t length(size_t i) const
			{
				return geometryFactory_.length(i);
			}

			size_t translate(size_t site,size_t dir,size_t amount) const
			{
				return geometryFactory_.translate(site,dir,amount);
			}

			void print(std::ostream& os,size_t linSize) const
			{
				size_t dofs = 1;
				for (size_t dof1=0;dof1<dofs;dof1++) {
					for (size_t dof2=0;dof2<dofs;dof2++) {
						os<<"dof1="<<dof1<<" dof2="<<dof2<<"\n";
						for (size_t i=0;i<linSize;i++) {
							for (size_t j=0;j<linSize;j++) {
								if (!connected(i,j)) {
									os<<0<<" ";
									continue;
								}
								os<<operator()(i,dof1,j,dof2)<<" ";
							}
							os<<"\n";
						}
						os<<"\n";
					}
				}
			}

			template<typename RealType_>	
			friend std::ostream& operator<<(std::ostream& os,const GeometryTerm<RealType_>& gt);
	
		private:

			void cacheValues()
			{
				size_t matrixRank = geometryFactory_.matrixRank(linSize_,maxEdof_);
				cachedValues_.resize(matrixRank,matrixRank);

				for (size_t i1=0;i1<linSize_;i1++) {
					for (size_t i2=0;i2<linSize_;i2++) {
						if (!geometryFactory_.connected(i1,i2)) continue;
						for (size_t edof1=0;edof1<maxEdof_;edof1++) {
							int k1 = geometryFactory_.index(i1,edof1,maxEdof_);
							if (k1<0) continue;
							for (size_t edof2=0;edof2<maxEdof_;edof2++) {
								int k2 = geometryFactory_.index(i2,edof2,maxEdof_);
								if (k2<0) continue;
								cachedValues_(k1,k2)=calcValue(i1,edof1,i2,edof2);
							}
						}
					}
				}
			}

			void findMaxEdof()
			{
				maxEdof_ = 0;
				for (size_t dir=0;dir<directions_.size();dir++) {
						maxEdof_ = directions_[dir].nRow();
						if (maxEdof_< directions_[dir].nCol())
							maxEdof_ = directions_[dir].nCol();
				}
			}

			RealType calcValue(size_t i1,size_t edof1,size_t i2,size_t edof2) const
			{
				if (!geometryFactory_.connected(i1,i2)) return 0.0;

				size_t dir = geometryFactory_.calcDir(i1,i2);
				assert(dir<directions_.size());
				return directions_[dir](i1,edof1,i2,edof2);
			}

			size_t linSize_;
			size_t maxEdof_;
			GeometryFactory geometryFactory_;
			typename Vector<GeometryDirectionType>::Type directions_;
			PsimagLite::Matrix<RealType> cachedValues_;
	}; // class GeometryTerm

	template<typename RealType>
	std::ostream& operator<<(std::ostream& os,const GeometryTerm<RealType>& gt)
	{
		os<<"#GeometryDirections="<<gt.directions_.size()<<"\n";
		for (size_t i=0;i<gt.directions_.size();i++) os<<gt.directions_[i];
		return os;
	}
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_TERM_H
