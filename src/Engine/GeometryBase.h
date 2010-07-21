// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

/*! \file GeometryBase.h
 *
 *  A class to isolate the geometry dependence of the Dmrg method.
 *  This is an abstract class and meant to simply provide a public interface
 */
#ifndef GEOMETRYBASE_HEADER_H
#define GEOMETRYBASE_HEADER_H

#include "Utils.h"
#include "GeometryImplementation.h"

namespace Dmrg {
	//! Interface to Geometry dependence of dmrg.
	//! Note that geometry proper is handled by the model connector's array (if any)
	//! To implement a new geometry inherit from this class and implement this interface
	template<typename RealType>
	class GeometryBase {
		public:
			typedef std::vector<size_t> BlockType;
			typedef GeometryImplementation<RealType> GeometryImplementationType;
			template<typename IoInputter>
			GeometryBase(IoInputter& io) : geometryImpl_(io)
			{
			}

			size_t connectionKind(size_t smax,size_t ind,size_t jnd) const
			{
				return geometryImpl_.connectionKind(smax,ind,jnd);
			}
			
			RealType operator()(size_t smax,size_t emin,
					  size_t ind,size_t edof1,size_t jnd,size_t edof2,size_t term) const
			{
				return geometryImpl_(smax,emin,ind,edof1,jnd,edof2,term);
			}
			
			RealType operator()(size_t ind,size_t edof1,size_t jnd,size_t edof2,size_t term) const
			{
				return geometryImpl_(ind,edof1,jnd,edof2,term);
			}
			
			size_t terms() const { return geometryImpl_.terms(); }
			
			void split(BlockType& S,std::vector<BlockType>& X,std::vector<BlockType>& Y,BlockType& E) const
			{
				geometryImpl_.split(S,X,Y,E);
			}
			
			bool connected(size_t i,size_t j) const
			{
				return geometryImpl_.connected(i,j);
			}
			
			size_t numberOfSites() const { return geometryImpl_.numberOfSites(); }
			
			
		private:
			GeometryImplementationType geometryImpl_;
	}; // class GeometryBase
} // namespace Dmrg 

/*@}*/
#endif
