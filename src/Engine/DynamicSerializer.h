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

/*! \file DynamicSerializer.h
 *
 *  Serialize dynamic DMRG data
 */
#ifndef DYN_SERIAL_H
#define DYN_SERIAL_H

#include "Utils.h"
#include "IoSimple.h"

namespace Dmrg {
	
	template<
		typename RealType,
		typename VectorType,
		typename ContinuedFractionType>
	class DynamicSerializer {
	public:

		DynamicSerializer(
				const ContinuedFractionType& cf,
				size_t site,
				const std::vector<VectorType>& targetVectors)
		: cf_(cf),
		  site_(site),
		  targetVectors_(targetVectors)
		{}

		template<typename IoInputType>
		DynamicSerializer(IoInputType& io,size_t lastInstance = 0)
		: cf_(io)
		{
			std::string s = "#DCENTRALSITE=";
			int xi=0;
			io.readline(xi,s);
			if (xi<0) throw std::runtime_error(
					"DynamicSerializer:: site cannot be negative\n");
			site_ = xi;

			s = "#DNUMBEROFVECTORS=";
			io.readline(xi,s);
			if (xi<=0) throw std::runtime_error(
					"DynamicSerializer:: n. of vectors must be positive\n");
			targetVectors_.resize(xi);
			for (size_t i=0;i<targetVectors_.size();i++) {
				s = "targetVector"+utils::ttos(i);
				targetVectors_[i].load(io,s);
			}
		}

		template<typename IoOutputter>
		void save(IoOutputter& io) const
		{
			cf_.save(io);

			std::string s = "#DCENTRALSITE=" + utils::ttos(site_);
			io.printline(s);
			s = "#DNUMBEROFVECTORS="+utils::ttos(targetVectors_.size());
			io.printline(s);
			for (size_t i=0;i<targetVectors_.size();i++) {
				std::string label = "targetVector"+utils::ttos(i);
				targetVectors_[i].save(io,label);
			}
		}

		size_t size(size_t i=0) const
		{
			return  targetVectors_[i].size();
		}

		size_t site() const
		{
			return  site_;
		}

		const VectorType& vector(size_t i=0) const
		{
			return targetVectors_[i];
		}

	private:
		const ContinuedFractionType& cf_;
		size_t site_;
		std::vector<VectorType> targetVectors_;
	}; // class TimeSerializer
} // namespace Dmrg 

/*@}*/
#endif
