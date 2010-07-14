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

/*! \file Connectors.h
 *
 *  FIXME
 *
 */
#ifndef CONNECTORS_H
#define CONNECTORS_H

#include "Connector.h"

namespace Dmrg {
	
	template<typename FieldType>
	class Connectors {
	public:
		
		typedef Connector<FieldType> ConnectorType;
		
		Connectors(size_t dof,size_t linSize) : dof_(dof), linSize_(linSize) 
		{
		}
		
		void push(const psimag::Matrix<FieldType>& connectorMatrix,FieldType defaultConnector) 
			
		{
			ConnectorType c(connectorMatrix,defaultConnector);
			connectors_.push_back(c);
		}
		
		void push(const psimag::Matrix<FieldType>& connectorMatrix,const std::vector<FieldType>& defaultConnector) 
			
		{
			ConnectorType c(connectorMatrix,defaultConnector);
			connectors_.push_back(c);
		}

		void push(	const std::vector<psimag::Matrix<FieldType> >& connectorsOneSite,
	  			size_t numberOfOrbitals,
				size_t leg) 
		{
		
			ConnectorType c(connectorsOneSite,linSize_,numberOfOrbitals,dof_,leg);
			connectors_.push_back(c);
		}
		
		size_t linSize() const { return linSize_; }
		
		size_t dof() const { return dof_; }
		
		FieldType getMatrix(size_t i,size_t j,size_t what = 0) const
		{
			return connectors_[what].getMatrix(i,j);
		}
		
		FieldType operator()(size_t i,size_t j,size_t what = 0) const
		{
			return connectors_[what](i,j);
		}
		
		FieldType defaultValue(size_t dir = 0, size_t a = 0, size_t b = 0,size_t what = 0) const
		{
			return connectors_[what].defaultValue(dir,a,b);
		}
		
		size_t n_row(size_t what = 0) const { return connectors_[what].n_row(); }
		
		size_t size() const { return connectors_.size(); }
		
		
	private:
		size_t dof_,linSize_;
		std::vector<ConnectorType> connectors_;
			
	}; //Connectors
	
// 	template<typename FieldType>
// 	std::ostream& operator<<(std::ostream& os,Connectors<FieldType>& connectors)
// 	{
// 		os<<connectors.connectors_;
// 		return os;
// 	}
} // namespace Dmrg
/*@}*/
#endif //TJ_ONEORBITAL_HEADER_H
