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

/*! \file TargetStructureParams.h
 *
 *  A class to isolate the geometry dependence of the Dmrg method.
 *  This is an abstract class and meant to simply provide a public interface
 */
#ifndef TARGET_STRUCT_PARAMS_H
#define TARGET_STRUCT_PARAMS_H

#include "Utils.h"

namespace Dmrg {
	//! Coordinates reading of TargetSTructure from input file
	template<typename TargettingStructureType,typename ModelType>
	class TargetStructureParams {
		public:
		typedef typename ModelType::RealType RealType;
		
		typedef typename ModelType::OperatorType OperatorType;
		typedef typename OperatorType::PairType PairType;
		typedef typename OperatorType::SparseMatrixType SparseMatrixType;
		typedef typename SparseMatrixType::value_type ComplexOrReal;
		typedef psimag::Matrix<ComplexOrReal> MatrixType;
		
		TargetStructureParams(TargettingStructureType& targetStruct,const ModelType& model)
			: targetStruct_(targetStruct),model_(model)
		{
		}
		
		void init(
				const std::string& filename,
    				RealType timeStep,
    				const std::vector<size_t>& sites,
			  	const std::vector<size_t>& loops)
		{
			targetStruct_.filename = filename;
			targetStruct_.timeStep = timeStep;
			targetStruct_.sites=sites;
			targetStruct_.startingLoops=loops;
			data_.resize(sites.size());
			targetStruct_.aOperators.resize(sites.size());
			typename ModelType::HilbertBasisType basis;
			model_.setNaturalBasis(basis,1);
			model_.findElectrons(targetStruct_.electrons,basis);
			
		}
		
		void setCookedData(size_t i,const std::string& s,const std::vector<size_t>& v)
		{
			data_[i]=model_.getOperator(s,v[0],v[1]);
		}
		
		void setRawData(size_t i,const MatrixType& m)
		{
			data_[i]=m;
		}
		
		void set(size_t i,int fermiSign,const PairType& jmValues,RealType angularFactor)
		{
			SparseMatrixType data(data_[i]);

			// FIXME: su2related needs to be set properly for when SU(2) is running: 
			typename OperatorType::Su2RelatedType su2Related; 
			OperatorType myOp(data,fermiSign, jmValues,angularFactor,su2Related);
			targetStruct_.aOperators[i] = myOp;
		}
		
		private:
		TargettingStructureType& targetStruct_;
		const ModelType& model_;
		std::vector<MatrixType> data_; 
	}; // class TargetStructureParams
} // namespace Dmrg 

/*@}*/
#endif
