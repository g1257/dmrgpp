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

/*! \file TargetParamsCommon.h
 *
 *  FIXME
 */
#ifndef TARGET_PARAMS_COMMON_H
#define TARGET_PARAMS_COMMON_H

#include "Utils.h"

namespace Dmrg {
	//! Coordinates reading of TargetSTructure from input file
	template<typename ModelType>
	class TargetParamsCommon {
		public:
			typedef typename ModelType::RealType RealType;
			
			typedef typename ModelType::OperatorType OperatorType;
			typedef typename OperatorType::PairType PairType;
			typedef typename OperatorType::SparseMatrixType SparseMatrixType;
			typedef typename SparseMatrixType::value_type ComplexOrReal;
			typedef PsimagLite::Matrix<ComplexOrReal> MatrixType;

			template<typename IoInputter>
			TargetParamsCommon(IoInputter& io,const ModelType& model)
				: filename("tst.txt"),sites(0),startingLoops(0),model_(model)
			{

				io.readline(filename,"TSPFilename="); // filename
				io.read(sites,"TSPSites");
				io.read(startingLoops,"TSPLoops");
			
				data_.resize(sites.size());
				aOperators.resize(sites.size());
				typename ModelType::HilbertBasisType basis;
				model_.setNaturalBasis(basis,1);
				model_.findElectrons(electrons,basis);
			
				for (size_t i=0;i<sites.size();i++) {
					std::string s;
					io.readline(s,"TSPOperator=");
					if (s == "cooked") {
						io.readline(s,"COOKED_OPERATOR=");
						std::vector<size_t> v;
						io.read(v,"COOKED_EXTRA");
						setCookedData(i,s,v);
					} else {
						PsimagLite::Matrix<ComplexOrReal> m;
						io.readMatrix(m,"RAW_MATRIX");
						setRawData(i,m);
					}
					int fermiSign=0;
					io.readline(fermiSign,"FERMIONSIGN=");
					std::pair<size_t,size_t> jmValues;
					std::vector<size_t> v(2);
					io.readKnownSize(v,"JMVALUES");
					jmValues.first = v[0]; jmValues.second = v[1];
					RealType angularFactor;
					io.readline(angularFactor,"AngularFactor=");
					//tsp.set(i,fermiSign,jmValues,angularFactor);
					SparseMatrixType data(data_[i]);
	
					// FIXME: su2related needs to be set properly for when SU(2) is running: 
					typename OperatorType::Su2RelatedType su2Related; 
					OperatorType myOp(data,fermiSign, jmValues,angularFactor,su2Related);
					aOperators[i] = myOp;
				}
			}
			
			// I know, there is public data here FIXME!!
			std::string filename;
			std::vector<size_t> sites;
			std::vector<size_t> startingLoops;
			std::vector<OperatorType> aOperators;
			std::vector<size_t> electrons;
		
		private:
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
			}
			
			const ModelType& model_;
			std::vector<MatrixType> data_; 
	}; // class TargetParamsCommon
	
	template<typename ModelType>
	inline std::ostream&
	operator<<(std::ostream& os,const TargetParamsCommon<ModelType>& t)
	{
		os<<"#TargetParams.operators="<<t.aOperators.size()<<"\n";
		for (size_t i=0;i<t.aOperators.size();i++) {
			os<<"#TargetParams.operator "<<i<<"\n";
			os<<t.aOperators[i];
		}
		os<<"#TargetParams.electrons\n";
		os<<t.electrons;
		os<<"#TargetParams.site="<<t.sites;
		os<<"#TargetParams.startingLoop="<<t.startingLoops<<"\n";
		os<<"#TargetParams.filename="<<t.filename<<"\n";
		return os;
	}
} // namespace Dmrg 

/*@}*/
#endif // TARGET_PARAMS_COMMON_H
