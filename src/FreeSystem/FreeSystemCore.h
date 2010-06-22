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

/*! \file FreeSystemCore.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef FREE_SYSTEM_CORE_H
#define FREE_SYSTEM_CORE_H

#include "Utils.h"
#include "HilbertVector.h"
#include "FreeOperator.h"

namespace Dmrg {
	// All interactions == 0
	template<typename FieldType>
	class FreeSystemCore {
			
			typedef psimag::Matrix<FieldType> MatrixType;
	
		public:
			
			typedef HilbertVector<FieldType> HilbertVectorType;
			typedef FreeOperator<MatrixType,HilbertVectorType> FreeOperatorType;
			typedef typename HilbertVectorType::HilbertTermType HilbertTermType;
			
			FreeSystemCore(const MatrixType& t,size_t dof,bool verbose=false) :
				t_(t),dof_(dof),verbose_(verbose)
			{
				diagonalize();
				std::cerr<<"Created core "<<eigenvectors_.n_row()<<"  times "<<eigenvectors_.n_col()<<"\n";
			}
			
			HilbertVectorType newState()
			{
				HilbertVectorType tmp(t_.n_row(),dof_);
				return tmp;	
			}
			
			HilbertVectorType newGroundState(const std::vector<size_t>& ne)
			{
				HilbertVectorType tmp(t_.n_row(),dof_);
				tmp.fill(ne);
				return tmp;	
			}
			
			FreeOperatorType newSimpleOperator(const std::string& label,size_t site,size_t flavor)
			{
				FreeOperatorType tmp(eigenvectors_,label,site,flavor,dof_);
				return tmp;
			}
			
		private:
		
			void diagonalize()
			{
				eigenvectors_ = t_;
				
				if (verbose_) {
					std::cerr<<"Matrix\n";
					std::cerr<<eigenvectors_;
				}
				utils::diag(eigenvectors_,eigenvalues_,'V');
				
				if (verbose_) {
					utils::vectorPrint(eigenvalues_,"eigenvalues",std::cerr);
					std::cerr<<"*************\n";
					std::cerr<<"Eigenvectors:\n";
					std::cerr<<eigenvectors_;
				}
			}

			const MatrixType& t_;
			size_t dof_;
			bool verbose_;
			psimag::Matrix<FieldType> eigenvectors_;
			std::vector<FieldType> eigenvalues_;
	}; // FreeSystemCore
} // namespace Dmrg 

/*@}*/
#endif
