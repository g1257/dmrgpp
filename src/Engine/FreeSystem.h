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

/*! \file FreeSystem.h
 *
 *  Calculates and returns static correlations and energies for a free system
 *  The system is specified by the model and geometry
 *
 */
#ifndef FREE_SYSTEM_H
#define FREE_SYSTEM_H

#include "Utils.h"

namespace Dmrg {
	// All interactions == 0
	// This is to test the DMRG code 
	template<typename FieldType,typename BlockType,typename ParametersModelType,typename GeometryType>
	class FreeSystem {
			
			static size_t const SPIN_UP=0,SPIN_DOWN=1;
		public:
			FreeSystem(const ParametersModelType& mp,size_t leg,size_t orbitals,size_t dof,bool verbose=false) : 
				mp_(mp),leg_(leg),orbitals_(orbitals),dof_(dof),verbose_(verbose),currentSize_(0),numberOfElectrons_(0) { }
			
			FieldType energy(size_t thisSize)
			{
				diagonalize(thisSize);
				FieldType sum=0;
				for (size_t i=0;i<numberOfElectrons_;i++) sum += eigenvalues_[i];
				return sum;
			}
			
			FieldType deltaCorrelation(size_t ind,size_t gamma,size_t jnd,size_t gamma2,size_t thisSize)
			{
				diagonalize(thisSize);
				return bAux(ind,gamma,jnd,gamma2)-aAux(ind,gamma,jnd,gamma2);
			}
			
			FieldType szCorrelation(size_t ind,size_t jnd,size_t thisSize)
			{
				return ninj(ind,0,jnd,0,thisSize)-ninj(ind,0,jnd,1,thisSize)
						-ninj(ind,1,jnd,0,thisSize)+ninj(ind,0,jnd,0,thisSize);
			}
			
			FieldType cCorrelation(size_t ind,size_t gamma,size_t sigma1,size_t jnd,size_t gamma2,size_t sigma2,size_t thisSize)
			{
				return cCorrelation(ind,gamma,sigma1,jnd,gamma2,sigma2,thisSize,numberOfElectrons_);
			}
			
			FieldType cCorrelationTotal(size_t ind,size_t gamma,size_t sigma1,size_t jnd,size_t gamma2,size_t sigma2,size_t thisSize)
			{
				return cCorrelation(ind,gamma,sigma1,jnd,gamma2,sigma2,thisSize,eigenvectors_.n_row());
			}
			
			FieldType ninj(size_t ind,size_t jnd,size_t thisSize)
			{
				FieldType tmp = 0;
				for (size_t spin=0;spin<2;spin++) 
					for (size_t spin2=0;spin2<2;spin2++) 
						tmp +=  ninj(ind,spin,jnd,spin2,thisSize);
				return tmp;
			}
			
			FieldType ninj(size_t ind,size_t spin,size_t jnd,size_t spin2,size_t thisSize)
			{
				FieldType sum=0;
				for (size_t gamma=0;gamma<orbitals_;gamma++) {
					for (size_t gamma2=0;gamma2<orbitals_;gamma2++) {
						sum += cCorrelation(ind,gamma,spin,ind,gamma,spin,thisSize)*
								cCorrelation(jnd,gamma2,spin2,jnd,gamma2,spin2,thisSize);
						if (spin==spin2 && gamma==gamma2) {
							sum -=  cCorrelation(ind,gamma,spin,jnd,gamma2,spin2,thisSize)*
								cCorrelation(jnd,gamma2,spin2,ind,gamma,spin,thisSize);
							if (ind==jnd) sum += cCorrelation(ind,gamma,spin,jnd,gamma2,spin2,thisSize)*
										cCorrelationTotal(jnd,gamma2,spin2,ind,gamma,spin,thisSize);
						}
						
					}
				}
				return sum;
			}
		private:
		
			void diagonalize(size_t thisSize)
			{
				if (thisSize==currentSize_) return;
				typedef typename GeometryType::ConnectorsType ConnectorsType;
				ConnectorsType connectors(dof_,thisSize);
			        std::vector<FieldType> defaultConnectors(2);
    				defaultConnectors[0] = mp_.hoppings(0,2); // x-direction
        			defaultConnectors[1] = mp_.hoppings(0,1); // y-direction
        			connectors.push(mp_.hoppings,defaultConnectors);

				//connectors.push(mp_.hoppingsOneSite,orbitals_,leg_);
      				
				GeometryType geometry_(connectors,1,leg_);
				eigenvectors_.resize(connectors.n_row()*dof_,connectors.n_row()*dof_);
				
				for (size_t i=0;i<connectors.n_row();i++) {
					for (size_t j=0;j<connectors.n_row();j++) {
						for (size_t s1=0;s1<dof_;s1++) { 
							for (size_t s2=0;s2<dof_;s2++) {  
								eigenvectors_(s1+i*dof_,s2+j*dof_)=  connectors.getMatrix(i,j);
								if (s1!=s2) eigenvectors_(s1+i*dof_,s2+j*dof_)=0;
							}
						}
					}
				}
				
				if (verbose_) {
					std::cerr<<"Matrix\n";
					std::cerr<<eigenvectors_;
				}
				utils::diag(eigenvectors_,eigenvalues_,'V');
				numberOfElectrons_ = size_t(thisSize * mp_.density);
				std::cerr<<"Diagonalizing matrix of rank="<<eigenvectors_.n_row()<<" electrons="<<numberOfElectrons_<<"\n";
				if (eigenvalues_.size()<numberOfElectrons_) throw std::runtime_error("freeSystem::diagonalize(): electrons > total levels\n");
				if (verbose_) {
					utils::vectorPrint(eigenvalues_,"eigenvalues",std::cerr);
					std::cerr<<"*************\n";
					std::cerr<<"Eigenvectors:\n";
					std::cerr<<eigenvectors_;
				}
				//utils::vectorPrint(eigenvalues_,"eigenvalues",std::cerr);
				//permutation_.resize(eigenvalues_.size());
				//utils::sort<double,double>(eigenvalues_,permutation_);
				currentSize_=thisSize;
				
			}
			
			FieldType aAux(size_t ind,size_t gamma,size_t jnd,size_t gamma2)
			{
				return nAux(ind,gamma,jnd,gamma2,SPIN_UP,SPIN_DOWN)* nAux(ind,gamma,jnd,gamma2,SPIN_DOWN,SPIN_UP);
			}
			
			FieldType bAux(size_t ind,size_t gamma,size_t jnd,size_t gamma2)
			{
				return nAux(ind,gamma,jnd,gamma2,SPIN_UP,SPIN_UP)* nAux(ind,gamma,jnd,gamma2,SPIN_DOWN,SPIN_DOWN);
			}
			
			FieldType nAux(size_t ind,size_t gamma,size_t jnd,size_t gamma2,size_t sigma1,size_t sigma2)
			{
				const psimag::Matrix<FieldType>& U=eigenvectors_;
				//size_t n = U.n_col();
				FieldType sum=0;
				for (size_t i=0;i<numberOfElectrons_;i++) {
					sum += conj(U(gamma+sigma1*orbitals_+ind*dof_,i))*U(gamma2+sigma2*orbitals_+jnd*dof_,i);
				}
				return sum;
			}
			
			FieldType cCorrelation(size_t ind,size_t gamma,size_t sigma1,size_t jnd,size_t gamma2,size_t sigma2,size_t thisSize,
					      size_t upto)
			{
				diagonalize(thisSize);
				const psimag::Matrix<FieldType>& U=eigenvectors_;
				FieldType sum=0;
				
				for (size_t i=0;i<upto;i++) {
					sum += conj(U(gamma+sigma1*orbitals_+ind*dof_,i))*U(gamma2+sigma2*orbitals_+jnd*dof_,i);
				}
				
				return sum;
			}
			
			
			psimag::Matrix<FieldType> eigenvectors_;
			std::vector<FieldType> eigenvalues_;
			ParametersModelType mp_;
			size_t leg_,orbitals_,dof_;
			bool verbose_;
			size_t currentSize_;
			size_t numberOfElectrons_;
	}; // FreeSystem
} // namespace Dmrg 

/*@}*/
#endif
