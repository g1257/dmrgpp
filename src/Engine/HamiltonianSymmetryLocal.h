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

/*! \file HamiltonianSymmetryLocal.h
 *
 *  This class contains the implementation of local symmetries
 *  An object of this class is contained by DmrgBasisImplementation
 *
 */
#ifndef HAM_SYMM_LOCAL_H
#define HAM_SYMM_LOCAL_H
#include "Sort.h" // in PsimagLite
#include "BasisData.h"

namespace Dmrg {
	template<typename RealType,typename SparseMatrixType>
	class	HamiltonianSymmetryLocal {
		typedef  BasisData<std::pair<size_t,size_t> > BasisDataType;
		typedef PsimagLite::CrsMatrix<RealType> FactorsType;
		
		public:
			static int const MAX = 100;

			static size_t encodeQuantumNumber(const std::vector<size_t>& v)
			{
				size_t x= v[0] + v[1]*MAX;
				if (v.size()==3) x += v[2]*MAX*MAX;
				return x;
			}

			static std::vector<size_t> decodeQuantumNumber(size_t q)
			{
				std::vector<size_t> v(2);
				size_t tmp = q ;
				v[1] = size_t(tmp/MAX);
				v[0] = tmp % MAX;
				return v;
			}

			static size_t pseudoQuantumNumber(const std::vector<size_t>& targets)
			{
				return encodeQuantumNumber(targets);
			}

			size_t getFlavor(size_t i) const
			{
				return 0; // meaningless
			}

			//! find quantum numbers for each state of this basis, 
			//! considered symmetries for this model are: n_up and n_down
			static  void findQuantumNumbers(std::vector<size_t> &q,const BasisDataType& basisData) 
			{
				q.clear();
				std::vector<size_t> qn(2);
				for (size_t i=0;i<basisData.electronsUp.size();i++) {
					// nup
					qn[0] = basisData.electronsUp[i];
					// ndown
					qn[1] = basisData.electronsDown[i];
					// total spin
					q.push_back(encodeQuantumNumber(qn));
				}
			}

			template<typename SolverParametersType>
			void calcRemovedIndices(std::vector<size_t>& removedIndices,
						std::vector<RealType>& eigs,
						size_t kept,
						const SolverParametersType& solverParams) const
			{
				if (eigs.size()<=kept) return;
				// we sort the eigenvalues
				// note: eigenvalues are not ordered because DensityMatrix is diagonalized in blocks
				std::vector<size_t> perm(eigs.size());
				Sort<std::vector<RealType> > sort;
				sort.sort(eigs,perm);
				std::vector<size_t> permInverse(perm.size());
				for (size_t i=0;i<permInverse.size();i++) permInverse[perm[i]]=i;
				
				size_t target = eigs.size()-kept;
				
				removedIndices.clear();
				for (size_t i=0;i<target;i++) {
					if (removedIndices.size()>=target) break;
					if (PsimagLite::isInVector(removedIndices,perm[i])>=0) continue;
					removedIndices.push_back(perm[i]);
					
					
				}
				//std::cerr<<"target="<<target<<" size="<<eigs.size();
				//std::cerr<<" kept="<<kept<<" achieved="<<removedIndices.size()<<"\n";
			}

			const FactorsType& getFactors() const 
			{
				return factors_;
			}

			void createDummyFactors(size_t ns, size_t ne)
			{
				factors_.makeDiagonal(ns*ne,1);
			}

			template<typename IoInputter>
			void load(IoInputter& io) 
			{
				size_t tmp=0;
				io.readline(tmp,"#FACTORSSIZE=");
				factors_.makeDiagonal(tmp,1);
			}

			template<typename IoOutputter>
			void save(IoOutputter& io) const
			{
				// don't print factors since they're the identity anywaysfactors_
//				io.print("#FACTORSSIZE=",factors_.row());
				std::string tmp = ttos(factors_.row());
				std::string s="#FACTORSSIZE="+tmp;
				io.printline(s);
			}

		private:
			FactorsType factors_;
	}; //class HamiltonianSymmetryLocal
} // namespace Dmrg

/*@}*/
#endif
