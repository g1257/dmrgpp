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

/*! \file ExtendedHubbard1Orb.h
 *
 *  Hubbard + V_{ij} n_i n_j
 *
 */
#ifndef EXTENDED_HUBBARD_1ORB_H
#define EXTENDED_HUBBARD_1ORB_H
#include "ModelHubbard.h"
#include "LinkProdExtendedHubbard1Orb.h"

namespace Dmrg {
	//! Extended Hubbard for DMRG solver, uses ModelHubbard by containment
	template<typename ModelBaseType>
	class ExtendedHubbard1Orb : public ModelBaseType {

	public:

		typedef ModelHubbard<ModelBaseType> ModelHubbardType;
		typedef typename ModelBaseType::ModelHelperType ModelHelperType;
		typedef typename ModelBaseType::GeometryType GeometryType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef LinkProdExtendedHubbard1Orb<ModelHelperType> LinkProductType;
		typedef	typename ModelBaseType::MyBasis MyBasis;
		typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
		typedef typename ModelHelperType::BlockType Block;
		typedef typename ModelHubbardType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
		typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
		typedef typename ModelBaseType::InputValidatorType InputValidatorType;
		typedef typename ModelBaseType::SolverParamsType SolverParamsType;

		ExtendedHubbard1Orb(const SolverParamsType& solverParams,
		                    InputValidatorType& io,
		                    GeometryType const &dmrgGeometry)
		: ModelBaseType(solverParams,io,dmrgGeometry),
		  modelParameters_(io),
		  dmrgGeometry_(dmrgGeometry),
		  modelHubbard_(solverParams,io,dmrgGeometry)
		{}

		SizeType hilbertSize(SizeType site) const
		{
			return modelHubbard_.hilbertSize(site);
		}

		//! find creation operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
		                     SparseMatrixType &hamiltonian,
		                     BasisDataType &q,
				     Block const &block,
				     RealType time) const
		{

			modelHubbard_.setNaturalBasis(creationMatrix,hamiltonian,q,block,time);

			// add ni to creationMatrix
			setNi(creationMatrix,block);

			// add V_{ij} n_i n_j to hamiltonian
			addNiNj(hamiltonian,creationMatrix,block);
		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
		                         Block const &block) const
		{
			modelHubbard_.setOperatorMatrices(creationMatrix,block);
			// add ni to creationMatrix
			setNi(creationMatrix,block);
		}

		PsimagLite::Matrix<SparseElementType> naturalOperator(const PsimagLite::String& what,
								      SizeType site,
								      SizeType dof) const
		{
			Block block;
			block.resize(1);
			block[0]=site;
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="n") {
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[2].data);
				return tmp;
			} else {
				return modelHubbard_.naturalOperator(what,site,dof);
			}
		}

		//! find total number of electrons for each state in the basis
		void findElectrons(typename PsimagLite::Vector<SizeType> ::Type&electrons,
				   const typename PsimagLite::Vector<HilbertState>::Type& basis,
		                   SizeType site) const
		{
			modelHubbard_.findElectrons(electrons,basis,site);
		}

		//! find all states in the natural basis for a block of n sites
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(HilbertBasisType  &basis,
		                     typename PsimagLite::Vector<SizeType>::Type& q,
		                     const typename PsimagLite::Vector<SizeType>::Type& block) const
		{
			modelHubbard_.setNaturalBasis(basis,q,block);
		}

		void print(std::ostream& os) const
		{
			modelHubbard_.print(os);
		}

		//! Full hamiltonian from creation matrices cm
		void calcHamiltonian(SparseMatrixType &hmatrix,
		                     const typename PsimagLite::Vector<OperatorType>::Type& cm,
		                     Block const &block,
		                     RealType time,
		                     RealType factorForDiagonals=1.0)  const
		{
			hmatrix.makeDiagonal(cm[0].data.row());

			this->addConnectionsInNaturalBasis(hmatrix,cm,block);

			modelHubbard_.addDiagonalsInNaturalBasis(hmatrix,
			                                         cm,
			                                         block,
			                                         time,
			                                         factorForDiagonals);
		}

	private:

		ParametersModelHubbard<RealType>  modelParameters_;
		const GeometryType &dmrgGeometry_;
		ModelHubbardType modelHubbard_;

		//! Find n_i in the natural basis natBasis
		SparseMatrixType findOperatorMatrices(int i,
		                                      const typename PsimagLite::Vector<HilbertState>::Type& natBasis) const
		{

			SizeType n = natBasis.size();
			PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

			for (SizeType ii=0;ii<natBasis.size();ii++) {
				HilbertState ket=natBasis[ii];
				cm(ii,ii) = 0.0;
				for (SizeType sigma=0;sigma<2;sigma++)
					if (HilbertSpaceHubbardType::isNonZero(ket,i,sigma))
						cm(ii,ii) += 1.0;
			}
// 			std::cout<<cm;
// 			std::cout<<"******************************\n";
			SparseMatrixType creationMatrix(cm);
			return creationMatrix;
		}

		//! Full hamiltonian from creation matrices cm
		void addNiNj(SparseMatrixType &hmatrix,
		             const typename PsimagLite::Vector<OperatorType>::Type& cm,
		             Block const &block) const
		{
			//Assume block.size()==1 and then problem solved!! there are no connection if there's only one site ;-)
			assert(block.size()==1);
// 			for (SizeType sigma=0;sigma<DEGREES_OF_FREEDOM;sigma++) 
// 				for (SizeType sigma2=0;sigma2<DEGREES_OF_FREEDOM;sigma2++) 
// 					addNiNj(hmatrix,cm,block,sigma,sigma2);
		}

		void setNi(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
		           Block const &block) const
		{
			assert(block.size()==1);
			typename PsimagLite::Vector<HilbertState>::Type natBasis;
			typename PsimagLite::Vector<SizeType>::Type q;
			modelHubbard_.setNaturalBasis(natBasis,q,block);

			SparseMatrixType tmpMatrix = findOperatorMatrices(0,natBasis);
			RealType angularFactor= 1;
			typename OperatorType::Su2RelatedType su2related;
			su2related.offset = 1; //check FIXME
			OperatorType myOp(tmpMatrix,1,typename OperatorType::PairType(0,0),angularFactor,su2related);

			creationMatrix.push_back(myOp);
		}
	};	//class ExtendedHubbard1Orb

	template<typename ModelBaseType>
	std::ostream &operator<<(std::ostream &os,
	                         const ExtendedHubbard1Orb<ModelBaseType>& model)
	{
		model.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif // EXTENDED_HUBBARD_1ORB_H

