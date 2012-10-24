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
/** \ingroup DMRG */
/*@{*/

/*! \file FeBasedScExtedned.h
 *
 *  An implementation of a Hubbard model for Fe-based superconductors
 *  to use with the DmrgSolver
 *  This extends the FeAsBasedSc model to include JNN and JNNN couplings
 *
 */
#ifndef FEAS_BASED_SC_EX
#define FEAS_BASED_SC_EX
#include "ModelFeBasedSc.h"
#include "LinkProductFeAsExtended.h"

namespace Dmrg {
	template<
		typename ModelHelperType_,
		typename SparseMatrixType,
		typename GeometryType,
		template<typename> class SharedMemoryTemplate>
	class FeAsBasedScExtended : public ModelBase<
			ModelHelperType_,SparseMatrixType,GeometryType,
			LinkProductFeAsExtended<ModelHelperType_>,
			SharedMemoryTemplate> {
		
	public:
		typedef ModelFeBasedSc<ModelHelperType_,SparseMatrixType,
			GeometryType,SharedMemoryTemplate> ModelFeAsType;
		typedef typename ModelFeAsType::HilbertState HilbertState;
		typedef typename ModelFeAsType::HilbertBasisType HilbertBasisType;
		typedef ModelHelperType_ ModelHelperType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename OperatorType::Su2RelatedType Su2RelatedType;

		typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
		typedef LinkProductFeAsExtended<ModelHelperType> LinkProductType;
		typedef ModelBase<ModelHelperType,SparseMatrixType,GeometryType,
				LinkProductType,SharedMemoryTemplate> ModelBaseType;
		typedef	 typename ModelBaseType::MyBasis MyBasis;
		typedef	 typename ModelBaseType::BasisWithOperatorsType
				MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef typename MyBasis::BlockType BlockType;
		typedef typename ModelBaseType::InputValidatorType InputValidatorType;

		static const size_t NUMBER_OF_ORBITALS =
				ModelFeAsType::NUMBER_OF_ORBITALS;
		static const size_t SPIN_UP = ModelFeAsType::SPIN_UP;
		static const size_t SPIN_DOWN = ModelFeAsType::SPIN_DOWN;

		FeAsBasedScExtended(InputValidatorType& io,GeometryType const &geometry,ConcurrencyType& concurrency)
			: ModelBaseType(geometry,concurrency),modelParameters_(io), geometry_(geometry),
			  modelFeAs_(io,geometry,concurrency)
		{}

		size_t orbitals() const { return modelFeAs_.orbitals(); }

		size_t hilbertSize(size_t site) const { return modelFeAs_.hilbertSize(site); }

		void print(std::ostream& os) const { modelFeAs_.print(os); }

		//! find creation operator matrices for (i,sigma) in the natural basis,
		//! find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(std::vector<OperatorType> &creationMatrix,
				     SparseMatrixType &hamiltonian,
				     BasisDataType &q,
				     BlockType const &block,
				     const RealType& time)  const
		{
			blockIsSize1OrThrow(block);

			modelFeAs_.setNaturalBasis(creationMatrix,hamiltonian,q,block,time);

			// add S^+_i to creationMatrix
			setSplus(creationMatrix,block);

			// add S^z_i to creationMatrix
			setSz(creationMatrix,block);

			// add J_{ij} S^+_i S^-_j + S^-_i S^+_j to Hamiltonia
			addSplusSminus(hamiltonian,creationMatrix,block);

			// add J_{ij} S^z_i S^z_j to Hamiltonian
			addSzSz(hamiltonian,creationMatrix,block);

		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(
				std::vector<OperatorType> &creationMatrix,
				BlockType const &block) const
		{
			blockIsSize1OrThrow(block);

			modelFeAs_.setOperatorMatrices(creationMatrix,block);

			// add S^+_i to creationMatrix
			setSplus(creationMatrix,block);

			// add S^z_i to creationMatrix
			setSz(creationMatrix,block);
		}

		PsimagLite::Matrix<SparseElementType> naturalOperator(const std::string& what,
								      size_t site,
								      size_t dof) const
		{
			BlockType block;
			block.resize(1);
			block[0]=site;
			std::vector<OperatorType> creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="z") {
				PsimagLite::Matrix<SparseElementType> tmp;
				size_t x = 2*NUMBER_OF_ORBITALS+1;
				crsMatrixToFullMatrix(tmp,creationMatrix[x].data);
				return tmp;
			}

			if (what=="+") {
				PsimagLite::Matrix<SparseElementType> tmp;
				size_t x = 2*NUMBER_OF_ORBITALS;
				crsMatrixToFullMatrix(tmp,creationMatrix[x].data);
				return tmp;
			}

			if (what=="-") { // delta = c^\dagger * c^dagger
				PsimagLite::Matrix<SparseElementType> tmp;
				size_t x = 2*NUMBER_OF_ORBITALS;
				SparseMatrixType tmp2;
				transposeConjugate(tmp2,creationMatrix[x].data);
				crsMatrixToFullMatrix(tmp,tmp2);
				return tmp;
			}
			return modelFeAs_.naturalOperator(what,site,dof);
		}

		//! find all states in the natural basis for a block of n sites
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(std::vector<HilbertState>  &basis,
		                     std::vector<size_t>& q,
		                     const std::vector<size_t>& block) const
		{
			modelFeAs_.setNaturalBasis(basis,q,block);
		}
		
		void findElectrons(std::vector<size_t>& electrons,const std::vector<HilbertState>  &basis,size_t site) const
		{
			modelFeAs_.findElectrons(electrons,basis,site);
		}

	private:

		// add S^+_i to creationMatrix
		void setSplus(
				std::vector<OperatorType> &creationMatrix,
				const BlockType& block) const
		{
			SparseMatrixType m;
			cDaggerC(m,creationMatrix,block,1.0,SPIN_UP,SPIN_DOWN);
			Su2RelatedType su2related;
			size_t offset = 2*NUMBER_OF_ORBITALS;
			su2related.source.push_back(offset);
			su2related.source.push_back(offset+1);
			su2related.source.push_back(offset);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(1);
			su2related.offset = 1;

			OperatorType sPlus(m,1,typename OperatorType::PairType(2,2),-1,
					su2related);
			creationMatrix.push_back(sPlus);
		}

		// add S^z_i to creationMatrix
		void setSz(
				std::vector<OperatorType> &creationMatrix,
				const BlockType& block) const
		{
			SparseMatrixType m1,m2;
			cDaggerC(m1,creationMatrix,block,0.5,SPIN_UP,SPIN_UP);
			cDaggerC(m2,creationMatrix,block,-0.5,SPIN_DOWN,SPIN_DOWN);
			Su2RelatedType su2related2;
			SparseMatrixType m = m1;
		       	m += m2;
			OperatorType sz(m,1,typename OperatorType::PairType(2,1),
					1.0/sqrt(2.0),su2related2);
			creationMatrix.push_back(sz);
		}

		// add S^+_i to creationMatrix
		void cDaggerC(
				SparseMatrixType& sum,
				const std::vector<OperatorType> &creationMatrix,
				const BlockType& block,
				RealType value,
				size_t spin1,
				size_t spin2) const
		{
			SparseMatrixType tmpMatrix,tmpMatrix2;
			for (size_t orbital=0;orbital<NUMBER_OF_ORBITALS;orbital++) {
				transposeConjugate(tmpMatrix2,
						creationMatrix[orbital+spin2*NUMBER_OF_ORBITALS].data);
				multiply(tmpMatrix,
						creationMatrix[orbital+spin1*NUMBER_OF_ORBITALS].data,
						tmpMatrix2);
				multiplyScalar(tmpMatrix2,tmpMatrix,value);
				if (orbital == 0) sum = tmpMatrix2;
				else sum += tmpMatrix2;
			}
		}

		// add J_{ij} S^+_i S^-_j + S^-_i S^+_j to Hamiltonia
		void addSplusSminus(
				SparseMatrixType &hamiltonian,
				const std::vector<OperatorType> &creationMatrix,
				const BlockType& block) const
		{
			// nothing if block.size == 1
		}

		// add J_{ij} S^z_i S^z_j to Hamiltonian
		void addSzSz(
				SparseMatrixType &hamiltonian,
				const std::vector<OperatorType> &creationMatrix,
				const BlockType& block) const
		{
			// nothing if block.size == 1
		}

		void blockIsSize1OrThrow(const BlockType& block) const
		{
			if (block.size()==1) return;
			throw std::runtime_error(
				"FeAsBasedExtended:: Added blocks must be of size 1"
					"or is unimplemented otherwise\n");
		}

		ParametersModelFeAs<RealType>  modelParameters_;
		GeometryType const &geometry_;
		ModelFeAsType modelFeAs_;
	};     //class FeAsBasedScExtended

	template<
		typename ModelHelperType,
		typename SparseMatrixType,
		typename GeometryType,
		template<typename> class SharedMemoryTemplate
		>
	std::ostream &operator<<(std::ostream &os,const FeAsBasedScExtended<
		ModelHelperType,
		SparseMatrixType,
		GeometryType,
		SharedMemoryTemplate
		>& model)
	{
		model.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif // FEAS_BASED_SC_EX
