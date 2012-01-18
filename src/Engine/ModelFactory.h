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

/*! \file ModelFactory.h
 *
 * From a string, chooses the model
 *
 */
 
#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include "ModelHubbard.h"
#include "ModelHeisenberg.h"

namespace Dmrg {
	
	template<typename ModelHelperType_,
	typename SparseMatrixType,
	typename GeometryType_,
	template<typename> class SharedMemoryTemplate>
	class ModelFactory  {

		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef VerySparseMatrix<SparseElementType> VerySparseMatrixType;

		// start models here:
		typedef ModelHubbard<ModelHelperType_,SparseMatrixType,GeometryType_,
				     SharedMemoryTemplate> ModelHubbardType;
		typedef ModelHeisenberg<ModelHelperType_,SparseMatrixType,GeometryType_,
					SharedMemoryTemplate> ModelHeisenbergType;
		// end models

		enum {HUBBARD_ONE_BAND,HEISENBERG_SPIN_ONEHALF};

	public:

		typedef std::vector<unsigned int long long> HilbertBasisType;
		typedef ModelHelperType_ ModelHelperType;
		typedef GeometryType_ GeometryType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename ModelHelperType::BlockType Block;
		typedef typename ModelHelperType::RealType RealType;

		typedef typename ModelHelperType::ConcurrencyType
				ConcurrencyType;
		typedef typename ModelHelperType::BasisType MyBasis;
		typedef typename ModelHelperType::BasisWithOperatorsType BasisWithOperatorsType;
//		typedef HamiltonianConnection<GeometryType,ModelHelperType,LinkProductType> HamiltonianConnectionType;
//		typedef SharedMemoryTemplate<HamiltonianConnectionType> SharedMemoryType;
//		typedef typename HamiltonianConnectionType::LinkProductStructType LinkProductStructType;
		typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename MyBasis::BasisDataType BasisDataType;

		template<typename SomeParametersType>
		ModelFactory(const SomeParametersType& params,PsimagLite::IoSimple::In& io,const GeometryType& geometry)
		: geometry_(geometry),
		  modelHubbard_(0),
		  modelHeisenberg_(0)
		{

			std::string name = params.model;
			if (name=="HubbardOneBand") {
				modelHubbard_ = new ModelHubbardType(io,geometry);
				ModelHubbardType::SharedMemoryType::setThreads(params.nthreads);
				model_=HUBBARD_ONE_BAND;
			} else if (name=="HeisenbergSpinOneHalf") {
				modelHeisenberg_ = new ModelHeisenbergType(io,geometry);
				ModelHeisenbergType::SharedMemoryType::setThreads(params.nthreads);
				model_=HEISENBERG_SPIN_ONEHALF;
			} else {
				std::string s(__FILE__);
				s += " Unknown model " + name + "\n";
				throw std::runtime_error(s.c_str());
			}
		}

		~ModelFactory()
		{
			switch (model_) {
			case HUBBARD_ONE_BAND:
				delete modelHubbard_;
				break;
			case HEISENBERG_SPIN_ONEHALF:
				delete modelHeisenberg_;
				break;
			}
		}

		size_t hilbertSize(size_t site) const
		{
			switch (model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->hilbertSize(site);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->hilbertSize(site);
			}
			assert(false);
			return 0;
		}

		void setNaturalBasis(std::vector<OperatorType> &creationMatrix,
				     SparseMatrixType &hamiltonian,
				     BasisDataType& q,
				     Block const &block,
				     RealType time) const
		{
			switch (model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
			}
		}

		void setOperatorMatrices(std::vector<OperatorType> &creationMatrix,
					 Block const &block)
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->setOperatorMatrices(creationMatrix,block);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->setOperatorMatrices(creationMatrix,block);
			}
		}

		PsimagLite::Matrix<SparseElementType> naturalOperator(const std::string& what,
								      size_t site,
								      size_t dof) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->naturalOperator(what,site,dof);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->naturalOperator(what,site,dof);
			}
			std::cerr<<__FILE__<<" Unknown model "<<model_<<"\n";
			throw std::runtime_error("naturalOperator\n");
		}

		void findElectrons(std::vector<size_t> &electrons,
				   const HilbertBasisType& basis,
				   size_t site) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->findElectrons(electrons,basis,site);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->findElectrons(electrons,basis,site);
			}
		}

		void setNaturalBasis(HilbertBasisType& basis,
				     std::vector<size_t>& q,
				     const std::vector<size_t>& block) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->setNaturalBasis(basis,q,block);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->setNaturalBasis(basis,q,block);
			}
		}

		void print(std::ostream& os) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->print(os);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->print(os);
			}
		}

		template<typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType& x,
					 const SomeVectorType& y,
					 ModelHelperType const &modelHelper) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->matrixVectorProduct(x,y,modelHelper);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->matrixVectorProduct(x,y,modelHelper);
			}
		}

		void addHamiltonianConnection(SparseMatrixType &matrix,const LeftRightSuperType& lrs) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->addHamiltonianConnection(matrix,lrs);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->addHamiltonianConnection(matrix,lrs);
			}
		}
		
		void hamiltonianConnectionProduct(std::vector<SparseElementType> &x,
						  std::vector<SparseElementType> const &y,
						  ModelHelperType const &modelHelper) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->hamiltonianConnectionProduct(x,y,modelHelper);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->hamiltonianConnectionProduct(x,y,modelHelper);
			}
		}

		void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->fullHamiltonian(matrix,modelHelper);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->fullHamiltonian(matrix,modelHelper);
			}
		}

		size_t maxConnections() const
		{
			return geometry_.maxConnections();
		}

	private:

		const GeometryType& geometry_;
		ModelHubbardType* modelHubbard_;
		ModelHeisenbergType* modelHeisenberg_;
		size_t model_;
	};     //class ModelFactory

	template<typename ModelHelperType,
		 typename SparseMatrixType,
		 typename GeometryType,
		 template<typename> class SharedMemoryTemplate>
	std::ostream& operator<<(std::ostream& os,const ModelFactory<ModelHelperType,
				 SparseMatrixType,GeometryType,SharedMemoryTemplate>& mf)
	{
		mf.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif
