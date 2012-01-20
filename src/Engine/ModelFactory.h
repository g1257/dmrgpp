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
 * Note: Due to BaseModel needing the LinkProductType it is
 *       currently not possible to use virtual inheritance
 *       to dispatch the calls.
 *       The workaround in this class has been to do a manual
 *       dispatch, which is kind of error prone.
 *       To mitigate the problem, some things have been cached.
 *
 */
 
#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include "ModelHubbard.h"
#include "ModelHeisenberg.h"
#include "ExtendedHubbard1Orb.h"
#include "ModelFeBasedSc.h"
#include "FeAsBasedScExtended.h"
#include "Immm.h"

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
		typedef ExtendedHubbard1Orb<ModelHelperType_,SparseMatrixType,GeometryType_,
					    SharedMemoryTemplate> ModelHubbardExtType;
		typedef ModelFeBasedSc<ModelHelperType_,SparseMatrixType,GeometryType_,
				       SharedMemoryTemplate> FeBasedScType;
		typedef FeAsBasedScExtended<ModelHelperType_,SparseMatrixType,GeometryType_,
					    SharedMemoryTemplate> FeBasedScExtType;
		typedef Immm<ModelHelperType_,SparseMatrixType,GeometryType_,
			     SharedMemoryTemplate> ImmmType;
		// end models

		enum {HUBBARD_ONE_BAND,HEISENBERG_SPIN_ONEHALF,
		      HUBBARD_ONE_BAND_EXT,FEAS,FEAS_EXT,IMMM};

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
		: name_(params.model),
		  geometry_(geometry),
		  hilbertSize_(geometry.numberOfSites()),
		  q_(hilbertSize_.size()),
		  basis_(q_.size()),
		  modelHubbard_(0),
		  modelHeisenberg_(0),
		  modelHubbardExt_(0),
		  modelFeAs_(0),
		  modelFeAsExt_(0),
		  modelImmm_(0)
		{
			if (name_=="HubbardOneBand") {
				modelHubbard_ = new ModelHubbardType(io,geometry);
				ModelHubbardType::SharedMemoryType::setThreads(params.nthreads);
				model_=HUBBARD_ONE_BAND;
				init(modelHubbard_);
			} else if (name_=="HeisenbergSpinOneHalf") {
				modelHeisenberg_ = new ModelHeisenbergType(io,geometry);
				ModelHeisenbergType::SharedMemoryType::setThreads(params.nthreads);
				model_=HEISENBERG_SPIN_ONEHALF;
				init(modelHeisenberg_);
			} else if (name_=="HubbardOneBandExtended") {
				modelHubbardExt_ = new ModelHubbardExtType(io,geometry);
				ModelHubbardExtType::SharedMemoryType::setThreads(params.nthreads);
				model_=HUBBARD_ONE_BAND_EXT;
				init(modelHubbardExt_);
			} else  if (name_=="FeAsBasedSc") {
				modelFeAs_ = new FeBasedScType(io,geometry);
				FeBasedScType::SharedMemoryType::setThreads(params.nthreads);
				model_=FEAS;
				init(modelFeAs_);
			} else if (name_=="FeAsBasedScExtended") {
				modelFeAsExt_ = new FeBasedScExtType(io,geometry);
				FeBasedScExtType::SharedMemoryType::setThreads(params.nthreads);
				model_=FEAS_EXT;
				init(modelFeAsExt_);
			} else if (name_=="Immm") {
				modelImmm_ = new ImmmType(io,geometry);
				ImmmType::SharedMemoryType::setThreads(params.nthreads);
				model_=IMMM;
				init(modelImmm_);
			} else {
				std::string s(__FILE__);
				s += " Unknown model " + name_ + "\n";
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
			case HUBBARD_ONE_BAND_EXT:
				delete modelHubbardExt_;
				break;
			case FEAS:
				delete modelFeAs_;
				break;
			case FEAS_EXT:
				delete modelFeAsExt_;
				break;
			case IMMM:
				delete modelImmm_;
				break;
			}
		}

		void setNaturalBasis(std::vector<OperatorType> &creationMatrix,
				     SparseMatrixType &hamiltonian,
				     BasisDataType& q,
				     const std::vector<size_t>& block,
				     const RealType& time) const
		{
			switch (model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
			case FEAS:
				return modelFeAs_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
			case FEAS_EXT:
				return modelFeAsExt_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
			case IMMM:
				return modelImmm_->setNaturalBasis(creationMatrix,hamiltonian,q,block,time);
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
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->naturalOperator(what,site,dof);
			case FEAS:
				return modelFeAs_->naturalOperator(what,site,dof);
			case FEAS_EXT:
				return modelFeAsExt_->naturalOperator(what,site,dof);
			case IMMM:
				return modelImmm_->naturalOperator(what,site,dof);
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
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->findElectrons(electrons,basis,site);
			case FEAS:
				return modelFeAs_->findElectrons(electrons,basis,site);
			case FEAS_EXT:
				return modelFeAsExt_->findElectrons(electrons,basis,site);
			case IMMM:
				return modelImmm_->findElectrons(electrons,basis,site);
			}
		}

		void print(std::ostream& os) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->print(os);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->print(os);
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->print(os);
			case FEAS:
				return modelFeAs_->print(os);
			case FEAS_EXT:
				return modelFeAsExt_->print(os);
			case IMMM:
				return modelImmm_->print(os);
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
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->matrixVectorProduct(x,y,modelHelper);
			case FEAS:
				return modelFeAs_->matrixVectorProduct(x,y,modelHelper);
			case FEAS_EXT:
				return modelFeAsExt_->matrixVectorProduct(x,y,modelHelper);
			case IMMM:
				return modelImmm_->matrixVectorProduct(x,y,modelHelper);
			}
		}

		void addHamiltonianConnection(SparseMatrixType &matrix,const LeftRightSuperType& lrs) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->addHamiltonianConnection(matrix,lrs);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->addHamiltonianConnection(matrix,lrs);
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->addHamiltonianConnection(matrix,lrs);
			case FEAS:
				return modelFeAs_->addHamiltonianConnection(matrix,lrs);
			case FEAS_EXT:
				return modelFeAsExt_->addHamiltonianConnection(matrix,lrs);
			case IMMM:
				return modelImmm_->addHamiltonianConnection(matrix,lrs);
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
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->hamiltonianConnectionProduct(x,y,modelHelper);
			case FEAS:
				return modelFeAs_->hamiltonianConnectionProduct(x,y,modelHelper);
			case FEAS_EXT:
				return modelFeAsExt_->hamiltonianConnectionProduct(x,y,modelHelper);
			case IMMM:
				return modelImmm_->hamiltonianConnectionProduct(x,y,modelHelper);
			}
		}

		void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
		{
			switch(model_) {
			case HUBBARD_ONE_BAND:
				return modelHubbard_->fullHamiltonian(matrix,modelHelper);
			case HEISENBERG_SPIN_ONEHALF:
				return modelHeisenberg_->fullHamiltonian(matrix,modelHelper);
			case HUBBARD_ONE_BAND_EXT:
				return modelHubbardExt_->fullHamiltonian(matrix,modelHelper);
			case FEAS:
				return modelFeAs_->fullHamiltonian(matrix,modelHelper);
			case FEAS_EXT:
				return modelFeAsExt_->fullHamiltonian(matrix,modelHelper);
			case IMMM:
				return modelImmm_->fullHamiltonian(matrix,modelHelper);
			}
		}

		size_t hilbertSize(size_t site) const
		{
			return hilbertSize_[site];
		}

		void setOperatorMatrices(std::vector<OperatorType> &creationMatrix,
					 Block const &block)
		{
			assert(block.size()==1);
			return creationMatrix[block[0]];
		}

		void setNaturalBasis(HilbertBasisType& basis,
				     std::vector<size_t>& q,
				     const std::vector<size_t>& block) const
		{
			assert(block.size()==1);
			size_t index=block[0];
			basis=basis_[index];
			q=q_[index];
		}

		size_t maxConnections() const
		{
			return geometry_.maxConnections();
		}

		const std::string& name() const
		{
			return name_;
		}

	private:

		template<typename SomeModelType>
		void init(SomeModelType* model)
		{
			for (size_t i=0;i<hilbertSize_.size();i++) {
				std::vector<size_t> block(1,i);
				hilbertSize_[i] = model->hilbertSize(i);
				model->setNaturalBasis(basis_[i],q_[i],block);
			}
		}

		const std::string& name_;
		const GeometryType& geometry_;
		std::vector<size_t> hilbertSize_;
		std::vector<std::vector<size_t> > q_;
		std::vector<HilbertBasisType> basis_;
		// models start
		ModelHubbardType* modelHubbard_;
		ModelHeisenbergType* modelHeisenberg_;
		ModelHubbardExtType* modelHubbardExt_;
		FeBasedScType* modelFeAs_ ;
		FeBasedScExtType* modelFeAsExt_;
		ImmmType* modelImmm_;
		// models end
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
