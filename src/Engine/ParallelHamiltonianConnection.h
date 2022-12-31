#ifndef PARALLELHAMILTONIANCONNECTION_H
#define PARALLELHAMILTONIANCONNECTION_H
#include "Concurrency.h"
#include "Vector.h"
#include "OpsForLink.hh"

namespace Dmrg {

template<typename HamiltonianConnectionType>
class ParallelHamiltonianConnection {

	typedef typename HamiltonianConnectionType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelHelperType::OperatorStorageType OperatorStorageType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename HamiltonianConnectionType::VectorType VectorType;
	typedef typename HamiltonianConnectionType::LinkType LinkType;
	typedef typename ModelHelperType::Aux AuxType;

public:

	ParallelHamiltonianConnection(VectorType& x,
	                              const VectorType& y,
	                              const HamiltonianConnectionType& hc,
	                              const AuxType& aux)
	    : x_(x),
	      y_(y),
	      hc_(hc),
	      aux_(aux),
	      xtemp_(ConcurrencyType::storageSize(ConcurrencyType::codeSectionParams.npthreads))
	{
		hc_.clearThreadSelves();
	}

	void doTask(SizeType taskNumber ,SizeType threadNum)
	{
		if (xtemp_[threadNum].size() != x_.size())
			xtemp_[threadNum].resize(x_.size(),0.0);

		if (taskNumber == 0) {
			hc_.modelHelper().hamiltonianLeftProduct(xtemp_[threadNum], y_, aux_);
//			const SparseMatrixType& hamiltonian = hc_.modelHelper().leftRightSuper().
//			        left().hamiltonian().getCRS();
//			hc_.kroneckerDumper().push(true, hamiltonian, y_);
			return;
		}

		if (taskNumber == 1) {
			hc_.modelHelper().hamiltonianRightProduct(xtemp_[threadNum], y_, aux_);
//			const SparseMatrixType& hamiltonian = hc_.modelHelper().leftRightSuper().
//			        right().hamiltonian().getCRS();
	//		hc_.kroneckerDumper().push(false, hamiltonian, y_);
			return;
		}

		assert(taskNumber > 1);
		taskNumber -= 2;

		OpsForLink<LeftRightSuperType> opsForLink = hc_.opsForLink();

		opsForLink.setPointer(taskNumber);

		hc_.modelHelper().fastOpProdInter(xtemp_[threadNum],
		                                  y_,
		                                  opsForLink.A().getCRS(),
		                                  opsForLink.B().getCRS(),
		                                  opsForLink.link(),
		                                  aux_);

//		hc_.kroneckerDumper().push(A->getCRS(),
//		                           B->getCRS(),
//		                           link2.value,
//		                           link2.fermionOrBoson,
//		                           y_);
	}

	SizeType tasks() const { return hc_.tasks() + 2; }

	void sync()
	{
		SizeType total = 0;
		for (SizeType threadNum = 0; threadNum < xtemp_.size(); threadNum++)
			if (xtemp_[threadNum].size() == x_.size()) total++;

		typename PsimagLite::Vector<ComplexOrRealType>::Type x(x_.size(),0);
		for (SizeType threadNum = 0; threadNum < total; threadNum++)
			for (SizeType i=0;i<x_.size();i++)
				x[i]+=xtemp_[threadNum][i];

		if (!ConcurrencyType::isMpiDisabled("HamiltonianConnection"))
			PsimagLite::MPI::allReduce(x);

		for (SizeType i=0;i<x_.size();i++)
			x_[i] += x[i];
	}

	template<typename SomeConcurrencyType,typename SomeOtherConcurrencyType>
	void sync(SomeConcurrencyType& conc,SomeOtherConcurrencyType& conc2)
	{
		conc.reduce(x_,conc2);
	}

private:

	VectorType& x_;
	const VectorType& y_;
	const HamiltonianConnectionType& hc_;
	const AuxType& aux_;
	typename PsimagLite::Vector<VectorType>::Type xtemp_;
};
}
#endif // PARALLELHAMILTONIANCONNECTION_H
