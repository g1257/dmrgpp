#ifndef PARALLELHAMILTONIANCONNECTION_H
#define PARALLELHAMILTONIANCONNECTION_H
#include "Concurrency.h"
#include "Vector.h"

namespace Dmrg {

template<typename HamiltonianConnectionType>
class ParallelHamiltonianConnection {

	typedef typename HamiltonianConnectionType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename HamiltonianConnectionType::GeometryType GeometryType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename HamiltonianConnectionType::VectorType VectorType;
	typedef typename HamiltonianConnectionType::LinkType LinkType;

public:

	ParallelHamiltonianConnection(VectorType& x,
	                              const VectorType& y,
	                              const HamiltonianConnectionType& hc)
	    : x_(x),
	      y_(y),
	      hc_(hc),
	      xtemp_(ConcurrencyType::storageSize(ConcurrencyType::codeSectionParams.npthreads))
	{}

	      void doTask(SizeType taskNumber ,SizeType threadNum)
	{
		if (xtemp_[threadNum].size() != x_.size())
			xtemp_[threadNum].resize(x_.size(),0.0);

		ComplexOrRealType tmp = 0.0;

		if (taskNumber == 0) {
			hc_.modelHelper().hamiltonianLeftProduct(xtemp_[threadNum],y_);
			return;
		}

		if (taskNumber == 1) {
			hc_.modelHelper().hamiltonianRightProduct(xtemp_[threadNum],y_);
			return;
		}

		assert(taskNumber > 1);
		taskNumber -= 2;
		AdditionalDataType additionalData;
		SizeType xx = 0;
		ProgramGlobals::ConnectionEnum type;
		SizeType term = 0;
		SizeType dofs =0;
		hc_.prepare(xx,type,tmp,term,dofs,additionalData,taskNumber);

		linkProduct(xtemp_[threadNum],y_,xx,type,tmp,term,dofs,additionalData);
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

	//! Computes x+=H_{ij}y where H_{ij} is a Hamiltonian that connects system and environment
	void linkProduct(typename PsimagLite::Vector<ComplexOrRealType>::Type& x,
	                 const typename PsimagLite::Vector<ComplexOrRealType>::Type& y,
	                 SizeType xx,
	                 ProgramGlobals::ConnectionEnum type,
	                 const ComplexOrRealType &valuec,
	                 SizeType term,
	                 SizeType dofs,
	                 const AdditionalDataType& additionalData) const
	{
		SparseMatrixType const* A = 0;
		SparseMatrixType const* B = 0;
		LinkType link2 = hc_.getKron(&A,&B,xx,type,valuec,term,dofs,additionalData);
		hc_.modelHelper().fastOpProdInter(x,y,*A,*B,link2);
	}

	VectorType& x_;
	const VectorType& y_;
	const HamiltonianConnectionType& hc_;
	typename PsimagLite::Vector<VectorType>::Type xtemp_;
};
}
#endif // PARALLELHAMILTONIANCONNECTION_H
