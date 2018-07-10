#ifndef PARALLELHAMILTONIANCONNECTION_H
#define PARALLELHAMILTONIANCONNECTION_H

namespace Dmrg {

class Loop {

public:

	Loop()
	    :
	      xtemp_(ConcurrencyType::storageSize(ConcurrencyType::codeSectionParams.npthreads))
	{}

	      void doTask(SizeType taskNumber ,SizeType threadNum)
	{
		if (xtemp_[threadNum].size() != x_.size())
			xtemp_[threadNum].resize(x_.size(),0.0);

		ComplexOrRealType tmp = 0.0;

		if (taskNumber == 0) {
			modelHelper_.hamiltonianLeftProduct(xtemp_[threadNum],y_);
			return;
		}

		if (taskNumber == 1) {
			modelHelper_.hamiltonianRightProduct(xtemp_[threadNum],y_);
			return;
		}

		assert(taskNumber > 1);
		taskNumber -= 2;
		AdditionalDataType additionalData;
		SizeType xx = 0;
		ProgramGlobals::ConnectionEnum type;
		SizeType term = 0;
		SizeType dofs =0;
		prepare(xx,type,tmp,term,dofs,additionalData,taskNumber);

		linkProduct(xtemp_[threadNum],y_,xx,type,tmp,term,dofs,additionalData);
	}

	SizeType tasks() const { return total_; }

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
		LinkType link2 = getKron(&A,&B,xx,type,valuec,term,dofs,additionalData);
		modelHelper_.fastOpProdInter(x,y,*A,*B,link2);
	}

	VectorType& x_;
	const VectorType& y_;
	VectorVectorType xtemp_;
};
}
#endif // PARALLELHAMILTONIANCONNECTION_H
