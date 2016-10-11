#ifndef FEASJZSYMMETRY_H
#define FEASJZSYMMETRY_H
#include "Vector.h"
#include "CrsMatrix.h"
#include "Operator.h"
#include "ModelBase.h"
#include "HilbertSpaceFeAs.h"

namespace Dmrg {

template<typename HilbertBasisType, typename VectorOperatorType, bool>
class FeAsJzSymmetry {
public:

	typedef typename VectorOperatorType::value_type OperatorType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;


	FeAsJzSymmetry(bool) {}
	void init(HilbertBasisType&,
	          VectorOperatorType&) {}

	void setElectronsAndJz(VectorSizeType& ,
	                       VectorSizeType& ) const {}

	void findElectrons(VectorSizeType& ,
	                   const HilbertBasisType& ,
	                   SizeType) const {}

	void jzReinterpret(MatrixType& ) const {}

	bool isEnabled() const { return false; }

	bool isSet() const { return false; }

};

template<typename HilbertBasisType, typename VectorOperatorType>
class FeAsJzSymmetry<HilbertBasisType,VectorOperatorType, true> {

public:

	typedef typename VectorOperatorType::value_type OperatorType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef  typename HilbertBasisType::value_type WordType;
	typedef  HilbertSpaceFeAs<WordType> HilbertSpaceFeAsType;

	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

	FeAsJzSymmetry(bool isEnabled) :
	    isEnabled_(isEnabled), isSet_(false)
	{}

	void init(HilbertBasisType& natBasis,
	          VectorOperatorType& creationMatrix)
	{
		assert(!isSet_);
		// write operator Jz in first basis
		MatrixType Jz_opr = Jz_opr_original_basis(creationMatrix);



		// reorder Jz so that it is block diagonal in n and find P,and determine the electrons_ and save it --> permutation
		VectorRealType blockOffsets;

		MatrixType P(natBasis.size(),natBasis.size());
		Get_P_and_Blocks_and_electrons(natBasis,blockOffsets,P);

		Rotate(P,Jz_opr);
		// diagonalize each block --> U, jzEigs
		VectorRealType jzEigs;

		//DiagonalizeBlocks_GetU(Jz_opr, blockOffsets, jzEigs);

		// convertJzEigs(jzModifiedEigs_,jzEigs);

		// do U' = U*P^\dagger --> u_ and utranspose_
		//Get_UP(UP,P);
		//rotate all operators
		//Rotate(UP,creationMatrix);

		isSet_ = true;
	}

	void setElectronsAndJz(VectorSizeType& electrons,
	                       VectorSizeType& electronsUp) const
	{
		if (!isEnabled_) return;
		assert(isSet_);
		electrons = electrons_;
		electronsUp = jzModifiedEigs_;
	}

	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType) const
	{
		if (!isEnabled_) return;
		assert(isSet_);
		if (basis.size() != electrons_.size())
			throw PsimagLite::RuntimeError("FeAsJzSymmetry: findElectrons");
		electrons = electrons_;
	}

	void jzReinterpret(MatrixType& cm) const
	{
		if (!isEnabled_) return;
		assert(isSet_);
		MatrixType tmp = utranspose_*cm;
		cm = tmp*u_;
	}







	bool isEnabled() const { return isEnabled_; }

	bool isSet() const { return isSet_; }

private:

	MatrixType Jz_opr_original_basis(VectorOperatorType& creationMatrix)
	{
		SizeType nrow = creationMatrix[0].data.row();
		MatrixType tmp(nrow,nrow);
		//Works only for 3 orbital model
		//xy=0,zx=1,yz=2
		//nxy_up-nxy_dn
		tmp += multiplyTc(creationMatrix[0].data,creationMatrix[0].data);
		tmp -= multiplyTc(creationMatrix[3].data,creationMatrix[3].data)  ;

		//nzx_dn-nzx_up
		tmp += multiplyTc(creationMatrix[4].data,creationMatrix[4].data);
		tmp-= multiplyTc(creationMatrix[1].data,creationMatrix[1].data) ;

		//nyz_dn-nyz_up
		tmp += multiplyTc(creationMatrix[5].data,creationMatrix[5].data) ;
		tmp -= multiplyTc(creationMatrix[2].data,creationMatrix[2].data) ;


		//i(czx_dn_dagg*cyz_dn) + hc (first with i, sec with -i, std::complex<double>(0.0,1.0)*)
		ComplexOrRealType sqrtMinus1(0,1);

		tmp += sqrtMinus1*(multiplyTc(creationMatrix[4].data,creationMatrix[5].data) );

		tmp += -sqrtMinus1*(multiplyTc(creationMatrix[5].data,creationMatrix[4].data) ) ;

		//i(czx_up_dagg*cyz_up) + hc
		tmp += (multiplyTc(creationMatrix[1].data,creationMatrix[2].data));

		tmp += (multiplyTc(creationMatrix[2].data,creationMatrix[1].data) ) ;

		return tmp;

	}


	void Get_P_and_Blocks_and_electrons(const HilbertBasisType& natBasis, VectorRealType& blockOffsets, MatrixType& P)
	{
		//HilbertBasisType newBasis;
		electrons_.resize(natBasis.size());

		//Works only for 3 orbital model
		SizeType j=0;
		for(SizeType ne=0;ne<7;ne++){
			for(SizeType i=0;i<natBasis.size();i++){
				if(no_of_electrons(natBasis[i])==ne){
					//newBasis.pushback(natBasis[i]);
					P(j,i)=1;
					electrons_[j]=ne;
					j=j+1;
				}
			}
			blockOffsets.push_back(j);

		}






	}



	void DiagonalizeBlocks_GetU(const MatrixType& Jz_opr,
	                            const VectorRealType& blockOffsets,
	                            VectorRealType& jzEigs)
	{
		//        SizeType nrow = Jz_opr.row();

		//    u_.resize(nrow,nrow);
		//    utranspose_.resize(nrow,nrow);



	}


	SizeType no_of_electrons(WordType basis_i)
	{

		SizeType tmp_e;
		// nup
		SizeType nup = HilbertSpaceFeAsType::electronsWithGivenSpin(basis_i,
		                                                            SPIN_UP);
		// ndown
		SizeType ndown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis_i,
		                                                              SPIN_DOWN);
		tmp_e = nup + ndown;
		return tmp_e;
	}


	void Rotate(const MatrixType& R, MatrixType& O)
	{
		SizeType nrow = R.n_row();
		MatrixType tmp(nrow,nrow);
		tmp = O*R;
		MatrixType R_dagg(nrow,nrow);
		R_dagg=transposeConjugate(R);
		O = R_dagg*tmp;
	}


	void convertJzEigs(VectorSizeType& electronselectronsUp,
	                   const VectorRealType& jzEigs) const
	{

	}

	bool isEnabled_;
	bool isSet_;
	MatrixType u_;
	MatrixType utranspose_;
	VectorSizeType jzModifiedEigs_;
	VectorSizeType electrons_;
}; // class FeAsJzSymmetry
} // namespace Dmrg
#endif // FEASJZSYMMETRY_H
