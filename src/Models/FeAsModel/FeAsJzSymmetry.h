#ifndef FEASJZSYMMETRY_H
#define FEASJZSYMMETRY_H
#include "Vector.h"
#include "CrsMatrix.h"
#include "Operator.h"
#include "ModelBase.h"
#include "HilbertSpaceFeAs.h"
#include "Matrix.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename HilbertBasisType, typename VectorOperatorType, bool>
class FeAsJzSymmetry {
public:

	typedef typename VectorOperatorType::value_type OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef typename OperatorStorageType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;


	FeAsJzSymmetry(bool) {}
	void init(HilbertBasisType&,
	          VectorOperatorType&) {}

	void setElectronsAndJz(SizeType& electrons,
	                       SizeType& electronsUp,
	                       SizeType ind) const {}

	bool isEnabled() const { return false; }

	bool isSet() const { return false; }

	void write(PsimagLite::String,
	           PsimagLite::IoNg::Out::Serializer&) const
	{}
};

template<typename HilbertBasisType, typename VectorOperatorType>
class FeAsJzSymmetry<HilbertBasisType,VectorOperatorType, true> {

public:

	typedef typename VectorOperatorType::value_type OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef typename OperatorStorageType::value_type ComplexOrRealType;
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
	{
		if (!isEnabled_) return;
		ProgramGlobals::oldChangeOfBasis = true;
	}

	void init(const HilbertBasisType& natBasis,
	          VectorOperatorType& creationMatrix)
	{
		assert(!isSet_);
		// write operator Jz in first basis
		MatrixType Jz_opr = Jz_opr_original_basis(creationMatrix);

		// reorder Jz so that it is block diagonal in n and find P,
		// and determine the electrons_ and save it --> permutation
		VectorRealType blockOffsets;

		MatrixType P(natBasis.size(),natBasis.size());
		MatrixType P_dagg(natBasis.size(),natBasis.size());
		Get_P_and_Blocks_and_electrons(natBasis,blockOffsets,P);
		transposeConjugate(P_dagg,P);

		//Rotate(P,Jz_opr);
		// diagonalize each block --> U, jzEigs
		VectorRealType jzEigs;
		jzEigs.resize(natBasis.size());

		DiagonalizeBlocks_GetU(Jz_opr, blockOffsets, jzEigs);

		convertJzEigs(jzModifiedEigs_,jzEigs);

		MatrixType UP(natBasis.size(),natBasis.size());
		UP=P_dagg*u_;

		//rotate all operators
		Rotate_all(u_,creationMatrix);

		isSet_ = true;
	}

	void setElectronsAndJz(SizeType& electrons,
	                       SizeType& electronsUp,
	                       SizeType ind) const
	{
		if (!isEnabled_) return;
		assert(isSet_);
		assert(ind < electrons_.size());
		electrons = electrons_[ind];
		assert(ind < jzModifiedEigs_.size());
		electronsUp = jzModifiedEigs_[ind];
	}

	bool isEnabled() const { return isEnabled_; }

	bool isSet() const { return isSet_; }

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/FeAsJzSymmetry";
		io.createGroup(label);
		io.write(label + "/isEnabled_", isEnabled_);
		if (!isEnabled_) return;
		io.write(label + "/isSet_", isSet_);
		if (!isSet_) return;
		u_.write(label + "/u_", io);
		utranspose_.write(label + "/utranspose_", io);
		Hamil_onsite_.write(label + "/Hamil_onsite_", io);
		io.write(label + "/jzModifiedEigs_", jzModifiedEigs_);
		io.write(label + "/electrons_", electrons_);
	}

private:

	MatrixType Jz_opr_original_basis(VectorOperatorType& creationMatrix)
	{
		SizeType nrow = creationMatrix[0].getStorage().rows();
		MatrixType tmp(nrow,nrow);
		//Works only for 3 orbital model
		//xy=0,zx=1,yz=2
		//nxy_up-nxy_dn
		tmp += 0.5*(multiplyTc(creationMatrix[0].getStorage(),creationMatrix[0].getStorage()));
		tmp += -0.5*(multiplyTc(creationMatrix[3].getStorage(),creationMatrix[3].getStorage()));


		//        //nzx_up-nzx_dn
		tmp += -0.5*(multiplyTc(creationMatrix[4].getStorage(),creationMatrix[4].getStorage()));
		tmp += 0.5*(multiplyTc(creationMatrix[1].getStorage(),creationMatrix[1].getStorage())) ;

		//        //nyz_up-nyz_dn
		tmp += -0.5*(multiplyTc(creationMatrix[5].getStorage(),creationMatrix[5].getStorage())) ;
		tmp += 0.5*(multiplyTc(creationMatrix[2].getStorage(),creationMatrix[2].getStorage())) ;

		//i(czx_dn_dagg*cyz_dn) + hc (first with i, sec with -i,
		ComplexOrRealType sqrtMinus1(0,1);

		tmp += +sqrtMinus1*(multiplyTc(creationMatrix[4].getStorage(),creationMatrix[5].getStorage()));

		tmp += -sqrtMinus1*(multiplyTc(creationMatrix[5].getStorage(),creationMatrix[4].getStorage()));

		//i(czx_up_dagg*cyz_up) + hc
		tmp += +sqrtMinus1*(multiplyTc(creationMatrix[1].getStorage(),creationMatrix[2].getStorage()));

		tmp += -sqrtMinus1*(multiplyTc(creationMatrix[2].getStorage(),creationMatrix[1].getStorage()) ) ;

		return tmp;

	}

	MatrixType Calculate_onsite_Hamiltonian(VectorOperatorType& creationMatrix){
		SizeType nrow = creationMatrix[0].getStorage().row();
		MatrixType tmp(nrow,nrow);
		MatrixType tmp2(nrow,nrow),tmp3(nrow,nrow),tmp4(nrow,nrow);
		RealType _U,_J,_Up;
		ComplexOrRealType sqrtMinus1(0,1);

		//    tmp += -0.5*sqrtMinus1*(multiplyTc(creationMatrix[1].getStorage(),creationMatrix[2].getStorage()));
		//    tmp +=  0.5*sqrtMinus1*(multiplyTc(creationMatrix[1].getStorage(),creationMatrix[3].getStorage()));
		//    tmp +=  0.5*sqrtMinus1*(multiplyTc(creationMatrix[2].getStorage(),creationMatrix[1].getStorage()));
		//    tmp +=  -0.5*(multiplyTc(creationMatrix[2].getStorage(),creationMatrix[3].getStorage()));
		//    tmp += -0.5*sqrtMinus1*(multiplyTc(creationMatrix[0].getStorage(),creationMatrix[4].getStorage()));
		//    tmp += 0.5*(multiplyTc(creationMatrix[0].getStorage(),creationMatrix[5].getStorage()));

		//    tmp += 0.5*sqrtMinus1*(multiplyTc(creationMatrix[4].getStorage(),creationMatrix[0].getStorage()));
		//    tmp += 0.5*sqrtMinus1*(multiplyTc(creationMatrix[4].getStorage(),creationMatrix[5].getStorage()));
		//    tmp += 0.5*(multiplyTc(creationMatrix[5].getStorage(),creationMatrix[0].getStorage()));
		//    tmp += -0.5*sqrtMinus1*(multiplyTc(creationMatrix[5].getStorage(),creationMatrix[4].getStorage()));
		//    tmp += -0.5*sqrtMinus1*(multiplyTc(creationMatrix[3].getStorage(),creationMatrix[1].getStorage()));
		//    tmp += -0.5*(multiplyTc(creationMatrix[3].getStorage(),creationMatrix[2].getStorage()));


		//U=1.0, J=0.25 ,Up=U-2J
		_U=2.0;_J=0.5;_Up=_U-2*_J;


		//U terms
		if(true){
			//n_up_xy*n_dn_xy
			tmp2 =  ( (multiplyTc(creationMatrix[0].getStorage(),creationMatrix[0].getStorage())) * ( multiplyTc(creationMatrix[3].getStorage(),creationMatrix[3].getStorage())) );
			tmp += _U*tmp2;
			//n_up_xz*n_dn_xz
			tmp2 = ( (multiplyTc(creationMatrix[1].getStorage(),creationMatrix[1].getStorage())) * ( multiplyTc(creationMatrix[4].getStorage(),creationMatrix[4].getStorage())) );
			tmp += _U*tmp2;
			//n_up_yz*n_dn_yz
			tmp2 = ( (multiplyTc(creationMatrix[2].getStorage(),creationMatrix[2].getStorage())) * ( multiplyTc(creationMatrix[5].getStorage(),creationMatrix[5].getStorage())) );
			tmp += _U*tmp2;
		}

		//Up-J/2 terms
		if(true){
			//nxy*nxz
			tmp2 = (multiplyTc(creationMatrix[0].getStorage(),creationMatrix[0].getStorage()));
			tmp2 += (multiplyTc(creationMatrix[3].getStorage(),creationMatrix[3].getStorage()));
			tmp3 = (multiplyTc(creationMatrix[1].getStorage(),creationMatrix[1].getStorage()));
			tmp3 += ( multiplyTc(creationMatrix[4].getStorage(),creationMatrix[4].getStorage()));
			tmp4 = tmp2*tmp3;

			tmp += (_Up - _J*(0.5))*tmp4;
			//nxy*nyz
			tmp2 = (multiplyTc(creationMatrix[0].getStorage(),creationMatrix[0].getStorage()));
			tmp2 += (multiplyTc(creationMatrix[3].getStorage(),creationMatrix[3].getStorage()));
			tmp3 = (multiplyTc(creationMatrix[2].getStorage(),creationMatrix[2].getStorage()));
			tmp3 += ( multiplyTc(creationMatrix[5].getStorage(),creationMatrix[5].getStorage()));
			tmp4 = tmp2*tmp3;

			tmp += (_Up - _J*(0.5))*tmp4;


			//nxz*nyz
			tmp2 = (multiplyTc(creationMatrix[1].getStorage(),creationMatrix[1].getStorage()));
			tmp2 += (multiplyTc(creationMatrix[4].getStorage(),creationMatrix[4].getStorage()));
			tmp3 = (multiplyTc(creationMatrix[2].getStorage(),creationMatrix[2].getStorage()));
			tmp3 += ( multiplyTc(creationMatrix[5].getStorage(),creationMatrix[5].getStorage()));
			tmp4 = tmp2*tmp3;

			tmp += (_Up - _J*(0.5))*tmp4;
		}



		//SzSz Hunds term
		if (true){
			//Sz_xy*Sz_xz
			tmp2 = (multiplyTc(creationMatrix[0].getStorage(),creationMatrix[0].getStorage()));
			tmp2 += (-1.0)*(multiplyTc(creationMatrix[3].getStorage(),creationMatrix[3].getStorage()));
			tmp3 = (multiplyTc(creationMatrix[1].getStorage(),creationMatrix[1].getStorage()));
			tmp3 += (-1.0)*( multiplyTc(creationMatrix[4].getStorage(),creationMatrix[4].getStorage()));
			tmp4 = tmp2*tmp3;

			tmp += (0.25)*(-2*_J)*tmp4;


			//Sz_xy*Sz_yz
			tmp2 = (multiplyTc(creationMatrix[0].getStorage(),creationMatrix[0].getStorage()));
			tmp2 += (-1.0)*(multiplyTc(creationMatrix[3].getStorage(),creationMatrix[3].getStorage()));
			tmp3 = (multiplyTc(creationMatrix[2].getStorage(),creationMatrix[2].getStorage()));
			tmp3 += (-1.0)*( multiplyTc(creationMatrix[5].getStorage(),creationMatrix[5].getStorage()));
			tmp4 = tmp2*tmp3;

			tmp += (0.25)*(-2*_J)*tmp4;


			//Sz_xz*Sz_yz
			tmp2 = (multiplyTc(creationMatrix[1].getStorage(),creationMatrix[1].getStorage()));
			tmp2 += (-1.0)*(multiplyTc(creationMatrix[4].getStorage(),creationMatrix[4].getStorage()));
			tmp3 = (multiplyTc(creationMatrix[2].getStorage(),creationMatrix[2].getStorage()));
			tmp3 += (-1.0)*( multiplyTc(creationMatrix[5].getStorage(),creationMatrix[5].getStorage()));
			tmp4 = tmp2*tmp3;

			tmp += (0.25)*(-2*_J)*tmp4;
		}

		//S+S- + S-S+ term

		if(true){

			//Sp_xy*Sm_xz + Sp_xz*Sm_xy
			tmp2 = multiplyTc(creationMatrix[0].getStorage(),creationMatrix[3].getStorage());
			tmp3 = multiplyTc(creationMatrix[4].getStorage(),creationMatrix[1].getStorage());
			tmp4 = tmp2*tmp3;
			tmp += (-_J)*tmp4;

			tmp2 = multiplyTc(creationMatrix[1].getStorage(),creationMatrix[4].getStorage());
			tmp3 = multiplyTc(creationMatrix[3].getStorage(),creationMatrix[0].getStorage());
			tmp4 = tmp2*tmp3;
			tmp += (-_J)*tmp4;

			//Sp_xy*Sm_yz + Sp_yz*Sm_xy
			tmp2 = multiplyTc(creationMatrix[0].getStorage(),creationMatrix[3].getStorage());
			tmp3 = multiplyTc(creationMatrix[5].getStorage(),creationMatrix[2].getStorage());
			tmp4 = tmp2*tmp3;
			tmp += (-_J)*tmp4;

			tmp2 = multiplyTc(creationMatrix[2].getStorage(),creationMatrix[5].getStorage());
			tmp3 = multiplyTc(creationMatrix[3].getStorage(),creationMatrix[0].getStorage());
			tmp4 = tmp2*tmp3;
			tmp += (-_J)*tmp4;

			//Sp_yz*Sm_xz + Sp_xz*Sm_yz
			tmp2 = multiplyTc(creationMatrix[2].getStorage(),creationMatrix[5].getStorage());
			tmp3 = multiplyTc(creationMatrix[4].getStorage(),creationMatrix[1].getStorage());
			tmp4 = tmp2*tmp3;
			tmp += (-_J)*tmp4;

			tmp2 = multiplyTc(creationMatrix[1].getStorage(),creationMatrix[4].getStorage());
			tmp3 = multiplyTc(creationMatrix[5].getStorage(),creationMatrix[2].getStorage());
			tmp4 = tmp2*tmp3;
			tmp += (-_J)*tmp4;
		}

		//Pair hopping term

		if (true){
			//P_xy_dagg P_xz + P_xz_dagg P_xy , gamma=xy,gamma_p=xz

			tmp2= (multiplyTc(creationMatrix[0].getStorage(),creationMatrix[4].getStorage()));
			tmp3= (multiplyTc(creationMatrix[1].getStorage(),creationMatrix[3].getStorage()));

			tmp4 = tmp2*transposeConjugate(tmp3);
			tmp += (-1.0)*(_J)*tmp4;

			tmp4 = tmp3*transposeConjugate(tmp2);
			tmp += (-1.0)*(_J)*tmp4;

			//gamma=xy,gamma_p=yz
			tmp2= (multiplyTc(creationMatrix[0].getStorage(),creationMatrix[5].getStorage()));
			tmp3= (multiplyTc(creationMatrix[2].getStorage(),creationMatrix[3].getStorage()));

			tmp4 = tmp2*transposeConjugate(tmp3);
			tmp += (-1.0)*(_J)*tmp4;

			tmp4 = tmp3*transposeConjugate(tmp2);
			tmp += (-1.0)*(_J)*tmp4;


			//gamma=xz,gamma_p=yz
			tmp2= (multiplyTc(creationMatrix[1].getStorage(),creationMatrix[5].getStorage()));
			tmp3= (multiplyTc(creationMatrix[2].getStorage(),creationMatrix[4].getStorage()));

			tmp4 = tmp2*transposeConjugate(tmp3);
			tmp += (-1.0)*(_J)*tmp4;

			tmp4 = tmp3*transposeConjugate(tmp2);
			tmp += (-1.0)*(_J)*tmp4;


		}

		return tmp;
	}

	void Get_P_and_Blocks_and_electrons(const HilbertBasisType& natBasis,
	                                    VectorRealType& blockOffsets,
	                                    MatrixType& P)
	{
		//HilbertBasisType newBasis;
		electrons_.resize(natBasis.size());

		//Works only for 3 orbital model
		SizeType j=0;
		for (SizeType ne=0;ne<7;ne++){
			for (SizeType i=0;i<natBasis.size();i++){
				if (no_of_electrons(natBasis[i])==ne){
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
		SizeType nrow = Jz_opr.n_row();

		u_.resize(nrow,nrow);
		utranspose_.resize(nrow,nrow);

		VectorRealType Jz_block_eigs;
		SizeType r_=0;
		for(SizeType i=0;i<7;i++){
			SizeType nrow_b=blockOffsets[i]-r_;

			MatrixType Jz_block(nrow_b, nrow_b);
			for(SizeType ir=0;ir<nrow_b;ir++){
				for(SizeType ic=0;ic<nrow_b;ic++){
					Jz_block(ir,ic)=Jz_opr(r_+ir,r_+ic);
				}
			}

			diag(Jz_block,Jz_block_eigs,'V');

			for(SizeType ir=0;ir<nrow_b;ir++){
				jzEigs[r_+ir]=Jz_block_eigs[ir];
				for(SizeType ic=0;ic<nrow_b;ic++){
					u_(r_+ir,r_+ic)=Jz_block(ir,ic);
					utranspose_(r_+ic,r_+ir)=Jz_block(ir,ic);

				}
			}

			r_=blockOffsets[i];
		}
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

	void Rotate_all(const MatrixType& R, VectorOperatorType& creationMatrix)
	{
		for(SizeType i = 0; i < creationMatrix.size(); ++i){
			MatrixType tmp = creationMatrix[i].getStorage().toDense()*R;
			MatrixType tmp2 = multiplyTransposeConjugate(R, tmp);
			creationMatrix[i].fromStorage(tmp2);
		}
	}

	void convertJzEigs(VectorSizeType& electronselectronsUp,
	                   const VectorRealType& jzEigs) const
	{

		electronselectronsUp.resize(jzEigs.size());
		SizeType tmp_Sztype;
		RealType tmp_doub;
		for(SizeType i=0;i<jzEigs.size();i++){
			tmp_doub=(2*jzEigs[i]);// + 5.0; //making all int value(by 2*) and positive(by +5.0)


			//Ne=1
			if(i>0 && i<7)  {
				tmp_doub += 3;
			}

			//Ne=2
			if(i>6 && i<22 ){
				tmp_doub += 6;
			}

			//Ne=3
			if(i>21 && i<42){
				tmp_doub += 9;
			}

			//Ne=4
			if(i>41 && i<57)  {
				tmp_doub += 12;
			}

			//Ne=5
			if(i>56 && i<63){
				tmp_doub += 15;
			}
			//Ne=6
			else if(i==63){
				tmp_doub += 18;
			}


			tmp_doub=tmp_doub*0.5;
			tmp_Sztype= (SizeType) (tmp_doub +0.5); //rounding off (like 0.999 to 1) before converting to SizeType


			electronselectronsUp[i]=tmp_Sztype;}



	}

	bool isEnabled_;
	bool isSet_;
	MatrixType u_;
	MatrixType utranspose_;
	MatrixType Hamil_onsite_;
	VectorSizeType jzModifiedEigs_;
	VectorSizeType electrons_;
}; // class FeAsJzSymmetry
} // namespace Dmrg
#endif // FEASJZSYMMETRY_H
