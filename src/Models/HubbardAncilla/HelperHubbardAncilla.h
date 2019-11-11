#ifndef HELPERHUBBARDANCILLA_H
#define HELPERHUBBARDANCILLA_H

#include "../FeAsModel/HilbertSpaceFeAs.h"
#include "Matrix.h"

namespace Dmrg {

template<typename ModelBaseType, typename ModelParametersType>
class HelperHubbardAncilla {

public:

	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef HilbertSpaceFeAs<HilbertState> HilbertSpaceFeAsType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename ModelHelperType::BlockType BlockType;

	static const SizeType ORBITALS = 2;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;
	static const int FERMION_SIGN = -1;

	HelperHubbardAncilla(const GeometryType& geometry, const ModelParametersType& modelParams)
	    : geometry_(geometry),
	      modelParameters_(modelParams),
	      hot_(geometry_.orbitals(0,0) > 1)
	{
		if (!hot_) return;

		PsimagLite::String msg("HubbardAncilla[Extended]: Hot ancilla mode is on");
		msg += " (EXPERIMENTAL feature)\n";
		std::cout<<msg;
		std::cerr<<msg;
	}

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io,
	           PsimagLite::String modelName) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + modelName;
		io.createGroup(label);
		modelParameters_.write(label, io);
		io.write(label + "/hot_", hot_);
	}

	bool isHot() const { return hot_; }

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType) const
	{
		SizeType n=block.size();
		HilbertBasisType natBasis;
		setBasis(natBasis, block);

		for (SizeType i=0;i<n;i++) {
			VectorSparseMatrixType cm;
			findAllMatrices(cm,i,natBasis);

			addInteraction(hmatrix, cm, block[i]);

			addPotentialV(hmatrix, cm, block[i]);
		}
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	static void setBasis(HilbertBasisType& basis,
	                     const VectorSizeType& block)
	{
		SizeType n = block.size();
		HilbertState total = (1<<(2*ORBITALS));
		total = pow(total,n);

		basis.resize(total);
		for (HilbertState a = 0; a < total; ++a) basis[a] = a;
	}


	//! set creation matrices for sites in block
	static void setLambdaMatrices(OpsLabelType& d,
	                              const VectorSparseMatrixType& vm)
	{
		typename OperatorType::Su2RelatedType su2related;
		for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
			SizeType spin2 = 1 - spin1;
			SparseMatrixType lambda;
			assert(1+spin2*ORBITALS < vm.size());
			multiply(lambda,vm[0+spin1*ORBITALS],vm[1+spin2*ORBITALS]);
			MatrixType dlambda;
			crsMatrixToFullMatrix(dlambda,lambda);
			correctLambda(dlambda, spin1, vm);

			OperatorType myOp(SparseMatrixType(dlambda),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  typename OperatorType::PairType(0,0),
			                  1,
			                  su2related);
			d.push(myOp);
		}
	}

	static void setSymmetryRelated(VectorQnType& qns,
	                               const HilbertBasisType& basis)
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		VectorSizeType other(4, 0);
		SizeType offset = basis.size();
		qns.resize(offset, QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = PairType(0,0);

			SizeType naUp = HilbertSpaceFeAsType::calcNofElectrons(basis[i],
			                                                       ORBITALS*SPIN_UP);
			SizeType naDown = HilbertSpaceFeAsType::calcNofElectrons(basis[i],
			                                                         ORBITALS*SPIN_DOWN);

			SizeType flavor = 0;

			// nup
			other[1] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                        SPIN_UP);
			// ntotal
			other[0] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                        SPIN_DOWN) + other[1];

			// up ancilla
			other[2] = naUp;

			// down ancilla
			other[3] = naDown;

			bool sign = other[0] & 1;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	static SparseMatrixType n(const SparseMatrixType& c)
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

		return tmpMatrix;
	}

	//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	static void findOperatorMatrices(MatrixType& creationMatrix,
	                                 int i,
	                                 int sigma,
	                                 const HilbertBasisType& natBasis)
	{
		HilbertState bra,ket;
		SizeType n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<n;ii++) {
			bra=ket=natBasis[ii];

			if (HilbertSpaceFeAsType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceFeAsType::create(bra,i,sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices: error\n");
				if (ii==SizeType(jj)) {
					std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra;
					std::cerr<<" sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("Creation op. cannot be diagonal\n");
				}

				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		transposeConjugate(creationMatrix,cm);
	}

	static void findAllMatrices(VectorSparseMatrixType& vm,
	                            SizeType i,
	                            const HilbertBasisType& natBasis)
	{
		for (SizeType sigma = 0; sigma < 2*ORBITALS; ++sigma) {
			MatrixType m;
			findOperatorMatrices(m,i,sigma,natBasis);
			vm.push_back(SparseMatrixType(m));
		}
	}

private:

	//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteraction(SparseMatrixType &hmatrix,
	                    const VectorSparseMatrixType& cm,
	                    SizeType actualSite) const
	{
		SparseMatrixType tmpMatrix;
		SizeType nsites = geometry_.numberOfSites();
		SizeType factor = (hot_) ? 2 : 1;
		if (modelParameters_.hubbardU.size() != factor*nsites)
			err("Number of Us is incorrect\n");

		for (SizeType alpha = 0; alpha < factor; ++alpha) {// real sites and ancillas
			SparseMatrixType m1 = cm[alpha + SPIN_UP*ORBITALS];
			SparseMatrixType m2 = cm[alpha + SPIN_DOWN*ORBITALS];

			multiply(tmpMatrix, n(m1), n(m2));
			assert(actualSite + nsites*alpha < modelParameters_.hubbardU.size());
			hmatrix += modelParameters_.hubbardU[actualSite + nsites*alpha]*tmpMatrix;
		}
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorSparseMatrixType& cm,
	                   SizeType actualIndexOfSite) const
	{
		SizeType factor = (hot_) ? 2 : 1;
		SizeType nsites = geometry_.numberOfSites();
		if (modelParameters_.potentialV.size() != factor*2*nsites)
			err("Number of Vs is incorrect\n");

		for (SizeType orbital = 0; orbital < factor; ++orbital) {
			SparseMatrixType nup = n(cm[orbital + SPIN_UP*ORBITALS]);
			SparseMatrixType ndown = n(cm[orbital + SPIN_DOWN*ORBITALS]);

			SizeType iUp = actualIndexOfSite + (orbital + 0*ORBITALS)*nsites;
			assert(iUp < modelParameters_.potentialV.size());
			hmatrix += modelParameters_.potentialV[iUp] * nup;
			SizeType iDown = actualIndexOfSite + (orbital + 1*ORBITALS)*nsites;
			assert(iDown < modelParameters_.potentialV.size());
			hmatrix += modelParameters_.potentialV[iDown] * ndown;
		}
	}

	static void correctLambda(MatrixType& dlambda,
	                          SizeType spin1,
	                          const VectorSparseMatrixType& vm)
	{
		SizeType n = dlambda.rows();
		MatrixType corrector(n,n);
		computeCorrector(corrector,spin1,vm);

		MatrixType dlambda2 = dlambda;

		dlambda = dlambda2 * corrector;
	}

	static void computeCorrector(MatrixType& corrector,
	                             SizeType spin1,
	                             const VectorSparseMatrixType& vm)
	{
		SizeType spin2 = 1 - spin1;
		SparseMatrixType cm1(vm[0+spin2*ORBITALS]);
		SparseMatrixType cm2(vm[1+spin1*ORBITALS]);
		SparseMatrixType n1 = n(cm1);
		SparseMatrixType n2 = n(cm2);
		MatrixType dn1;
		MatrixType dn2;
		crsMatrixToFullMatrix(dn1,n1);
		crsMatrixToFullMatrix(dn2,n2);

		SizeType n = corrector.rows();
		ComplexOrRealType f1 = (-1.0);
		for (SizeType i = 0; i < n; ++i)
			corrector(i,i) = std::abs(dn1(i,i) + dn2(i,i) + f1);
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to
	//! basis state ket
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	static RealType sign(HilbertState const &ket, int i,SizeType sigma)
	{
		int value=0;
		SizeType dofs = 2*ORBITALS;
		for (SizeType alpha=0;alpha<dofs;alpha++)
			value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);
		// add electron on site 0 if needed
		if (i>0) value += HilbertSpaceFeAsType::electrons(ket);

		//order for sign is: a up, a down, b up, b down, etc
		unsigned int x = HilbertSpaceFeAsType::get(ket,i);
		int spin = sigma/ORBITALS;
		SizeType orb = sigma % ORBITALS;

		for (SizeType j=0;j<orb;j++) {
			for (SizeType k=0;k<2;k++) {
				SizeType ind = j + k * ORBITALS;
				int mask = (1<<ind);
				if (x & mask) value++;
			}
		}

		if (spin==SPIN_DOWN) {
			int mask = (1<<orb);
			if (x & mask) value++;
		}

		return (value==0 || value%2==0) ? 1.0 : FERMION_SIGN;
	}

	const GeometryType& geometry_;
	const ModelParametersType& modelParameters_;
	bool hot_;
};
}
#endif // HELPERHUBBARDANCILLA_H
