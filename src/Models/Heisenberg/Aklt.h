#ifndef AKLT_H
#define AKLT_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ModelBaseType>
class Aklt {

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef std::pair<SizeType, char> PairSizeCharType;

	Aklt(ModelBaseType& modelBase, PsimagLite::String additional)
	    : modelBase_(modelBase), enabled_(additional == "Aklt")
	{}

	void fillLabeledOperators(SizeType site,
	                          const SparseMatrixType& splus,
	                          const SparseMatrixType& sz)
	{
		if (!enabled_) return;

		assert(site == 0);

		OpsLabelType& aklt = modelBase_.createOpsLabel("aklt");
		modelBase_.makeTrackable("aklt");

		SparseMatrixType sminus;
		transposeConjugate(sminus, splus);

		SparseMatrixType tmpMatrix = splus*splus;
		pushOneOperator(aklt, tmpMatrix);

		tmpMatrix = splus*sminus;
		pushOneOperator(aklt, tmpMatrix);

		tmpMatrix = splus*sz;
		pushOneOperator(aklt, tmpMatrix);

		tmpMatrix = sminus*splus;
		pushOneOperator(aklt, tmpMatrix);

		tmpMatrix = sminus*sminus;
		pushOneOperator(aklt, tmpMatrix);

		tmpMatrix = sminus*sz;
		pushOneOperator(aklt, tmpMatrix);

//		tmpMatrix = sz*splus; == transpose conjugate of aklt5

//		tmpMatrix = sz*sminus; == transpose conjugate of aklt2

		tmpMatrix = sz*sz;
		pushOneOperator(aklt, tmpMatrix);
	}

	void fillModelLinks()
	{
		if (!enabled_) return;
		ModelTermType& aklt = ModelBaseType::createTerm("Aklt", false);
		const auto su2prop = typename ModelTermType::Su2Properties(1, 0);
		for (SizeType mu = 0; mu < 3; ++mu) { // mu = 0 is S+, mu = 1 is S-, mu=2 is Sz
			for (SizeType mup = 0; mup < 3; ++mup) {
				const RealType factor = findFactor(mu)*findFactor(mup)/3.0;
				auto valueModifier = [factor](ComplexOrRealType& value) { value *= factor;};
				SizeType index1 = indexFor(mu, mup);
				SizeType index2 = indexFor(barOf(mu), barOf(mup));

				PairSizeCharType pair1 = operatorForIndex(index1);
				PairSizeCharType pair2 = operatorForIndex(index2);

				OpForLinkType a("aklt", pair1.first);
				OpForLinkType b("aklt", pair2.first);
				aklt.push(a, pair1.second, b, pair2.second, valueModifier, su2prop);
			}
		}
	}

private:

	static RealType findFactor(SizeType mu)
	{
		assert(mu < 3);
		return (mu == 2) ? 1 : 0.5;
	}

	static SizeType indexFor(SizeType mu1, SizeType mu2)
	{
		return mu1*3 + mu2;
	}

	static SizeType barOf(SizeType mu)
	{
		assert(mu < 3);
		if (mu == 2) return mu;
		return 1 - mu;
	}

	static PairSizeCharType operatorForIndex(SizeType index)
	{
		assert(index < 9);

		if (index < 6) return PairSizeCharType(index, 'N');

		if (index == 8) return PairSizeCharType(6, 'N');

		if (index == 6) return PairSizeCharType(5, 'C');

		assert(index == 7);
		return PairSizeCharType(2, 'C');
	}

	void pushOneOperator(OpsLabelType& aklt, const SparseMatrixType& matrix)
	{
		typename OperatorType::Su2RelatedType su2related;
		OperatorType myOp(matrix,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  PairType(0, 0),
		                  1.0,
		                  su2related);
		aklt.push(myOp);
	}

	ModelBaseType& modelBase_;
	bool enabled_;
};
}
#endif // AKLT_H
