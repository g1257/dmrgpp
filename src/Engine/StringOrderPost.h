#ifndef STRINGORDERPOST_H
#define STRINGORDERPOST_H
#include "Vector.h"
#include "Braket.h"

namespace Dmrg {

template<typename ObserverType>
class StringOrderPost {

public:

	typedef typename ObserverType::ModelType ModelType;
	typedef Braket<ModelType> BraketType;
	typedef typename BraketType::VectorStringType VectorStringType;
	typedef typename BraketType::AlgebraType AlgebraType;
	typedef typename BraketType::OperatorSpecType OperatorSpecType;
	typedef typename BraketType::VectorAlgebraType VectorAlgebraType;
	typedef typename BraketType::VectorIntType VectorIntType;
	typedef typename BraketType::MatrixType MatrixType;
	typedef typename BraketType::SparseMatrixType SparseMatrixType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename BraketType::ComplexOrRealType ComplexOrRealType;

	StringOrderPost(const BraketType& braket, const ObserverType& observe)
	    : braket_(braket), observe_(observe)
	{
		static const PsimagLite::String stringop = "!stringorder=";
		const PsimagLite::String special = braket.opName(0);

		const SizeType l = stringop.length();
		if (special.substr(0, l) == stringop) {
			PsimagLite::String rest = special.substr(l, special.length() - l);
			PsimagLite::split(ops_, rest, ":");
			return;
		}

		err("FATAL: Braket: special op " + special + " not understood\n");
	}

	void computeMatrix(MatrixType& m)
	{
		const SizeType rows = m.rows();
		const SizeType cols = m.cols();
		if (rows != cols)
			err("StringOrderPost::computeMatrix() must be square\n");

		for (SizeType i = 0; i < rows; ++i) {
			for (SizeType k = i + 2; k < cols; ++k) {
				m(i, k) = computeOneMatrixElement(i, k);
			}
		}
	}

private:

	ComplexOrRealType computeOneMatrixElement(SizeType i, SizeType k) const
	{
		const SizeType n = k - i - 1;
		if (n == 0)
			err("StringOrderPost::computeOneElement\n");

		const PsimagLite::String bra = braket_.bra().toString();
		const PsimagLite::String ket = braket_.ket().toString();
		const PsimagLite::String op1 = "exp_i_pi_" + ops_[1];
		const PsimagLite::String left = "<" + bra + "|" + ops_[0] + "[" + ttos(i) + "]";
		const PsimagLite::String right = ops_[2] + "[" + ttos(k) + "]|" + ket + ">";
		if (n == 1) {
			if (k != i + 2)
				err("StringOrderPost::computeOneElement k != i + 2\n");
			// braket = op_[0]  * exp(i*pi*op_[1]) * op_[2]
			// sites are i, i + 1, k
			PsimagLite::String str(left + ";" + op1 + "[" + ttos(i + 1) + "];" + right);
			BraketType braket(braket_.model(), str);
			return observe_.threePoint(braket, 0, 0, false);
		}

		// braket = op_[0]  * exp(i*pi*op_[1]) * exp(i*pi*op_[1]) * ... exp(i*pi*op_[1]) * op_[2]
		// sites are i, i + 1, i + 2, ..., k - 1, k
		PsimagLite::String op1copies = "";
		for (SizeType j = i + 1; j < k; ++j)
			op1copies += op1 + "[" + ttos(j) + "];";

		PsimagLite::String str(left + ";" + op1copies + right); // op1copies carries a last ;
		BraketType braket(braket_.model(), str);
		return observe_.anyPoint(braket, false);
	}

	const BraketType& braket_;
	VectorStringType ops_;
	const ObserverType& observe_;
};
}
#endif // STRINGORDERPOST_H
