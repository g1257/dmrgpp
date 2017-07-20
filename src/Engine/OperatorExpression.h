#ifndef OPERATOREXPRESSION_H
#define OPERATOREXPRESSION_H
#include "OperatorSpec.h"
#include "PsimagLite.h"

namespace Dmrg {

template<typename ModelType>
class OperatorExpression {

	typedef OperatorSpec<ModelType> OperatorSpecType;
	typedef typename ModelType::MyBasis BasisType;
	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	OperatorExpression(const ModelType& model)
	    : operatorSpec_(model)
	{}

	// make sure that if a site is specified in an opsec, it
	// is the same in all others
	OperatorType operator()(PsimagLite::String opLabel, int& site2)
	{
		bool isSu2 = BasisType::useSu2Symmetry();
		// canonical expressions only for now
		// opLabel --> opLabelCanonical
		PsimagLite::String opLabelCanonical = opLabel;
		VectorStringType vecStr;
		OperatorType op;

		PsimagLite::split(vecStr, opLabelCanonical, "+");
		if (isSu2 && vecStr.size() > 1)
			err("Expressions not supported when SU(2) is in use\n");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			OperatorType opTerm;
			procCanonicalTerm(opTerm, vecStr[i], site2);
			if (checkMetaEqual(opTerm, op) != 0)
				err("OperatorExpression: metas not equal\n");

			if (i == 0)
				op = opTerm;
			else
				op.data += opTerm.data; // no changes to metas
		}

		bool isEmpty = (op.data.rows() == 0);
		if (isEmpty)
			err("OperatorExpression: expression has no opSec\n");

		return op;
	}

private:

	void procCanonicalTerm(OperatorType& opTerm,
	                       PsimagLite::String opTermCanonical,
	                       int& site2)
	{
		bool isSu2 = BasisType::useSu2Symmetry();
		ComplexOrRealType factor = 1.0;
		VectorStringType vecStr;
		PsimagLite::split(vecStr, opTermCanonical, "*");
		if (isSu2 && vecStr.size() > 1)
			err("Expressions not supported when SU(2) is in use\n");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			procCanonicalFactor(opTerm, factor, vecStr[i], site2);
		}

		bool isEmpty = (opTerm.data.rows() == 0);
		if (isEmpty)
			err("OperatorExpression: term has no opSec\n");

		opTerm.data *= factor;
	}

	void procCanonicalFactor(OperatorType& op,
	                         ComplexOrRealType& factor,
	                         PsimagLite::String termStr,
	                         int& site2)
	{
		bool isEmpty = (op.data.rows() == 0);
		bool isScalar =  isCanonicalScalar(termStr);
		if (isScalar) {
			ComplexOrRealType f = findFactor(termStr,
			                                 false,
			                                 static_cast<ComplexOrRealType*>(0));
			factor *= f;
			return;
		}

		OperatorType term = operatorSpec_(termStr, site2);
		if (isEmpty) {
			op = term;
			return;
		}

		if (checkMetaEqual(term, op) > 1)
			err("OperatorExpression: metas not equal\n");

		SparseMatrixType crs;
		multiply(crs, op.data, term.data);
		op.data = crs;
		op.fermionSign *= term.fermionSign;
	}

	bool isCanonicalScalar(PsimagLite::String termStr) const
	{
		if (termStr.length() == 0)
			err("OperatorExpression: operator term must no be empty\n");

		char c = termStr[0];
		bool isDigit = (c >= '0' && c <= '9');
		return (c == '.' || c == '(' || isDigit);
	}

	// Deal with complex here, FIXME
	ComplexOrRealType findFactor(PsimagLite::String termStr,
	                             bool prevHadParens,
	                             RealType* dummy) const
	{
		SizeType l = termStr.length();
		if (l == 0)
			err("OperatorExpression: scalar must no be empty\n");

		char c = termStr[0];
		bool isDigit = (c >= '0' && c <= '9');
		if (c == '.') termStr  = "0" + termStr;
		if (isDigit || c == '.') return atof(termStr.c_str());

		if (c == '-') {
			if (!prevHadParens)
				err("Negative scalar must have enclosing parens\n");

			return atof(termStr.c_str());
		}

		if (c != '(' || termStr[l-1] != ')')
			err("OperatorExpression: expected enclosing parens\n");

		PsimagLite::String tmp = termStr.substr(1, l-2);
		return findFactor(tmp, true, dummy);
	}

	ComplexOrRealType findFactor(PsimagLite::String termStr,
	                             bool prevHadParens,
	                             std::complex<RealType>*) const
	{
		std::cerr<<"WARNING: OperatorExpression: ";
		std::cerr<<"Complex scalars not yet implemented (sorry)\n";
		RealType *x = 0;
		return findFactor(termStr, prevHadParens, x);
	}

	SizeType checkMetaEqual(const OperatorType& op1, const OperatorType& op2) const
	{
		SizeType code = 0;
		VectorBoolType b(4, false);
		b[0] = (op1.fermionSign != op2.fermionSign);
		b[1] = (op1.angularFactor != op2.angularFactor);
		b[2] = (op1.jm != op2.jm);
		//b[3] = (op1.su2Related != op2.su2Related);

		SizeType orFactor = 0;
		for (SizeType i = 0; i < b.size(); ++i) {
			if (b[i]) code |= orFactor;
			orFactor <<= 1;
		}

		return code;
	}

	OperatorSpecType operatorSpec_;
};
}
#endif // OPERATOREXPRESSION_H
