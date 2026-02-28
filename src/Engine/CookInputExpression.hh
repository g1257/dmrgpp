#ifndef COOKINPUTEXPRESSION_HH
#define COOKINPUTEXPRESSION_HH
#include "InputCheck.h"
#include "InputNg.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include <string>
#include <vector>

namespace Dmrg {

template <typename ComplexOrRealType> class CookInputExpression {
public:

	using RealType             = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType          = PsimagLite::InputNg<Dmrg::InputCheck>;
	using InputNgValidatorType = InputNgType::Readable;

	CookInputExpression(const InputNgValidatorType& io)
	    : io_(io)
	{ }

	std::string operator()(const std::string& expr)
	{
		std::string label = "!readTable";
		if (expr.substr(0, label.size()) == label) {
			std::string str = expr.substr(label.size(), std::string::npos);
			// delete spaces FIXME TODO
			// split on comma
			std::vector<std::string> args;
			PsimagLite::split(args, str, ",");
			if (args.size() != 2) {
				err("readTable expects two arguments\n");
			}

			// delete PARENS FIXME TODO
			// delete double quotes FIXME TODO
			PsimagLite::Matrix<RealType> matrix;
			InputNgValidatorType& io_non_const = const_cast<InputNgValidatorType&>(io_);
			io_non_const.read(matrix, args[0]);
			RealType value = findValueFor(matrix, PsimagLite::atof(args[1]));
			return ttos(value);
		} else {
			return expr;
		}
	}

private:

	static RealType findValueFor(const PsimagLite::Matrix<RealType>& matrix, const RealType& t)
	{
		SizeType rows = matrix.rows();
		if (matrix.cols() != 2) {
			err("findValueFor(): not a table\n");
		}

		for (SizeType i = 0; i < rows; ++i) {
			if (matrix(i, 0) == t) {
				return matrix(i, 1);
			}
		}

		throw std::runtime_error("Value not found in table\n");
	}

	const InputNgValidatorType& io_;
};
}
#endif // COOKINPUTEXPRESSION_HH
