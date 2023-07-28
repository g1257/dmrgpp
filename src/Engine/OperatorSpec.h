#ifndef OPERATORSPEC_H
#define OPERATORSPEC_H
#include "InputCheck.h"
#include "InputNg.h"
#include "LazyAlgebra.h"
#include "OneOperatorSpec.h"

namespace Dmrg
{

template <typename ModelType,
    typename AlgebraType = LazyAlgebra<typename ModelType::OperatorType>>
class OperatorSpec
{

public:

	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef LazyAlgebra<typename ModelType::OperatorType> LazyAlgebraType;
	typedef PsimagLite::OneOperatorSpec OneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;
	typedef AlgebraType ResultType;
	typedef typename OperatorStorageType::value_type ComplexOrRealType;
	typedef int AuxiliaryType;

	OperatorSpec(const ModelType& model)
	    : model_(model)
	{
	}

	ResultType operator()(PsimagLite::String opLabel, int& site2) const
	{
		PsimagLite::String copyOfOpLabel = opLabel;
		SiteSplitType site3Split = OneOperatorSpecType::extractSiteIfAny(opLabel);

		if (site2 >= 0 && site3Split.hasSiteString && static_cast<SizeType>(site2) != OneOperatorSpecType::strToNumberOrFail(site3Split.siteString))
			err(PsimagLite::String(__FILE__) + " FATAL , delete site from " + copyOfOpLabel + "\n");

		opLabel = site3Split.root;

		if (site2 < 0 && site3Split.hasSiteString)
			site2 = OneOperatorSpecType::strToNumberOrFail(site3Split.siteString);

		SizeType site = (site2 < 0) ? 0 : site2;

		if (opLabel == "_1" || opLabel == "identity")
			return specialOperator(site, 1.0);

		if (opLabel == "_0" || opLabel == "zero")
			return specialOperator(site, 0.0);

		OperatorType nup;

		try {
			nup = findOperator(opLabel, site);
		} catch (std::exception& e) {
			if (opLabel[0] == ':') {
				std::cerr << e.what();
				throw e;
			}

			OneOperatorSpecType nos(opLabel);
			nup = model_.naturalOperator(nos.label, site, nos.dof);
			if (nos.transpose)
				nup.dagger();
		}

		return nup;
	}

	static bool metaEqual(const ResultType& op1, const ResultType& op2)
	{
		return (op1.metaDiff(op2) == 0);
	}

	static bool isEmpty(const ResultType& op)
	{
		return op.isEmpty();
	}

private:

	OperatorType specialOperator(SizeType site, ComplexOrRealType value) const
	{
		SizeType rows = model_.hilbertSize(site);
		OperatorStorageType tmp;
		tmp.makeDiagonal(rows, value);
		typename OperatorType::Su2RelatedType su2Related;
		return OperatorType(tmp,
		    ProgramGlobals::FermionOrBosonEnum::BOSON,
		    typename OperatorType::PairType(0, 0),
		    1.0,
		    su2Related);
	}

	OperatorType findOperator(const PsimagLite::String& name,
	    SizeType site) const
	{
		if (name.length() < 2 || name[0] != ':') {
			PsimagLite::String str("OperatorInterpreter: syntax error for ");
			str += name + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String label = name.substr(1, name.length() - 1);

		replaceString(label, ttos(site));
		InputCheck inputCheck;
		PsimagLite::InputNg<InputCheck>::Writeable ioWriteable(label, inputCheck);
		PsimagLite::InputNg<InputCheck>::Readable io(ioWriteable);

		PsimagLite::String prefix = "";
		return OperatorType(io, model_, OperatorType::MUST_BE_NONZERO, prefix);
	}

	void replaceString(PsimagLite::String& str,
	    PsimagLite::String substr) const
	{
		/* Locate the substring to replace. */
		size_t index = str.find('$');
		if (index == PsimagLite::String::npos)
			return;

		PsimagLite::String str1 = str.substr(0, index);
		++index;
		PsimagLite::String str2 = str.substr(index);

		str = str1 + substr + str2;
	}

	const ModelType& model_;
};
}
#endif // OPERATORSPEC_H
