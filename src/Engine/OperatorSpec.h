#ifndef OPERATORSPEC_H
#define OPERATORSPEC_H
#include "Io/IoSelector.h"

namespace Dmrg {

template<typename ModelType>
class OperatorSpec {

	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;

	struct NaturalOpStruct {
		NaturalOpStruct(PsimagLite::String label_)
		    : dof(0),label(label_),transpose(false)
		{
			SizeType lastIndex = label.length();
			if (lastIndex > 0) lastIndex--;
			if (label[lastIndex] == '\'') {
				label = label.substr(0,lastIndex);
				transpose = true;
			}

			label_ = label;

			SizeType i = 0;
			for (; i < label.length(); ++i) {
				if (label[i] == '?') break;
			}

			if (i == label.length()) return;

			if (i + 1 == label.length())
				err("WRONG op. spec. " + label_ + ", nothing after ?\n");

			label = label_.substr(0, i);
			dof = atoi(label_.substr(i + 1, label_.length()).c_str());
		}

		SizeType dof;
		PsimagLite::String label;
		bool transpose;
	}; // struct NaturalOpStruct

public:

	typedef OperatorType ResultType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef int AuxiliaryType;

	OperatorSpec(const ModelType& model)
	    : model_(model)
	{}

	OperatorType operator()(PsimagLite::String opLabel, int& site2) const
	{
		if (site2 < 0) site2 = extractSiteIfAny(opLabel);

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
				std::cerr<<e.what();
				throw e;
			}

			NaturalOpStruct nos(opLabel);
			nup = model_.naturalOperator(nos.label,site,nos.dof);
			if (nos.transpose)
				nup.conjugate();
		}

		SizeType foundSize = nup.data.rows();
		SizeType expectedSize = model_.hilbertSize(site);
		if (foundSize != expectedSize) {
			PsimagLite::String str("getOperatorForTest ");
			str += " Expected size " + ttos(expectedSize);
			str += " but found size " + ttos(foundSize);
			str += " for operator " + opLabel + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		return nup;
	}

	static bool metaEqual(const OperatorType& op1, const OperatorType& op2)
	{
		return (op1.metaDiff(op2) == 0);
	}

	static bool isEmpty(const OperatorType& op)
	{
		return (op.data.rows() == 0);
	}

private:

	OperatorType specialOperator(SizeType site, ComplexOrRealType value) const
	{
		SizeType rows = model_.hilbertSize(site);
		SparseMatrixType tmp(rows,rows);
		tmp.makeDiagonal(rows, value);
		typename OperatorType::Su2RelatedType su2Related;
		return OperatorType(tmp,
		                    1.0,
		                    typename OperatorType::PairType(0,0),
		                    1.0,
		                    su2Related);
	}

	int extractSiteIfAny(PsimagLite::String& name) const
	{
		int firstIndex = -1;
		int lastIndex = -1;
		for (SizeType i = 0; i < name.length(); ++i) {
			if (name[i] == '[') {
				firstIndex = i;
				continue;
			}

			if (name[i] == ']') {
				lastIndex = i;
				continue;
			}
		}

		if (firstIndex < 0 && lastIndex < 0) return -1;

		bool b1 = (firstIndex < 0 && lastIndex >= 0);
		bool b2 = (firstIndex >= 0 && lastIndex < 0);
		if (b1 || b2) {
			PsimagLite::String str("Braket operator ");
			str += name + " has unmatched [ or ]\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (static_cast<SizeType>(lastIndex) != name.length() - 1) {
			PsimagLite::String str("Braket operator ");
			str += name + " has [] but does not end in ]\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String str = name.substr(0,firstIndex);
		int site = atoi(name.substr(firstIndex+1,lastIndex-1).c_str());
		name = str;
		return site;
	}

	OperatorType findOperator(const PsimagLite::String& name,
	                          SizeType site) const
	{
		if (name.length()<2 || name[0]!=':') {
			PsimagLite::String str("OperatorInterpreter: syntax error for ");
			str += name + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String label = name.substr(1,name.length()-1);

		replaceString(label, ttos(site));
		PsimagLite::IoSelector::In io(label);

		PsimagLite::String prefix = "";
		return OperatorType(io,model_,OperatorType::MUST_BE_NONZERO, prefix);
	}

	void replaceString(PsimagLite::String& str,
	                   PsimagLite::String substr) const
	{
		/* Locate the substring to replace. */
		size_t index = str.find('$');
		if (index == PsimagLite::String::npos) return;

		PsimagLite::String str1 = str.substr(0, index);
		++index;
		PsimagLite::String str2 = str.substr(index);

		str = str1 + substr + str2;
	}

	const ModelType& model_;
};
}
#endif // OPERATORSPEC_H
