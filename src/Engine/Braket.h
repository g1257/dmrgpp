#ifndef DMRG_braket_H
#define DMRG_braket_H
#include "Vector.h"
#include "Tokenizer.h"
#include "Matrix.h"

template<typename ModelType>
class Braket {

	typedef typename ModelType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<int>::Type VectorIntType;
	typedef typename OperatorType::PairType PairType;
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	struct NaturalOpStruct {
		NaturalOpStruct(PsimagLite::String label_)
		    : dof(0),label(label_),transpose(false)
		{
			SizeType i = 0;
			for (; i < label_.length(); ++i) {
				if (label_[i] == '?') break;
			}

			if (i == label_.length()) return;
			SizeType j = i;
			label = label_.substr(0,j);
			for (; i < label_.length(); ++i) {
				if (label_[i] == '-') break;
			}

			SizeType lastIndex = label.length();
			if (lastIndex > 0) lastIndex--;
			if (label[lastIndex] == '\'') {
				label = label.substr(0,lastIndex);
				transpose = true;
			}

			dof = atoi(label_.substr(j+1,i).c_str());
		}

		SizeType dof;
		PsimagLite::String label;
		bool transpose;
	}; // struct NaturalOpStruct

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	Braket(const ModelType& model,const PsimagLite::String& braket)
	    : model_(model), braket_(2,""),savedString_(braket)
	{		
		VectorStringType vecStr;
		PsimagLite::tokenizer(braket,vecStr,"|");

		if (vecStr.size() != 3) {
			PsimagLite::String str("ObserverInterpreter: syntax error for ");
			str += braket + " is not a Braket\n";
			throw PsimagLite::RuntimeError(str);
		}

		braket_[0] = vecStr[0].substr(1,vecStr[0].length()-1);
		if (!isBraket(0)) {
			PsimagLite::String str("ObserverInterpreter: syntax error: ");
			str += braket_[0] + " must be <gs or <time or <P\\d+\n";
			throw PsimagLite::RuntimeError(str);
		}

		braket_[1] = vecStr[2].substr(0,vecStr[2].length()-1);
		if (!isBraket(1)) {
			PsimagLite::String str("ObserverInterpreter: syntax error: ");
			str += braket_[1] + " must be gs> or time>  or <P\\d+\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::tokenizer(vecStr[1],name_,";");

		sites_.resize(name_.size(),-1);

		for (SizeType i = 0; i < name_.size(); ++i) {
			op_.push_back(internal_(name_[i],sites_[i]));
		}
	}

	const OperatorType& op(SizeType ind) const
	{
		assert(ind < op_.size());
		return op_[ind];
	}

	PsimagLite::String opName(SizeType ind) const
	{
		assert(ind < name_.size());
		return name_[ind];
	}

	PsimagLite::String bra() const
	{
		assert(isBraket(0));
		return braket_[0];
	}

	PsimagLite::String ket() const
	{
		assert(isBraket(1));
		return braket_[1];
	}

	SizeType points() const { return name_.size(); }

	SizeType site(SizeType ind) const
	{
		if (sites_[ind] >= 0) return sites_[ind];
		throw PsimagLite::RuntimeError("site is negative\n");
	}

	PsimagLite::String toString() const { return savedString_; }

	static int getPtype(PsimagLite::String str)
	{
		// str == P\d+
		if (str.length() < 2) return -1;
		if (str[0] != 'P') return -1;
		PsimagLite::String number("");
		for (SizeType i = 1; i < str.length(); ++i) {
			number += str[i];
			unsigned char x = str[i];
			if (x < 48 || x > 57) return -1;
		}

		return atoi(number.c_str()) + 1;
	}

private:

	bool isBraket(SizeType ind) const
	{
		if (ind >= braket_.size()) return false;
		int pType = getPtype(braket_[ind]);

		return (braket_[ind] == "gs" ||
		        braket_[ind] == "time" ||
		        pType >= 0);
	}

	OperatorType internal_(PsimagLite::String& opLabel, int& site2) const
	{
		if (site2 < 0) site2 = extractSiteIfAny(opLabel);

		SizeType site = (site2 < 0) ? 0 : site2;

		OperatorType nup;
		try {
			nup = findOperator(opLabel);
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

		SizeType foundSize = nup.data.row();
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

	OperatorType findOperator(const PsimagLite::String& name) const
	{
		if (name.length()<2 || name[0]!=':') {
			PsimagLite::String str("OperatorInterpreter: syntax error for ");
			str += name + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String label = name.substr(1,name.length()-1);

		PsimagLite::IoSimple::In io(label);

		return OperatorType(io,model_,OperatorType::MUST_BE_NONZERO);
	}

	const ModelType& model_;
	VectorStringType braket_;
	PsimagLite::String savedString_;
	VectorStringType name_;
	SizeType type_;
	VectorOperatorType op_;
	VectorIntType sites_;
}; // class Braket

#endif // DMRG_braket_H

