#ifndef AINURCONVERT_HH
#define AINURCONVERT_HH
#include "../Matrix.h"
#include "../PsimagLite.h"
#include "AinurMacros.hh"

namespace PsimagLite {

struct AinurVariable {
	std::string key;
	std::string value;
	std::string type;
	std::string opaque;
};

class AinurConvert {

	template<typename T>
	struct Action {

		Action(String name, std::vector<T>& t, const AinurMacros& ainurMacros)
		    : name_(name), t_(t), ainurMacros_(ainurMacros)
		{}

		template <typename A, typename ContextType>
		void operator()(A& attr,
		                ContextType&,
		                bool&) const;

	private:

		String name_;
		std::vector<T>& t_;
		const AinurMacros& ainurMacros_;
	}; // struct Action

	template<typename T>
	struct ActionMatrix {

		ActionMatrix(String name, Matrix<T>& t, const AinurMacros& ainurMacros)
		    : name_(name), t_(t), ainurMacros_(ainurMacros)
		{}

		template <typename A, typename ContextType>
		void operator()(A& attr,
		                ContextType&,
		                bool&) const;

	private:

		String name_;
		Matrix<T>& t_;
		const AinurMacros& ainurMacros_;
	}; // struct ActionMatrix

public:

	AinurConvert(const AinurMacros& ainurMacros) : ainurMacros_(ainurMacros)
	{}

	template<typename T>
	void convert(std::vector<T>& t,
	             const AinurVariable& ainurVariable,
	             typename EnableIf<Loki::TypeTraits<T>::isArith ||
	             IsComplexNumber<T>::True ||
	             TypesEqual<T, String>::True,
	             int>::Type = 0);

	template<typename T>
	void convert(Matrix<T>& t,
	             const AinurVariable& ainurVariable);

	template<typename T>
	void convert(T& t,
	             const AinurVariable& ainurVariable,
	             typename EnableIf<Loki::TypeTraits<T>::isIntegral,
	             int>::Type = 0)
	{
		String label = ainurMacros_.valueFromFunction(ainurVariable.value);

		try {
			t = PsimagLite::atoi(label.c_str());
		} catch (std::exception& e) {
			std::cerr<<"FATAL: AinurState: Label " + label + " must be an integer\n";
			throw e.what();
		}
	}

	template<typename T>
	void convert(T& t,
	             const AinurVariable& ainurVariable,
	             typename EnableIf<Loki::TypeTraits<T>::isFloat,
	             int>::Type = 0)
	{
		String label = ainurMacros_.valueFromFunction(ainurVariable.value);

		try {
			t = PsimagLite::atof(label.c_str());
		} catch (...) {
			err("FATAL: AinurState: Label " + label + " must be a real number\n");
		}
	}

	void convert(String& t, const AinurVariable& ainurVariable)
	{
		String label = ainurMacros_.valueFromFunction(ainurVariable.value);
		SizeType l = label.size();
		if (l > 1 && label[0] == '"' && label[l - 1] == '"') {
			t = (l == 2) ? "" : label.substr(1,l - 2);
			return;
		}

		t = label;
	}

private:

	static String stringContext(std::string::iterator it,
	                            std::string::iterator start,
	                            std::string::iterator end,
	                            SizeType before = 5,
	                            SizeType after = 10)
	{
		std::string::iterator alpha = it;
		SizeType counter = 0;
		while (alpha != start && counter++ < before)
			--alpha;

		std::string::iterator omega = it;
		counter = 0;
		while (omega != end && counter++ < after)
			++omega;

		return String(alpha, omega);
	}

	const AinurMacros& ainurMacros_;
};
}
#endif // AINURCONVERT_HH
