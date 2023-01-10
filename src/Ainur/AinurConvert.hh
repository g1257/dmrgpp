#ifndef AINURCONVERT_HH
#define AINURCONVERT_HH
#include "../Matrix.h"
#include "../PsimagLite.h"

namespace PsimagLite {

class AinurConvert {

	template<typename T>
	struct Action {

		Action(String name, std::vector<T>& t)
		    : name_(name), t_(t)
		{}

		template <typename A, typename ContextType>
		void operator()(A& attr,
		                ContextType&,
		                bool&) const;

	private:

		String name_;
		std::vector<T>& t_;
	}; // struct Action

	template<typename T>
	struct ActionMatrix {

		ActionMatrix(String name, Matrix<T>& t)
		    : name_(name), t_(t)
		{}

		template <typename A, typename ContextType>
		void operator()(A& attr,
		                ContextType&,
		                bool&) const;

	private:

		String name_;
		Matrix<T>& t_;
	}; // struct ActionMatrix

public:

	template<typename T>
	static void convert(std::vector<T>& t,
	                    String value,
	                    typename EnableIf<Loki::TypeTraits<T>::isArith ||
	                    IsComplexNumber<T>::True ||
	                    TypesEqual<T, String>::True,
	                    int>::Type = 0);

	template<typename T>
	static void convert(Matrix<T>& t,
	                    String value);

	template<typename T>
	static void convert(T& t,
	                    String label,
	                    typename EnableIf<Loki::TypeTraits<T>::isIntegral,
	                    int>::Type = 0)
	{
		try {
			t = PsimagLite::atoi(label.c_str());
		} catch (std::exception& e) {
			std::cerr<<"FATAL: AinurState: Label " + label + " must be an integer\n";
			throw e.what();
		}
	}

	template<typename T>
	static void convert(T& t,
	                    String label,
	                    typename EnableIf<Loki::TypeTraits<T>::isFloat,
	                    int>::Type = 0)
	{
		try {
			t = PsimagLite::atof(label.c_str());
		} catch (...) {
			err("FATAL: AinurState: Label " + label + " must be a real number\n");
		}
	}

	static void convert(String& t, String label)
	{
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
};
}
#endif // AINURCONVERT_HH
