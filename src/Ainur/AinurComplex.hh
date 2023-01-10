#ifndef AINURCOMPLEX_HH
#define AINURCOMPLEX_HH

#include "../PsimagLite.h"

namespace PsimagLite {

struct AinurComplex {

	static void convert(long unsigned int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(unsigned int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(long int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(double& t, std::string str)
	{
		t = PsimagLite::atof(str);
	}

	static void convert(float& t, std::string str)
	{
		t = PsimagLite::atof(str);
	}

	template<typename T>
	static void convert(std::complex<T>& t, std::string str)
	{
		t = toComplex<T>(str);
	}

	static void convert(String& t, std::string str)
	{
		t = str;
	}

	template<typename T>
	static void convert(T& t, std::string str)
	{
		String msg("Unknown type ");
		throw RuntimeError("convert(): " + msg + typeid(t).name() + " for " + str + "\n");
	}

private:

	template<typename RealType>
	static std::complex<RealType> toComplex(std::string str)
	{
		typedef std::complex<RealType> ComplexType;

		if (str == "i") return ComplexType(0., 1.);
		if (str == "-i") return ComplexType(0., -1.);

		String buffer;
		bool flag = false;
		const SizeType n = str.length();
		RealType real1 = 0;
		for (SizeType i = 0; i < n; ++i) {
			bool isSqrtMinus1 = (str[i] == 'i');
			if (isSqrtMinus1 && flag)
				throw RuntimeError("Error parsing number " + str + "\n");

			if (isSqrtMinus1) {
				flag = true;
				real1 = atof(buffer.c_str());
				buffer = "";
				continue;
			}

			buffer += str[i];
		}

		return (flag) ? ComplexType(real1, atof(buffer.c_str())) :
		                ComplexType(atof(buffer.c_str()), 0);
	}

};

}
#endif // AINURCOMPLEX_HH
