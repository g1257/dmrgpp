#ifndef PSI_PSIMAGLITE_H
#define PSI_PSIMAGLITE_H

#include "MicroArchitecture.h"
#include <iostream>
#include <utility>
#include "Concurrency.h"
#include "AnsiColors.h"
#include "TypeToString.h"
#include "Vector.h"
#include "Random48.h"
#include "PsiBase64.h"
#include <fstream>
#include <sstream>

namespace PsimagLite {

std::ostream& operator<<(std::ostream&,const std::pair<SizeType,SizeType>&);

std::istream& operator>>(std::istream&,std::pair<SizeType,SizeType>&);

SizeType log2Integer(SizeType x);

int atoi(String);

int atof(String);

void err(String);

struct MatchPathSeparator {
    bool operator()(char ch) const
    {
        return (ch == '/');
    }
};

template<typename T>
void fillRandom(T& v, typename EnableIf<IsVectorLike<T>::True, int>::Type = 0)
{
	SizeType n = v.size();
	if (n == 0)
		throw std::runtime_error("fillRandom must be called with size > 0\n");

	Random48<typename T::value_type> myrng(time(0));
	typename PsimagLite::Real<typename T::value_type>::Type sum = 0;
	const typename T::value_type zeroPointFive = 0.5;
	for (SizeType i = 0; i < n; ++i) {
		v[i] = myrng() - zeroPointFive;
		sum += PsimagLite::real(v[i]*PsimagLite::conj(v[i]));
	}

	sum = 1.0/sqrt(sum);
	for (SizeType i = 0; i < n; ++i) v[i] *= sum;
}

void split(Vector<String>::Type& tokens, String str, String delimiters = " ");

String basename(const String&);

class PsiApp {

public:

	PsiApp(String appName, int* argc, char*** argv, int nthreads)
	    : concurrency_(argc, argv, nthreads),
	      appName_(basename(appName)),
	      microArch_(MicroArchitecture().vendorId())
	{
		chekSizeType();

		SizeType n = *argc;
		char** temp = *argv;
		for (SizeType i = 0; i < n; ++i)
			cmdLine_ += String(temp[i]) + " ";
	}

	void checkMicroArch(std::ostream& os, PsimagLite::String compiledArch) const
	{
		os<<"Compiled MicroArchitecture is "<<compiledArch<<"\n";
		os<<"Running on MicroArchitecture "<<microArch_<<"\n";
		if (compiledArch == microArch_) return;
		os<<"WARNING: Compiled MicroArchitecture is DIFFERENT than Running one\n";
	}

	const String& name() const { return appName_; }

	void printCmdLine(std::ostream& os) const
	{
		os<<"PsiApp: CmdLine: "<<cmdLine_<<"\n";
	}

	static void echoBase64(std::ostream& os, String filename)
	{
		os<<"PsiApp::echoBase64: Echo of "<<filename<<" in base64\n";
		PsiBase64::Encode base64(slurp(filename));
		os<<base64()<<"\n";
	}

	static String slurp(String filename)
	{
		std::ifstream fin(filename.c_str());
		std::stringstream sstr;
	    sstr << fin.rdbuf();
	    return sstr.str();
	}

private:

	void chekSizeType()
	{
		if (sizeof(SizeType) == libSizeOfSizeType_) return;
		std::string msg("PsimagLite compiled with -DUSE_LONG but");
		msg += "application without. Or viceversa.\n";
		throw std::runtime_error(msg);
	}

	static const int libSizeOfSizeType_;

	Concurrency concurrency_;
	String appName_;
	String cmdLine_;
	String microArch_;
};
} // namespace PsimagLite

void err(PsimagLite::String);

#endif // PSI_PSIMAGLITE_H

