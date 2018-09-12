#ifndef PSI_PSIMAGLITE_H
#define PSI_PSIMAGLITE_H

#include <iostream>
#include <utility>
#include "Concurrency.h"
#include "AnsiColors.h"
#include "TypeToString.h"
#include "Vector.h"
#include "Random48.h"

namespace PsimagLite {

std::ostream& operator<<(std::ostream&,const std::pair<SizeType,SizeType>&);

std::istream& operator>>(std::istream&,std::pair<SizeType,SizeType>&);

SizeType log2Integer(SizeType x);

void err(String);

struct MatchPathSeparator {
    bool operator()(char ch) const
    {
        return (ch == '/');
    }
};

template<typename X,typename A>
void fillRandom(std::vector<X, A>& v)
{
	SizeType n = v.size();
	if (n == 0)
		throw std::runtime_error("fillRandom must be called with size > 0\n");

	Random48<X> myrng(time(0));
	for (SizeType i = 0; i < n; ++i)
		v[i] = myrng() - 0.5;
}

void split(Vector<String>::Type& tokens, String str, String delimiters = " ");

String basename(const String&);

class PsiApp {
public:

	PsiApp(String appName, int* argc, char*** argv, int nthreads)
	    : concurrency_(argc,argv,nthreads), appName_(basename(appName))
	{
		chekSizeType();
	}

	const String& name() const { return appName_; }

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
};
} // namespace PsimagLite

void err(PsimagLite::String);

#endif // PSI_PSIMAGLITE_H

