#include "Options.h"
#include <iostream>
#include "String.h"

int main(int argc,char* argv[])
{
	if (argc<2) {
		std::cerr<<"USAGE is "<<argv[0]<<" comma,separated,list,of,options\n";
		return 1;
	}
	PsimagLite::Vector<PsimagLite::String>::Type registerOpts;
	registerOpts.push_back("fast");
	registerOpts.push_back("verbose");
	registerOpts.push_back("hasthreads");
	PsimagLite::Options::Writeable optWriteable(registerOpts,PsimagLite::Options::Writeable::STRICT);
	
	PsimagLite::String myoptions(argv[1]);
	PsimagLite::Options::Readable optsReadable(optWriteable,myoptions);
	std::cout<<"fast="<<optsReadable.isSet("fast")<<"\n";
}

