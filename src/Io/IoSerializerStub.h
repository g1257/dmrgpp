#ifndef IOSERIALIZERSTUB_H
#define IOSERIALIZERSTUB_H

#ifdef USE_IO_NG
#include "IoNgSerializer.h"
namespace PsimagLite {
typedef IoNgSerializer IoSerializer;
}
#else
namespace PsimagLite {
class IoSerializer {

public:

	void createGroup(String root)
	{
		std::cerr<<"IoSerializer::createGroup("<<root<<"): I'm just a dummy!\n";
	}

	template<typename T>
	void write(String name2, const T& what)
	{
		std::cerr<<"IoSerializer::write("<<name2<<","<<what<<"): I'm just a dummy!\n";
	}
};
}
#endif
#endif // IOSERIALIZERSTUB_H
