#ifndef IOSERIALIZERSTUB_H
#define IOSERIALIZERSTUB_H

#ifndef USE_IO_SIMPLE
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
	void write(String name, const T& what)
	{
		std::cerr<<"IoSerializer::write("<<name<<","<<what<<"): I'm just a dummy!\n";
	}

	template<typename T>
	void read(T& what, String name)
	{
		std::cerr<<"IoSerializer::read("<<what<<","<<name<<"): I'm just a dummy!\n";
	}
};
}
#endif
#endif // IOSERIALIZERSTUB_H
