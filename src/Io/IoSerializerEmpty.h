#ifndef IO_SERIALIZER_EMPTY_H
#define IO_SERIALIZER_EMPTY_H
#include "../Vector.h"
#include "../TypeToString.h"
#include <cassert>
#include <stack>
#include "../Complex.h"

namespace PsimagLite {

class IoSerializerEmpty {

	typedef std::vector<unsigned char> VectorOfBoolInternalType;

public:


	enum WriteMode {NO_OVERWRITE, ALLOW_OVERWRITE};

	IoSerializerEmpty(String, unsigned int)
	{
		errorPrint("ctor");
	}

	void open(String, unsigned int)
	{
		errorPrint("open");
	}

	void close()
	{
		errorPrint("close");
	}

	void flush()
	{
		errorPrint("flush");
	}

	String filename() const
	{
		errorPrint("filename");
		return "";
	}

	void createGroup(String)
	{
		errorPrint("createGroup");
	}

	bool doesGroupExist(String)
	{
		errorPrint("doesGroupExist");
		return false;
	}

	/* write functions START */

	template<typename T>
	void write(String, const T&, WriteMode = NO_OVERWRITE)
	{
		errorPrint("write");
	}

	// Note: THIS WILL EMPTY THE STACK OBJECT!
	template<typename T>
	void write(String,
	           std::stack<T>&,
	           WriteMode = NO_OVERWRITE,
	           typename EnableIf<!Loki::TypeTraits<typename Real<T>::Type>::isArith,
	           int*>::Type = 0)
	{
		errorPrint("write");
	}

	template<typename T>
	void read(T&, String)
	{
		errorPrint("read");
	}

	static void unitTest(std::vector<bool>&)
	{
		errorPrint("unitTest");
	}

	// read functions END

private:

	IoSerializerEmpty(const IoSerializerEmpty&);

	IoSerializerEmpty& operator=(const IoSerializerEmpty&);

	static void errorPrint(String fname)
	{
		throw RuntimeError("FATAL: You called IoSerializer::" + fname
		                   + " but you compiled with USE_IO_SIMPLE."
		                   + " Please delete USE_IO_SIMPLE from Makefile and"
		                   + " enable HDF5 support.\n");
	}
};
}
#endif // IO_SERIALIZER_EMPTY_H
