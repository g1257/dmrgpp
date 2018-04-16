#include "Io/IoNg.h"

void test1()
{
	PsimagLite::IoNg::Out ioOut("hello.hdf5");

	std::vector<double> v(10, 42.0);
	ioOut.write(v, "MyVector");

	ioOut.close();

	sleep(1);

	PsimagLite::IoNg::In ioIn("hello.hdf5");

	std::vector<double> w;

	ioIn.read(w, "MyVector");

	std::cout<<w;
}


int main(int argc, char** argv)
{
	test1();
}
