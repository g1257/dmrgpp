#include "Concurrency.h"
#include "Io/IoNg.h"

void test1()
{
	PsimagLite::IoNg::Out ioOut("hello.hdf5", PsimagLite::IoNg::ACC_TRUNC);

	std::vector<double> v(10, 42.0);
	ioOut.write(v, "MyVector");

	ioOut.close();

	PsimagLite::IoNg::In ioIn("hello.hdf5");

	std::vector<double> w;

	ioIn.read(w, "MyVector");

	std::cout << w;
}

int main(int argc, char* argv[])
{
	constexpr unsigned int nthreads = 1;
	PsimagLite::Concurrency(&argc, &argv, nthreads);

	test1();
}
