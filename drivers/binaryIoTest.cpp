#include "IoBinary.h"
#include "Vector.h"

void readMe(const String& myfile)
{
	PsimagLite::IoBinary::In fin(myfile);
	typename Vector<double>::Type v;
	fin.read(v,"MyVector");
	std::cout<<"MyVector\n";
	std::cout<<v;
	std::cout<<"------------------\n";

	PsimagLite::Matrix<float> m;
	fin.read(m, "MyMatrix");
	std::cout<<"MyMatrix";
	std::cout<<m;
	std::cout<<"------------------\n";
}

int main()
{
	SizeType rank = 0;
	String myfile = "myfile.txt";
	PsimagLite::IoBinary::Out fout(myfile,rank);
	String s = "Hello World!";
	fout.print(s);
	typename Vector<double>::Type m(10);

	srand48(3490201);
	for (SizeType i=0;i<m.size();i++)
		m[i]=drand48();
	fout.printVector(m,"MyVector");
	std::cout<<"MyVector\n";
	std::cout<<m;
	std::cout<<"------------------\n";

	PsimagLite::Matrix<float> a(10,20);
	for (SizeType i=0;i<a.n_row();i++)
		for (SizeType j=0;j<a.n_col();j++)
			a(i,j)=drand48();

	fout.write(a, "MyMatrix");
	fout.close();

	readMe(myfile);

}
