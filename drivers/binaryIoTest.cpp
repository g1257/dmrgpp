#include "IoBinary.h"
#include "Vector.h"

void readMe(const std::string& myfile)
{
	PsimagLite::IoBinary::In fin(myfile);
	std::vector<double> v;
	fin.read(v,"MyVector");
	std::cout<<"MyVector\n";
	std::cout<<v;
	std::cout<<"------------------\n";

	PsimagLite::Matrix<float> m;
	fin.readMatrix(m,"MyMatrix");
	std::cout<<"MyMatrix";
	std::cout<<m;
	std::cout<<"------------------\n";
}

int main()
{
	size_t rank = 0;
	std::string myfile = "myfile.txt";
	PsimagLite::IoBinary::Out fout(myfile,rank);
	std::string s = "Hello World!";
	fout.print(s);
	std::vector<double> m(10);

	srand48(3490201);
	for (size_t i=0;i<m.size();i++)
		m[i]=drand48();
	fout.printVector(m,"MyVector");
	std::cout<<"MyVector\n";
	std::cout<<m;
	std::cout<<"------------------\n";

	PsimagLite::Matrix<float> a(10,20);
	for (size_t i=0;i<a.n_row();i++)
		for (size_t j=0;j<a.n_col();j++)
			a(i,j)=drand48();

	fout.printMatrix(a,"MyMatrix");
	fout.close();

	readMe(myfile);

}
