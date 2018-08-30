#include "Io/IoNg.h"

void createVector(std::vector<bool>& v)
{
	SizeType n = v.size();
	for (SizeType i = 0; i < n; ++i)
		v[i] = (drand48() < 0.5);
}

void saveAll(PsimagLite::String filename,
             const std::vector<bool>& v1,
             PsimagLite::String vname1,
             const std::vector<bool>& v2,
             PsimagLite::String vname2)
{
	PsimagLite::IoNg::Out io(filename, PsimagLite::IoNg::ACC_TRUNC);
	io.write(v1, vname1);
	io.write(v2, vname2);
}

void compare(const std::vector<bool>& v,
             const std::vector<bool>& w)
{
	SizeType n = v.size();
	if (n != w.size())
		std::cerr<<"DIFFER on SIZE!!!\n";

	for (SizeType i = 0; i < n; ++i) {
		if (v[i] == w[i]) continue;
		std::cerr<<"DIFFER!!!\n";
	}

	std::cerr<<"OK\n";

}
int main(int argc, char **argv)
{
	if (argc < 2) return 1;
	SizeType n = atoi(argv[1]);
	std::vector<bool> v1(n, false);
	createVector(v1);
	std::vector<bool> v2(n, false);
	createVector(v2);
	PsimagLite::String filename = "file.hd5";
	PsimagLite::String vname1 = "myv1";
	PsimagLite::String vname2 = "myv2";
	saveAll(filename, v1, vname1, v2, vname2);

	std::vector<bool> w1;
	std::vector<bool> w2;
	PsimagLite::IoNg::In io(filename);
	io.read(w1, vname1);
	io.read(w2, vname2);
	compare(v1, w1);
	compare(v2, w2);
}


