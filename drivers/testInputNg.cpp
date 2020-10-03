/* PSIDOC InputNg
  \section{InputNg}
  blah blah blah
 */

#include "InputNg.h"
#include "InputCheckBase.h"

class MyInputCheck : public PsimagLite::InputCheckBase {

public:

	std::string import() const
	{
		std::string str("integer myscalar;\n");
		str += "vector myvector;\n";
		str += "string mystring;\n";
		return str;
	}
};

int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::cerr<<"USAGE "<<argv[0]<<" filename\n";
		return 1;
	}

	typedef PsimagLite::InputNg<MyInputCheck> InputNgType;

	std::string filename(argv[1]);
	MyInputCheck myInputCheck;
	InputNgType::Writeable ioWriteable(filename, myInputCheck);
	InputNgType::Readable io(ioWriteable);

	int myscalar = 0;
	io.readline(myscalar, "myscalar=");

	std::cout<<"I've read label myscalar with value ";
	std::cout<<myscalar<<" from "<<io.filename()<<"\n";
	std::vector<double> v;
	io.read(v, "myvector");

	std::string mystr;
	try {
		io.readline(mystr, "mystring=");
	} catch (std::exception&)
	{}
}

