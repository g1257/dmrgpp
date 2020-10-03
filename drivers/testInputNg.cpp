/* PSIDOC InputNg_Intro
  \section{Lesson I: InputNg}
  In this lesson we'll learn how to read a
  properly formatted input file. This way when
  you write your own C++ programs you'll have
  an already available way to read parameters
  from an input file and use them in your program.

  You can go ahead an compile this example with
\begin{footnotesize}
\begin{verbatim}
g++ testInputNg.cpp -I ../src/ -I.. -DUSE_BOOST -L ../lib -lpsimaglite -o testInputNg
\end{verbatim}
\end{footnotesize}
Note that you need boost-dev or boost-devel and also you
must have lib/PsimagLite already compiled.
You can already go ahead an run it with the provided simple input
\begin{verbatim}
./testInputNg testInput.ain
\end{verbatim}
What does this program do? It reads two or three ``labels'' from
testInput.ain and prints some of them to the terminal.
The labels have the form \texttt{label=value}, and we have
a scalar, a vector, and a string as examples.

So, now let's discuss the short program we have here.
First, note that we need two includes from PsimagLite.
PsimagLite files are usually under src, but in some
cases subdirectories are used. Here are the includes.
*/

/* PSIDOC_CODE_START InputNg_Includes nocapture */
#include "InputNg.h"
#include "InputCheckBase.h"
/* PSIDOC_CODE_END */

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

