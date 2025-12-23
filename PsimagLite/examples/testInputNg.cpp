/* PSIDOC InputNg_Intro
  \section{Lesson I: InputNg}
  In this lesson we'll learn how to read a
  properly formatted input file. This way when
  you write your own C++ programs you'll have
  an already available way to read parameters
  from an input file and use them in your program.

  You can go ahead an compile this example with
\begin{tiny}
\begin{verbatim}
g++ testInputNg.cpp -std=c++11 -Wall -I ../src/ -I.. -DUSE_BOOST -L ../lib
-lpsimaglite -o testInputNg \end{verbatim} \end{tiny} Note that you need
boost-dev or boost-devel and also you must have lib/libsimaglite.a already
compiled. You can already go ahead an run it with the provided simple input
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
\begin{lstlisting}
PSIDOC_CONTINUE
*/

#include "InputCheckBase.h"
#include "InputNg.h"

/* PSIDOC_RESUME
\end{lstlisting}
InputNg has an option to check the inputs, so that
 you can define the labels that you expect to find
 in your program. This is optional but recommended.
 For this we create a class, say \texttt{MyInputCheck},
 and derive it from PsimagLite's \texttt{InputCheckBase}.
This inheritance isn't needed, but save us from
having to provide all functions, as the base class
implements defaults. For a more realistic use
case you can check DMRG++'s InputCheck.h under
dmrgpp/src/Engine. For now, here's our short
input checking class.
\begin{lstlisting}
PSIDOCCOPY InputNg_Class_MyInputCheck
\end{lstlisting}
In our example, we are defining a scalar
called myscalar, a vector called myvector, and
a string called mystring. This is what we
expect to read from the input file, even though
the writer of the input may add other labels.
*/

/* PSIDOC_CODE_START InputNg_Class_MyInputCheck nocapture */
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
/* PSIDOC_CODE_END */

/* PSIDOC InputNg_main1
Now for the actual reading of the input file, we'll use
\cppFile{InputNg}. We'll alias its type first with
\begin{verbatim}
typedef PsimagLite::InputNg<MyInputCheck> InputNgType;
\end{verbatim}
InputNg expects one template argument, our input checking class,
which we have just described.
We'll also have to use the actual name provided to this
program, which should be in \verb!argv[1]!, which we put
in C++ variable \texttt{filename}. We need to create an
object of our class MyInputCheck as well.
We then have the following code so far
\begin{lstlisting}
PSIDOC_CONTINUE
*/

int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::cerr << "USAGE " << argv[0] << " filename\n";
		return 1;
	}

	typedef PsimagLite::InputNg<MyInputCheck> InputNgType;

	std::string filename(argv[1]);
	MyInputCheck myInputCheck;

	/* PSIDOC_RESUME
	\end{lstlisting}
	We are now ready to read the file, which we do
	with the following two statements
	\begin{lstlisting}
	PSIDOC_CONTINUE
	*/

	InputNgType::Writeable ioWriteable(filename, myInputCheck);
	InputNgType::Readable io(ioWriteable);

	/* PSIDOC_RESUME
	\end{lstlisting}
	From now on, we can forget about the myInputCheck object,
	and the ioWriteable object as well, and consider
	only the io object, which we will use to read labels.
	The data from the file is now in memory, and the file
	does not have to be read or even present anymore.
	Let's now read some data from the file using the
	io object. (The data is now in memory, and it
	is \emph{not actually} read from the file, but we will use
	that terminology anyway.) We read the scalar first, and
	print it to the terminal with the following code.
	\begin{lstlisting}
	PSIDOC_CONTINUE
	*/

	int myscalar = 0;
	io.readline(myscalar, "myscalar=");

	std::cout << "I've read label myscalar with value ";
	std::cout << myscalar << " from " << io.filename() << "\n";

	/* PSIDOC_RESUME
	\end{lstlisting}
	The first argument to io.readline will be filled
	with the value from the file that follows the
	label myscalar. Even though the value will be filled,
	it's best practice to initialize it first anyway.

	Let's now read a vector. Note that we use just
	io.read to read vectors (and matrices), whereas we
	use io.readline to read scalars. The vector will
	be resized for you as needed, so you do not need
	to size it. If you choose to size it, then Ainur
	format will be able to use ellipsis to fill the vector
	and you may be able to enter in the input file something
	like myvector=[42, ...]; which will cause your vector
	to be filled with the number 42.
	\begin{lstlisting}
	PSIDOC_CONTINUE
	*/
	std::vector<double> v;
	io.read(v, "myvector");

	/* PSIDOC_RESUME
	\end{lstlisting}
	The two previous examples required the labels to
	be present in the input file. But what if we want
	to have an \emph{optional} label in the input file?
	To deal with that we put the io.readline statement
	inside a try and catch block, as follows.
	\begin{lstlisting}
	PSIDOCCOPY InputNg_code_main5
	\end{lstlisting}
	This way if the user provides the label mystring
	then the C++ variable mystr will have the value
	the user provided. Otherwise, the value of mystr
	will remain ``default'', and no error will be issued.
	*/

	/* PSIDOC_CODE_START InputNg_code_main5 nocapture */
	std::string mystr("default");
	try {
		io.readline(mystr, "mystring=");
	} catch (std::exception&) { }
}
/* PSIDOC_CODE_END */

/* PSIDOC InputNg_Recap

\subsection*{Recap}
In this lesson we have learned, blah blah.

*/
