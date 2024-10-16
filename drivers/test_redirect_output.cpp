#include "RedirectOutput.hh"

int main(int argc, char* argv[])
{
	PsimagLite::RedirectOutput::setAppName(argv[0]);

	if (argc == 2) {
		PsimagLite::RedirectOutput::doIt("test", std::ofstream::out, true);
	}
	std::cout << "Hello!\n";
	return 0;
}
