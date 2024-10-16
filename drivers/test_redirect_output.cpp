#include "RedirectOutput.hh"

int main(int argc, char* argv[])
{
	Gep::RedirectOutput::setAppName(argv[0]);

	if (argc == 2) {
		Gep::RedirectOutput::doIt("test", std::ofstream::out, true);
	}
	std::cout << "Hello!\n";
	return 0;
}
