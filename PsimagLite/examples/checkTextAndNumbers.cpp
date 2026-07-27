#include <PsimagLite/TextAndNumbersChecker.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

std::string parseIgnoredLinePrefix(const char* text)
{
	const std::string argument(text);
	const std::string option = "--ignore-line-prefix=";
	if (argument.compare(0, option.size(), option) != 0 || argument.size() == option.size())
		throw std::invalid_argument("Invalid option: \"" + argument + "\"");

	return argument.substr(option.size());
}

double parseTolerance(const char* text)
{
	const std::string token(text);
	std::size_t       parsed = 0;
	try {
		const double value = std::stod(token, &parsed);
		if (parsed != token.size() || !std::isfinite(value))
			throw std::invalid_argument("not a finite real number");
		return value;
	} catch (const std::invalid_argument&) {
		throw std::invalid_argument("Invalid tolerance: \"" + token + "\"");
	} catch (const std::out_of_range&) {
		throw std::invalid_argument("Tolerance is out of range: \"" + token + "\"");
	}
}

} // namespace

int main(int argc, char* argv[])
{
	if (argc != 4 && argc != 5) {
		std::cerr << "Usage: " << argv[0]
		          << " FILE1 FILE2 TOLERANCE [--ignore-line-prefix=PREFIX]\n";
		return EXIT_FAILURE;
	}

	std::string                                        file1;
	std::string                                        file2;
	std::unique_ptr<PsimagLite::TextAndNumbersChecker> checker;

	try {
		file1 = argv[1];
		file2 = argv[2];
		const std::string ignoredLinePrefix
		    = (argc == 5) ? parseIgnoredLinePrefix(argv[4]) : "";
		checker = std::make_unique<PsimagLite::TextAndNumbersChecker>(
		    parseTolerance(argv[3]), ignoredLinePrefix);
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return EXIT_FAILURE;
	}

	checker->run(file1, file2);
}
