#include <PsimagLite/TextAndNumbersChecker.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

enum class ComparisonMode
{
	Generic,
	IgnoreLinePrefix,
	Observable
};

struct ParsedOption {
	ComparisonMode mode;
	std::string    value;
};

ParsedOption parseOption(const char* text)
{
	const std::string argument(text);
	const std::string ignoredLinePrefixOption = "--ignore-line-prefix=";
	const std::string insituLabelOption        = "--insitu-label=";

	if (argument.compare(0, ignoredLinePrefixOption.size(), ignoredLinePrefixOption)
	    == 0) {
		if (argument.size() == ignoredLinePrefixOption.size())
			throw std::invalid_argument("Invalid option: \"" + argument + "\"");
		return { ComparisonMode::IgnoreLinePrefix,
		         argument.substr(ignoredLinePrefixOption.size()) };
	}

	if (argument.compare(0, insituLabelOption.size(), insituLabelOption) == 0) {
		if (argument.size() == insituLabelOption.size())
			throw std::invalid_argument("Invalid option: \"" + argument + "\"");
		return { ComparisonMode::Observable,
		         argument.substr(insituLabelOption.size()) };
	}

	throw std::invalid_argument("Invalid option: \"" + argument + "\"");
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
		          << " FILE1 FILE2 TOLERANCE "
		             "[--ignore-line-prefix=PREFIX|--insitu-label=LABEL]\n";
		return EXIT_FAILURE;
	}

	try {
		const ParsedOption option = (argc == 5)
		    ? parseOption(argv[4])
		    : ParsedOption { ComparisonMode::Generic, "" };
		const std::string ignoredLinePrefix
		    = (option.mode == ComparisonMode::IgnoreLinePrefix) ? option.value : "";
		const PsimagLite::TextAndNumbersChecker checker(parseTolerance(argv[3]),
		                                                ignoredLinePrefix);

		if (option.mode == ComparisonMode::Observable)
			checker.runObservables(argv[1], argv[2], option.value);
		else
			checker.run(argv[1], argv[2]);
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return EXIT_FAILURE;
	}
}
