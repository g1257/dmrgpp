#include <PsimagLite/TextAndNumbersChecker.hpp>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

std::string file1;
std::string file2;
std::unique_ptr<PsimagLite::TextAndNumbersChecker> checker;

double parseTolerance(const char* text)
{
	const std::string token(text);
	std::size_t parsed = 0;
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

TEST_CASE("text and numbers in two files match", "[TextAndNumbersChecker]")
{
	REQUIRE(static_cast<bool>(checker));
	REQUIRE_NOTHROW(checker->run(file1, file2));
}

int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cerr << "Usage: " << argv[0] << " FILE1 FILE2 TOLERANCE\n";
		return EXIT_FAILURE;
	}

	try {
		file1 = argv[1];
		file2 = argv[2];
		checker = std::make_unique<PsimagLite::TextAndNumbersChecker>(parseTolerance(argv[3]));
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return EXIT_FAILURE;
	}

	char* catchArguments[] = { argv[0] };
	return Catch::Session().run(1, catchArguments);
}
