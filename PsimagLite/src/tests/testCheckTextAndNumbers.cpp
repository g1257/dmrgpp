#include <PsimagLite/TextAndNumbersChecker.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace PsimagLite {
namespace {

class TemporaryFiles {
public:
	TemporaryFiles()
	{
		static unsigned long counter = 0;
		const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
		directory_ = std::filesystem::temp_directory_path()
		    / ("psimaglite-text-checker-" + std::to_string(stamp) + "-"
		       + std::to_string(counter++));
		std::filesystem::create_directory(directory_);
	}

	~TemporaryFiles()
	{
		std::error_code error;
		std::filesystem::remove_all(directory_, error);
	}

	std::string write(const std::string& name, const std::string& contents) const
	{
		const auto path = directory_ / name;
		std::ofstream stream(path, std::ios::binary);
		if (!stream)
			throw std::runtime_error("Cannot create test file: " + path.string());

		stream.write(contents.data(), static_cast<std::streamsize>(contents.size()));
		if (!stream)
			throw std::runtime_error("Cannot write test file: " + path.string());

		return path.string();
	}

	std::string path(const std::string& name) const { return (directory_ / name).string(); }

private:
	std::filesystem::path directory_;
};

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;

} // namespace

TEST_CASE("TextAndNumbersChecker accepts empty files", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "");
	const auto            file2 = files.write("second.txt", "");

	CHECK_NOTHROW(checker.run(file1, file2));
}

TEST_CASE("TextAndNumbersChecker compares words verbatim", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);

	SECTION("identical words match")
	{
		const auto file1 = files.write("first.txt", "Alpha beta word,");
		const auto file2 = files.write("second.txt", "Alpha beta word,");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("case differences fail")
	{
		const auto file1 = files.write("first.txt", "Word");
		const auto file2 = files.write("second.txt", "word");
		CHECK_THROWS_WITH(checker.run(file1, file2),
		                  "Token 1: words differ; first file has \"Word\", second file has \"word\"");
	}

	SECTION("punctuation differences fail")
	{
		const auto file1 = files.write("first.txt", "word,");
		const auto file2 = files.write("second.txt", "word");
		CHECK_THROWS_WITH(checker.run(file1, file2),
		                  "Token 1: words differ; first file has \"word,\", second file has \"word\"");
	}
}

TEST_CASE("TextAndNumbersChecker recognizes every ASCII whitespace", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const auto file1 = files.write("first.txt", "one two\tthree\nfour\rfive\fsix\vseven");
	const auto file2 = files.write("second.txt", "one\ttwo three\rfour\nfive\vsix\fseven");

	CHECK_NOTHROW(checker.run(file1, file2));
}

TEST_CASE("TextAndNumbersChecker compares parsed integers", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);

	SECTION("equivalent representations match")
	{
		const auto file1 = files.write("first.txt", "03 +3 -0 -0042");
		const auto file2 = files.write("second.txt", "3 3 0 -42");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("different values fail")
	{
		const auto file1 = files.write("first.txt", "same 3");
		const auto file2 = files.write("second.txt", "same 4");
		CHECK_THROWS_WITH(checker.run(file1, file2),
		                  "Token 2: integer values differ; first file has \"3\", second file has \"4\"");
	}
}

TEST_CASE("TextAndNumbersChecker rejects different token classes", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "3");
	const auto            file2 = files.write("second.txt", "3.0");

	CHECK_THROWS_WITH(
	    checker.run(file1, file2),
	    "Token 1: different token classes (integer and double); first file has \"3\", second file has \"3.0\"");
}

TEST_CASE("TextAndNumbersChecker recognizes strict double syntax", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const std::vector<std::string> doubles { "3.0", ".5", "5.", "1e3", "-2.5e-4" };

	for (const auto& token : doubles) {
		CAPTURE(token);
		const auto file1 = files.write("first.txt", token);
		const auto file2 = files.write("second.txt", token);
		CHECK_THROWS_WITH(checker.run(file1, file2), "checkRealNumbers is not implemented");
	}
}

TEST_CASE("TextAndNumbersChecker treats malformed numeric tokens as words",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);

	SECTION("identical tokens match")
	{
		const auto file1 = files.write("first.txt", "1e . + 3.0suffix");
		const auto file2 = files.write("second.txt", "1e . + 3.0suffix");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("different tokens fail as words")
	{
		const auto file1 = files.write("first.txt", "1e . + 3.0suffix");
		const auto file2 = files.write("second.txt", "1e . + 4.0suffix");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring("Token 4: words differ")));
	}
}

TEST_CASE("TextAndNumbersChecker detects unequal token counts", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);

	SECTION("first file is longer")
	{
		const auto file1 = files.write("first.txt", "one two");
		const auto file2 = files.write("second.txt", "one");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Token 2: the files contain different numbers of tokens")));
	}

	SECTION("second file is longer")
	{
		const auto file1 = files.write("first.txt", "one");
		const auto file2 = files.write("second.txt", "one two");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Token 2: the files contain different numbers of tokens")));
	}
}

TEST_CASE("TextAndNumbersChecker stops at the first mismatch", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "same first third");
	const auto            file2 = files.write("second.txt", "same second fourth");

	CHECK_THROWS_WITH(
	    checker.run(file1, file2),
	    "Token 2: words differ; first file has \"first\", second file has \"second\"");
}

TEST_CASE("TextAndNumbersChecker reports files that cannot be opened", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);

	SECTION("first file is missing")
	{
		const auto missing = files.path("missing-first.txt");
		const auto file2   = files.write("second.txt", "same");
		CHECK_THROWS_WITH(checker.run(missing, file2), "Cannot open first file: " + missing);
	}

	SECTION("second file is missing")
	{
		const auto file1   = files.write("first.txt", "same");
		const auto missing = files.path("missing-second.txt");
		CHECK_THROWS_WITH(checker.run(file1, missing), "Cannot open second file: " + missing);
	}
}

TEST_CASE("TextAndNumbersChecker reports integer overflow", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const std::string     value = "999999999999999999999999999999999999";
	const auto            file1 = files.write("first.txt", value);
	const auto            file2 = files.write("second.txt", value);

	CHECK_THROWS_WITH(checker.run(file1, file2),
	                  "Integer token is out of range: \"" + value + "\"");
}

TEST_CASE("TextAndNumbersChecker reports double overflow", "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "1e9999");
	const auto            file2 = files.write("second.txt", "1e9999");

	CHECK_THROWS_WITH(checker.run(file1, file2),
	                  "Double token is out of range: \"1e9999\"");
}

TEST_CASE("TextAndNumbersChecker rejects non-finite tolerances", "[TextAndNumbersChecker]")
{
	CHECK_THROWS_AS(TextAndNumbersChecker(std::numeric_limits<double>::infinity()),
	                std::invalid_argument);
	CHECK_THROWS_AS(TextAndNumbersChecker(-std::numeric_limits<double>::infinity()),
	                std::invalid_argument);
	CHECK_THROWS_AS(TextAndNumbersChecker(std::numeric_limits<double>::quiet_NaN()),
	                std::invalid_argument);
}

TEST_CASE("TextAndNumbersChecker real-number comparison is unimplemented",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles       files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "1.0");
	const auto            file2 = files.write("second.txt", "1.0");

	CHECK_THROWS_WITH(checker.run(file1, file2), "checkRealNumbers is not implemented");
}

} // namespace PsimagLite
