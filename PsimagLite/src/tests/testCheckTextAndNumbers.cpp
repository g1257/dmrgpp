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
			const auto           stamp
			    = std::chrono::steady_clock::now().time_since_epoch().count();
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

		std::filesystem::path write(const std::string& name,
		                            const std::string& contents) const
		{
			const auto    path = directory_ / name;
			std::ofstream stream(path, std::ios::binary);
			if (!stream)
				throw std::runtime_error("Cannot create test file: "
				                         + path.string());

			stream.write(contents.data(),
			             static_cast<std::streamsize>(contents.size()));
			if (!stream)
				throw std::runtime_error("Cannot write test file: "
				                         + path.string());

			return path;
		}

		std::filesystem::path path(const std::string& name) const
		{
			return (directory_ / name);
		}

	private:

		std::filesystem::path directory_;
	};

	using Catch::Matchers::ContainsSubstring;
	using Catch::Matchers::MessageMatches;

} // namespace

TEST_CASE("TextAndNumbersChecker accepts empty files", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "");
	const auto            file2 = files.write("second.txt", "");

	CHECK_NOTHROW(checker.run(file1, file2));
}

TEST_CASE("TextAndNumbersChecker compares words verbatim", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
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
		CHECK_THROWS_WITH(
		    checker.run(file1, file2),
		    "Token 1: words differ; first file has \"Word\", second file has \"word\"");
	}

	SECTION("punctuation differences fail")
	{
		const auto file1 = files.write("first.txt", "word,");
		const auto file2 = files.write("second.txt", "word");
		CHECK_THROWS_WITH(
		    checker.run(file1, file2),
		    "Token 1: words differ; first file has \"word,\", second file has \"word\"");
	}
}

TEST_CASE("TextAndNumbersChecker recognizes every ASCII whitespace", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const auto file1 = files.write("first.txt", "one two\tthree\nfour\rfive\fsix\vseven");
	const auto file2 = files.write("second.txt", "one\ttwo three\rfour\nfive\vsix\fseven");

	CHECK_NOTHROW(checker.run(file1, file2));
}

TEST_CASE("TextAndNumbersChecker compares parsed integers", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
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
		CHECK_THROWS_WITH(
		    checker.run(file1, file2),
		    "Token 2: integer values differ; first file has \"3\", second file has \"4\"");
	}
}

TEST_CASE("TextAndNumbersChecker rejects different token classes", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "3");
	const auto            file2 = files.write("second.txt", "3.0");

	CHECK_THROWS_WITH(checker.run(file1, file2),
	                  "Token 1: different token classes (integer and double); first file has "
	                  "\"3\", second file has \"3.0\"");
}

TEST_CASE("TextAndNumbersChecker recognizes strict double syntax", "[TextAndNumbersChecker]")
{
	TemporaryFiles                 files;
	TextAndNumbersChecker          checker(0.01);
	const std::vector<std::string> doubles { "3.0", ".5", "5.", "1e3", "-2.5e-4" };

	for (const auto& token : doubles) {
		CAPTURE(token);
		const auto file1 = files.write("first.txt", token);
		const auto file2 = files.write("second.txt", token);
		CHECK_NOTHROW(checker.run(file1, file2));
	}
}

TEST_CASE("TextAndNumbersChecker treats malformed numeric tokens as words",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
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
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);

	SECTION("first file is longer")
	{
		const auto file1 = files.write("first.txt", "one two");
		const auto file2 = files.write("second.txt", "one");
		CHECK_THROWS_MATCHES(
		    checker.run(file1, file2),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring(
		        "Token 2: the files contain different numbers of tokens")));
	}

	SECTION("second file is longer")
	{
		const auto file1 = files.write("first.txt", "one");
		const auto file2 = files.write("second.txt", "one two");
		CHECK_THROWS_MATCHES(
		    checker.run(file1, file2),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring(
		        "Token 2: the files contain different numbers of tokens")));
	}
}

TEST_CASE("TextAndNumbersChecker stops at the first mismatch", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "same first third");
	const auto            file2 = files.write("second.txt", "same second fourth");

	CHECK_THROWS_WITH(
	    checker.run(file1, file2),
	    "Token 2: words differ; first file has \"first\", second file has \"second\"");
}

TEST_CASE("TextAndNumbersChecker reports files that cannot be opened", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);

	SECTION("first file is missing")
	{
		const auto missing = files.path("missing-first.txt");
		const auto file2   = files.write("second.txt", "same");
		CHECK_THROWS_WITH(checker.run(missing, file2),
		                  "Cannot open first file: " + missing.string());
	}

	SECTION("second file is missing")
	{
		const auto file1   = files.write("first.txt", "same");
		const auto missing = files.path("missing-second.txt");
		CHECK_THROWS_WITH(checker.run(file1, missing),
		                  "Cannot open second file: " + missing.string());
	}
}

TEST_CASE("TextAndNumbersChecker reports generic stream read errors", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const auto            firstDirectory  = files.path("first-directory");
	const auto            secondDirectory = files.path("second-directory");
	std::filesystem::create_directory(firstDirectory);
	std::filesystem::create_directory(secondDirectory);

	CHECK_THROWS_MATCHES(checker.run(firstDirectory, secondDirectory),
	                     std::runtime_error,
	                     MessageMatches(ContainsSubstring("Error while reading first file")));
}

TEST_CASE("TextAndNumbersChecker reports integer overflow", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const std::string     value = "999999999999999999999999999999999999";
	const auto            file1 = files.write("first.txt", value);
	const auto            file2 = files.write("second.txt", value);

	CHECK_THROWS_WITH(checker.run(file1, file2),
	                  "Integer token is out of range: \"" + value + "\"");
}

TEST_CASE("TextAndNumbersChecker rejects an oversized generic integer token",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const std::string     value(200000, '9');
	const auto            file1 = files.write("first.txt", value);
	const auto            file2 = files.write("second.txt", value);

	CHECK_THROWS_MATCHES(checker.run(file1, file2),
	                     std::runtime_error,
	                     MessageMatches(ContainsSubstring("Integer token is out of range")));
}

TEST_CASE("TextAndNumbersChecker reports double overflow", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.01);
	const auto            file1 = files.write("first.txt", "1e9999");
	const auto            file2 = files.write("second.txt", "1e9999");

	CHECK_THROWS_WITH(checker.run(file1, file2), "Double token is out of range: \"1e9999\"");
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

TEST_CASE("TextAndNumbersChecker rejects negative tolerances", "[TextAndNumbersChecker]")
{
	CHECK_THROWS_AS(TextAndNumbersChecker(-0.01), std::invalid_argument);
}

TEST_CASE("TextAndNumbersChecker ignores configured line prefixes", "[TextAndNumbersChecker]")
{
	TemporaryFiles files;

	SECTION("marked lines are ignored independently")
	{
		TextAndNumbersChecker checker(0.001, "##");
		const auto            file1 = files.write(
                    "first.txt", "## generated at 1.25 seconds\nword\n## first only\n1.0\n");
		const auto file2
		    = files.write("second.txt", "## generated at 9.75 seconds\nword\n1.0005\n");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("the marker must start at column one")
	{
		TextAndNumbersChecker checker(0.001, "##");
		const auto            file1 = files.write("first.txt", " ## not ignored\nword\n");
		const auto            file2 = files.write("second.txt", "word\n");
		CHECK_THROWS(checker.run(file1, file2));
	}

	SECTION("a marker in the middle of a line is not ignored")
	{
		TextAndNumbersChecker checker(0.001, "##");
		const auto            file1 = files.write("first.txt", "word ## first\n");
		const auto            file2 = files.write("second.txt", "word ## second\n");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring("Token 3: words differ")));
	}

	SECTION("strict comparison remains the default")
	{
		TextAndNumbersChecker checker(0.001);
		const auto            file1 = files.write("first.txt", "## first\nword\n");
		const auto            file2 = files.write("second.txt", "## second\nword\n");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring("Token 2: words differ")));
	}
}

TEST_CASE("TextAndNumbersChecker compares real numbers with absolute tolerance",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles files;

	SECTION("equal values match")
	{
		TextAndNumbersChecker checker(0.0);
		const auto            file1 = files.write("first.txt", "1.0");
		const auto            file2 = files.write("second.txt", "1.0");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("a difference smaller than the tolerance matches")
	{
		TextAndNumbersChecker checker(0.125);
		const auto            file1 = files.write("first.txt", "-1.0");
		const auto            file2 = files.write("second.txt", "-1.0625");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("a difference equal to the tolerance matches")
	{
		TextAndNumbersChecker checker(0.125);
		const auto            file1 = files.write("first.txt", "1.0");
		const auto            file2 = files.write("second.txt", "1.125");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("a difference larger than the tolerance fails")
	{
		TextAndNumbersChecker checker(0.125);
		const auto            file1 = files.write("first.txt", "1.0");
		const auto            file2 = files.write("second.txt", "1.25");
		CHECK_THROWS_WITH(checker.run(file1, file2),
		                  "Token 1: double values differ; first file has \"1.0\", second "
		                  "file has \"1.25\"");
	}

	SECTION("exponent notation uses parsed values")
	{
		TextAndNumbersChecker checker(0.125);
		const auto            file1 = files.write("first.txt", "1e3");
		const auto            file2 = files.write("second.txt", "1000.125");
		CHECK_NOTHROW(checker.run(file1, file2));
	}
}

TEST_CASE("TextAndNumbersChecker compares complex numbers by magnitude", "[TextAndNumbersChecker]")
{
	TemporaryFiles files;

	SECTION("equivalent representations match")
	{
		TextAndNumbersChecker checker(0.0);
		const auto file1 = files.write("first.txt", "(1.0,-2.0) (1e0,-2e0) (1,0)");
		const auto file2 = files.write("second.txt", "(1.00,-2.00) (1.0,-2.0) (1.0,0.0)");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("a difference smaller than the tolerance matches")
	{
		TextAndNumbersChecker checker(0.6);
		const auto            file1 = files.write("first.txt", "(0,0)");
		const auto            file2 = files.write("second.txt", "(0.3,0.4)");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("a difference equal to the tolerance matches")
	{
		TextAndNumbersChecker checker(5.0);
		const auto            file1 = files.write("first.txt", "(0,0)");
		const auto            file2 = files.write("second.txt", "(3,4)");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("a difference larger than the tolerance fails")
	{
		TextAndNumbersChecker checker(0.49);
		const auto            file1 = files.write("first.txt", "(0,0)");
		const auto            file2 = files.write("second.txt", "(0.3,0.4)");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring("complex values differ")));
	}

	SECTION("differences in both components use one magnitude")
	{
		TextAndNumbersChecker checker(4.0);
		const auto            file1 = files.write("first.txt", "(0,0)");
		const auto            file2 = files.write("second.txt", "(3,3)");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring("complex values differ")));
	}
}

TEST_CASE("TextAndNumbersChecker classifies complex tokens strictly", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.0);

	SECTION("identical malformed tokens remain words")
	{
		const auto file1 = files.write("first.txt", "(1,) (1,2)suffix (1, 2)");
		const auto file2 = files.write("second.txt", "(1,) (1,2)suffix (1, 2)");
		CHECK_NOTHROW(checker.run(file1, file2));
	}

	SECTION("different malformed tokens fail as words")
	{
		const auto file1 = files.write("first.txt", "(1,)");
		const auto file2 = files.write("second.txt", "(2,)");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring("words differ")));
	}

	SECTION("a complex token and a word have different classes")
	{
		const auto file1 = files.write("first.txt", "(1,2)");
		const auto file2 = files.write("second.txt", "word");
		CHECK_THROWS_MATCHES(checker.run(file1, file2),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "different token classes (complex and word)")));
	}
}

TEST_CASE("TextAndNumbersChecker rejects invalid complex components", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.0);

	SECTION("overflow")
	{
		const auto file1 = files.write("first.txt", "(1e9999,0)");
		const auto file2 = files.write("second.txt", "(1e9999,0)");
		CHECK_THROWS_MATCHES(
		    checker.run(file1, file2),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring("Complex component is out of range")));
	}

	SECTION("non-finite values")
	{
		const auto file1 = files.write("first.txt", "(nan,0)");
		const auto file2 = files.write("second.txt", "(nan,0)");
		CHECK_THROWS_MATCHES(
		    checker.run(file1, file2),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring("Complex token has a non-finite component")));
	}

	SECTION("oversized components fail without exhausting the parser")
	{
		const std::string value = "(" + std::string(200000, '9') + ",0)";
		const auto        file1 = files.write("first.txt", value);
		const auto        file2 = files.write("second.txt", value);
		CHECK_THROWS_MATCHES(
		    checker.run(file1, file2),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring("Complex component is out of range")));
	}
}

TEST_CASE("TextAndNumbersChecker compares final indexed observable records",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(1e-6);
	const std::string     label     = "<gs|sz|P0>";
	const auto            reference = files.write("reference.txt",
                                           "DMRG++ reference log\n"
	                                              "5 (9,9) 0.300000 <gs|sz|P0> (8,8)\n"
	                                              "4 (0.25,0) 0.3 <gs|sz|P0> (0,0)\n"
	                                              "5 (-0.18,0.05) 0.300000 <gs|sz|P0> (0,0)\n");
	const auto            actual    = files.write("actual.txt",
                                        "PsiApp: CmdLine: dmrg -f input21.ain <gs|sz|P0>\n"
	                                              "5 (-7,-7) 0.3 <gs|sz|P0> (-6,-6)\n"
	                                              "5 (-0.18,0.05) 0.3 <gs|sz|P0> (0,0)\n"
	                                              "4 (0.250000,0.000000) 0.300000 <gs|sz|P0> (0.0,0.0)\n");

	CHECK_NOTHROW(checker.runObservables(reference, actual, label));
}

TEST_CASE("TextAndNumbersChecker compares both final observable data fields",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(0.11);
	const std::string     label = "<gs|sz|P0>";
	const auto reference = files.write("reference.txt", "5 (1,2) 0.3 <gs|sz|P0> (3,4)\n");

	SECTION("early agreement cannot hide a final value mismatch")
	{
		const auto actual = files.write("actual.txt",
		                                "5 (1,2) 0.3 <gs|sz|P0> (3,4)\n"
		                                "5 (1.2,2) 0.3 <gs|sz|P0> (3,4)\n");
		CHECK_THROWS_MATCHES(
		    checker.runObservables(reference, actual, label),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring(
		        "observable value mismatch at (site=5, time=0.3, label=<gs|sz|P0>)")));
	}

	SECTION("superdensity is compared independently")
	{
		const auto actual = files.write("actual.txt", "5 (1,2) 0.3 <gs|sz|P0> (3.2,4)\n");
		CHECK_THROWS_MATCHES(
		    checker.runObservables(reference, actual, label),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring(
		        "superdensity mismatch at (site=5, time=0.3, label=<gs|sz|P0>)")));
	}

	SECTION("complex magnitude tolerance applies to both fields")
	{
		const auto actual
		    = files.write("actual.txt", "5 (1.06,2.08) 0.3 <gs|sz|P0> (3.06,4.08)\n");
		CHECK_NOTHROW(checker.runObservables(reference, actual, label));
	}
}

TEST_CASE("TextAndNumbersChecker requires exact observable index sets", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(1.0);
	const std::string     label     = "<gs|sz|P0>";
	const auto            reference = files.write("reference.txt",
                                           "4 (1,0) 0.3 <gs|sz|P0> (0,0)\n"
	                                              "5 (2,0) 0.3 <gs|sz|P0> (0,0)\n");

	SECTION("a genuinely different time is a key mismatch")
	{
		const auto actual = files.write("actual.txt",
		                                "4 (1,0) 0.31 <gs|sz|P0> (0,0)\n"
		                                "5 (2,0) 0.3 <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(
		    checker.runObservables(reference, actual, label),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring("site=4, time=0.3, label=<gs|sz|P0>")));
	}

	SECTION("a missing index reports its full identity")
	{
		const auto actual = files.write("actual.txt", "4 (1,0) 0.3 <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(
		    checker.runObservables(reference, actual, label),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring(
		        "missing from actual file: (site=5, time=0.3, label=<gs|sz|P0>)")));
	}

	SECTION("an extra index reports its full identity")
	{
		const auto actual = files.write("actual.txt",
		                                "4 (1,0) 0.3 <gs|sz|P0> (0,0)\n"
		                                "5 (2,0) 0.3 <gs|sz|P0> (0,0)\n"
		                                "6 (3,0) 0.3 <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(
		    checker.runObservables(reference, actual, label),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring(
		        "extra in actual file: (site=6, time=0.3, label=<gs|sz|P0>)")));
	}
}

TEST_CASE("TextAndNumbersChecker validates selected observable records", "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(1e-6);
	const std::string     label = "<gs|sz|P0>";
	const auto reference = files.write("reference.txt", "5 (1,2) 0.3 <gs|sz|P0> (0,0)\n");

	SECTION("a malformed selected record fails")
	{
		const auto actual = files.write("actual.txt", "5 malformed 0.3 <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(checker.runObservables(reference, actual, label),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Malformed observable record in actual file at line 1")));
	}

	SECTION("a malformed selected record with no numeric site fails closed")
	{
		const auto actual = files.write("actual.txt",
		                                "x bad 0.3 <gs|sz|P0> bad\n"
		                                "5 (1,2) 0.3 <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(checker.runObservables(reference, actual, label),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Malformed observable record in actual file at line 1")));
	}

	SECTION("a missing value cannot shift the selected label out of position")
	{
		const auto actual = files.write("actual.txt",
		                                "5 0.3 <gs|sz|P0> (0,0)\n"
		                                "5 (1,2) 0.3 <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(checker.runObservables(reference, actual, label),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Malformed observable record in actual file at line 1")));
	}

	SECTION("a missing time cannot shift the selected label out of position")
	{
		const auto actual = files.write("actual.txt",
		                                "5 (1,2) <gs|sz|P0> (0,0)\n"
		                                "5 (1,2) 0.3 <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(checker.runObservables(reference, actual, label),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Malformed observable record in actual file at line 1")));
	}

	SECTION("trailing fields fail")
	{
		const auto actual
		    = files.write("actual.txt", "5 (1,2) 0.3 <gs|sz|P0> (0,0) trailing\n");
		CHECK_THROWS_MATCHES(checker.runObservables(reference, actual, label),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Malformed observable record in actual file at line 1")));
	}

	SECTION("no selected records fails explicitly")
	{
		const auto actual = files.write("actual.txt",
		                                "PsiApp: CmdLine: dmrg -f input21.ain <gs|sz|P0>\n"
		                                "5 (1,2) 0.3 <gs|other|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(
		    checker.runObservables(reference, actual, label),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring(
		        "No observable records for label <gs|sz|P0> in actual file")));
	}

	SECTION("non-finite time fails")
	{
		const auto actual = files.write("actual.txt", "5 (1,2) nan <gs|sz|P0> (0,0)\n");
		CHECK_THROWS_MATCHES(checker.runObservables(reference, actual, label),
		                     std::runtime_error,
		                     MessageMatches(ContainsSubstring(
		                         "Malformed observable record in actual file at line 1")));
	}

	SECTION("a read error is not treated as end of file")
	{
		const auto actual = files.path("actual-directory");
		std::filesystem::create_directory(actual);
		CHECK_THROWS_MATCHES(
		    checker.runObservables(reference, actual, label),
		    std::runtime_error,
		    MessageMatches(ContainsSubstring("Error while reading actual file")));
	}
}

TEST_CASE("TextAndNumbersChecker rejects an oversized observable site token",
          "[TextAndNumbersChecker]")
{
	TemporaryFiles        files;
	TextAndNumbersChecker checker(1e-6);
	const std::string     label = "<gs|sz|P0>";
	const auto reference = files.write("reference.txt", "5 (1,2) 0.3 <gs|sz|P0> (0,0)\n");
	const auto actual
	    = files.write("actual.txt", std::string(200000, '9') + " (1,2) 0.3 <gs|sz|P0> (0,0)\n");

	CHECK_THROWS_MATCHES(checker.runObservables(reference, actual, label),
	                     std::runtime_error,
	                     MessageMatches(ContainsSubstring("Site token is out of range")));
}

} // namespace PsimagLite
