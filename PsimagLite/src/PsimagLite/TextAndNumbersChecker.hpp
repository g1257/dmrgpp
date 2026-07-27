//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   PsimagLite/TextAndNumbersChecker.hpp
 * \brief  Compare text files containing words and numbers
 */
//---------------------------------------------------------------------------//

#ifndef TEXT_AND_NUMBERS_CHECKER_HPP
#define TEXT_AND_NUMBERS_CHECKER_HPP

#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace PsimagLite {

//===========================================================================//
/*!
 * \class TextAndNumbersChecker
 *
 * \brief Compare whitespace-separated words and numbers in two text files
 *
 * Words are compared verbatim, integers by value, and real numbers using an
 * absolute tolerance. Tokens in corresponding positions must have the same
 * class. Lines beginning with an optional prefix can be excluded from the
 * comparison.
 */
//===========================================================================//
class TextAndNumbersChecker {
public:

	using RealType = double;

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Constructor
	 *
	 * \param[in] tolerance          Maximum absolute difference allowed between
	 *                               corresponding real-number tokens
	 * \param[in] ignoredLinePrefix  Prefix identifying lines to omit; an empty
	 *                               prefix enables strict comparison of all lines
	 *
	 * 	hrows std::invalid_argument if the tolerance is negative or non-finite
	 */
	explicit TextAndNumbersChecker(RealType           tolerance,
	                               const std::string& ignoredLinePrefix = "")
	    : tolerance_(tolerance)
	    , ignoredLinePrefix_(ignoredLinePrefix)
	{
		if (!std::isfinite(tolerance_) || tolerance_ < 0)
			throw std::invalid_argument(
			    "The tolerance must be non-negative and finite");
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Compare two text files token by token
	 *
	 * \param[in] file1 Path to the first file
	 * \param[in] file2 Path to the second file
	 *
	 * 	hrows std::runtime_error if a file cannot be opened, a numeric token
	 *         cannot be represented, or the file contents do not match
	 */
	void run(const std::string& file1, const std::string& file2) const
	{
		std::ifstream stream1(file1);
		if (!stream1)
			throw std::runtime_error("Cannot open first file: " + file1);

		std::ifstream stream2(file2);
		if (!stream2)
			throw std::runtime_error("Cannot open second file: " + file2);

		TokenReader reader1(stream1, ignoredLinePrefix_);
		TokenReader reader2(stream2, ignoredLinePrefix_);
		std::size_t position = 1;
		for (;;) {
			std::string token1;
			std::string token2;
			const bool  hasToken1 = reader1.read(token1);
			const bool  hasToken2 = reader2.read(token2);

			if (!hasToken1 && !hasToken2)
				return;

			if (hasToken1 != hasToken2) {
				const std::string shown1
				    = hasToken1 ? quote(token1) : "<end of file>";
				const std::string shown2
				    = hasToken2 ? quote(token2) : "<end of file>";
				fail(position,
				     shown1,
				     shown2,
				     "the files contain different numbers of tokens");
			}

			compareTokens(token1, token2, position);
			++position;
		}
	}

private:

	enum class TokenClass
	{
		Word,
		Integer,
		Double
	};

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Read the next ASCII-whitespace-delimited token from a stream
	 *
	 * \param[in/out] stream Stream from which to read
	 * \param[out] token     Token read from the stream
	 *
	 * \returns True if a token was read, or false at end of stream
	 */
	static bool readToken(std::istream& stream, std::string& token)
	{
		token.clear();
		char c = 0;
		while (stream.get(c)) {
			if (!std::isspace(c)) {
				token.push_back(c);
				break;
			}
		}

		if (token.empty())
			return false;

		while (stream.get(c)) {
			if (std::isspace(c))
				break;
			token.push_back(c);
		}

		return true;
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \class TokenReader
	 *
	 * \brief Read tokens while omitting lines with a configured prefix
	 */
	class TokenReader {
	public:

		//---------------------------------------------------------------------------//
		/*!
		 * \brief Constructor
		 *
		 * \param[in/out] stream           Stream from which to read
		 * \param[in] ignoredLinePrefix    Prefix identifying lines to omit
		 */
		TokenReader(std::istream& stream, const std::string& ignoredLinePrefix)
		    : stream_(stream)
		    , ignoredLinePrefix_(ignoredLinePrefix)
		{ }

		//---------------------------------------------------------------------------//
		/*!
		 * \brief Read the next token from a non-ignored line
		 *
		 * \param[out] token Token read from the stream
		 *
		 * \returns True if a token was read, or false at end of stream
		 */
		bool read(std::string& token)
		{
			for (;;) {
				if (readToken(lineStream_, token))
					return true;

				if (!std::getline(stream_, line_))
					return false;

				if (!ignoredLinePrefix_.empty()
				    && line_.compare(
				           0, ignoredLinePrefix_.size(), ignoredLinePrefix_)
				        == 0)
					continue;

				lineStream_.clear();
				lineStream_.str(line_);
			}
		}

	private:

		std::istream&      stream_;
		const std::string& ignoredLinePrefix_;
		std::istringstream lineStream_;
		std::string        line_;
	};

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Classify a token as a word, integer, or real number
	 *
	 * \param[in] token Token to classify
	 *
	 * \returns The token class
	 */
	static TokenClass classify(const std::string& token)
	{
		static const std::regex integerPattern(R"([+-]?[0-9]+)");
		static const std::regex doublePattern(
		    R"([+-]?(([0-9]+\.[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?|[0-9]+[eE][+-]?[0-9]+))");

		if (std::regex_match(token, integerPattern))
			return TokenClass::Integer;
		if (std::regex_match(token, doublePattern))
			return TokenClass::Double;
		return TokenClass::Word;
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Return a printable name for a token class
	 *
	 * \param[in] tokenClass Token class to name
	 *
	 * \returns A static string naming the token class
	 */
	static const char* className(TokenClass tokenClass)
	{
		switch (tokenClass) {
		case TokenClass::Word:
			return "word";
		case TokenClass::Integer:
			return "integer";
		case TokenClass::Double:
			return "double";
		}

		return "unknown";
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Parse an integer token
	 *
	 * \param[in] token Token to parse
	 *
	 * \returns The represented integer value
	 * 	hrows std::runtime_error if the token is invalid or out of range
	 */
	static long long parseInteger(const std::string& token)
	{
		std::size_t parsed = 0;
		try {
			const long long value = std::stoll(token, &parsed, 10);
			if (parsed != token.size())
				throw std::runtime_error("Invalid integer token: " + quote(token));
			return value;
		} catch (const std::invalid_argument&) {
			throw std::runtime_error("Invalid integer token: " + quote(token));
		} catch (const std::out_of_range&) {
			throw std::runtime_error("Integer token is out of range: " + quote(token));
		}
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Parse a real-number token
	 *
	 * \param[in] token Token to parse
	 *
	 * \returns The represented finite real value
	 * 	hrows std::runtime_error if the token is invalid, non-finite, or out of range
	 */
	static RealType parseDouble(const std::string& token)
	{
		std::size_t parsed = 0;
		try {
			const RealType value = std::stod(token, &parsed);
			if (parsed != token.size() || !std::isfinite(value))
				throw std::runtime_error("Invalid double token: " + quote(token));
			return value;
		} catch (const std::invalid_argument&) {
			throw std::runtime_error("Invalid double token: " + quote(token));
		} catch (const std::out_of_range&) {
			throw std::runtime_error("Double token is out of range: " + quote(token));
		}
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Compare two tokens and report a mismatch
	 *
	 * \param[in] token1   Token from the first file
	 * \param[in] token2   Token from the second file
	 * \param[in] position One-based token position
	 *
	 * 	hrows std::runtime_error if the token classes or values differ
	 */
	void compareTokens(const std::string& token1,
	                   const std::string& token2,
	                   std::size_t        position) const
	{
		const TokenClass class1 = classify(token1);
		const TokenClass class2 = classify(token2);
		if (class1 != class2) {
			fail(position,
			     quote(token1),
			     quote(token2),
			     std::string("different token classes (") + className(class1) + " and "
			         + className(class2) + ")");
		}

		switch (class1) {
		case TokenClass::Word:
			if (token1 != token2)
				fail(position, quote(token1), quote(token2), "words differ");
			return;
		case TokenClass::Integer:
		{
			const long long value1 = parseInteger(token1);
			const long long value2 = parseInteger(token2);
			if (value1 != value2)
				fail(position,
				     quote(token1),
				     quote(token2),
				     "integer values differ");
			return;
		}
		case TokenClass::Double:
		{
			const RealType value1 = parseDouble(token1);
			const RealType value2 = parseDouble(token2);
			if (!checkRealNumbers(value1, value2))
				fail(
				    position, quote(token1), quote(token2), "double values differ");
			return;
		}
		}
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Compare two real values using the configured absolute tolerance
	 *
	 * \param[in] value1 Value from the first file
	 * \param[in] value2 Value from the second file
	 *
	 * \returns True if the absolute difference does not exceed the tolerance
	 */
	bool checkRealNumbers(RealType value1, RealType value2) const
	{
		return std::abs(value1 - value2) <= tolerance_;
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Enclose a token in double quotation marks
	 *
	 * \param[in] token Token to quote
	 *
	 * \returns The quoted token
	 */
	static std::string quote(const std::string& token) { return "\"" + token + "\""; }

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Throw a comparison-failure exception
	 *
	 * \param[in] position One-based token position
	 * \param[in] token1   Display text for the first file's token
	 * \param[in] token2   Display text for the second file's token
	 * \param[in] reason   Explanation of the mismatch
	 *
	 * 	hrows std::runtime_error unconditionally
	 */
	[[noreturn]] static void fail(std::size_t        position,
	                              const std::string& token1,
	                              const std::string& token2,
	                              const std::string& reason)
	{
		throw std::runtime_error("Token " + std::to_string(position) + ": " + reason
		                         + "; first file has " + token1 + ", second file has "
		                         + token2);
	}

	RealType    tolerance_;
	std::string ignoredLinePrefix_;
};

} // namespace PsimagLite

#endif // TEXT_AND_NUMBERS_CHECKER_HPP
