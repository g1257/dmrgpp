#ifndef TEXT_AND_NUMBERS_CHECKER_HPP
#define TEXT_AND_NUMBERS_CHECKER_HPP

#include <cmath>
#include <cstddef>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <string>

namespace PsimagLite {

class TextAndNumbersChecker {
public:
	using RealType = double;

	explicit TextAndNumbersChecker(RealType tolerance)
	    : tolerance_(tolerance)
	{
		if (!std::isfinite(tolerance_))
			throw std::invalid_argument("The tolerance must be finite");
	}

	void run(const std::string& file1, const std::string& file2) const
	{
		std::ifstream stream1(file1);
		if (!stream1)
			throw std::runtime_error("Cannot open first file: " + file1);

		std::ifstream stream2(file2);
		if (!stream2)
			throw std::runtime_error("Cannot open second file: " + file2);

		std::size_t position = 1;
		for (;;) {
			std::string token1;
			std::string token2;
			const bool hasToken1 = readToken(stream1, token1);
			const bool hasToken2 = readToken(stream2, token2);

			if (!hasToken1 && !hasToken2)
				return;

			if (hasToken1 != hasToken2) {
				const std::string shown1 = hasToken1 ? quote(token1) : "<end of file>";
				const std::string shown2 = hasToken2 ? quote(token2) : "<end of file>";
				fail(position, shown1, shown2, "the files contain different numbers of tokens");
			}

			compareTokens(token1, token2, position);
			++position;
		}
	}

private:
	enum class TokenClass { Word, Integer, Double };

	static bool isAsciiWhitespace(char c)
	{
		switch (c) {
		case ' ':
		case '\t':
		case '\n':
		case '\r':
		case '\f':
		case '\v':
			return true;
		default:
			return false;
		}
	}

	static bool readToken(std::istream& stream, std::string& token)
	{
		token.clear();
		char c = 0;
		while (stream.get(c)) {
			if (!isAsciiWhitespace(c)) {
				token.push_back(c);
				break;
			}
		}

		if (token.empty())
			return false;

		while (stream.get(c)) {
			if (isAsciiWhitespace(c))
				break;
			token.push_back(c);
		}

		return true;
	}

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

	void compareTokens(const std::string& token1,
	                   const std::string& token2,
	                   std::size_t position) const
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
		case TokenClass::Integer: {
			const long long value1 = parseInteger(token1);
			const long long value2 = parseInteger(token2);
			if (value1 != value2)
				fail(position, quote(token1), quote(token2), "integer values differ");
			return;
		}
		case TokenClass::Double: {
			const RealType value1 = parseDouble(token1);
			const RealType value2 = parseDouble(token2);
			if (!checkRealNumbers(value1, value2))
				fail(position, quote(token1), quote(token2), "double values differ");
			return;
		}
		}
	}

	bool checkRealNumbers(RealType, RealType) const
	{
		throw std::logic_error("checkRealNumbers is not implemented");
	}

	static std::string quote(const std::string& token) { return "\"" + token + "\""; }

	[[noreturn]] static void fail(std::size_t position,
	                              const std::string& token1,
	                              const std::string& token2,
	                              const std::string& reason)
	{
		throw std::runtime_error("Token " + std::to_string(position) + ": " + reason
		                         + "; first file has " + token1 + ", second file has " + token2);
	}

	RealType tolerance_;
};

} // namespace PsimagLite

#endif // TEXT_AND_NUMBERS_CHECKER_HPP
