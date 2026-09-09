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
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace PsimagLite {

//===========================================================================//
/*!
 * \class TextAndNumbersChecker
 *
 * \brief Compare whitespace-separated words and numbers in two text files
 *
 * Words are compared verbatim, integers by value, and real or complex numbers
 * using an absolute tolerance. Complex-number error is the magnitude of the
 * difference. Tokens in corresponding positions must have the same class.
 * Lines beginning with an optional prefix can be excluded from the comparison.
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
	 *                               corresponding real- or complex-number tokens
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
	void run(const std::filesystem::path& file1, const std::filesystem::path& file2) const
	{
		std::ifstream stream1(file1);
		if (!stream1)
			throw std::runtime_error("Cannot open first file: " + file1.string());

		std::ifstream stream2(file2);
		if (!stream2)
			throw std::runtime_error("Cannot open second file: " + file2.string());

		TokenReader reader1(stream1, ignoredLinePrefix_);
		TokenReader reader2(stream2, ignoredLinePrefix_);
		std::size_t position = 1;
		for (;;) {
			std::string token1;
			std::string token2;
			const bool  hasToken1 = reader1.read(token1);
			if (!hasToken1 && !stream1.eof())
				throw std::runtime_error("Error while reading first file: "
				                         + file1.string());
			const bool hasToken2 = reader2.read(token2);
			if (!hasToken2 && !stream2.eof())
				throw std::runtime_error("Error while reading second file: "
				                         + file2.string());

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

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Compare final indexed records for one observable label
	 *
	 * Each selected row has fields `(site, value, time, label, superdensity)`.
	 * Records are indexed by the exact parsed `(site, time, label)` tuple. If an
	 * index occurs repeatedly in either file, its last occurrence is authoritative.
	 * The final key sets must match exactly; both complex data fields are compared
	 * by error magnitude using the configured tolerance.
	 *
	 * \param[in] referenceFile  Path to the reference file
	 * \param[in] actualFile     Path to the generated output file
	 * \param[in] observableLabel Exact label selecting observable rows
	 *
	 * 	hrows std::runtime_error if a file cannot be opened, a selected row is
	 *         malformed, no selected records exist, or the final maps differ
	 */
	void runObservables(const std::filesystem::path& referenceFile,
	                    const std::filesystem::path& actualFile,
	                    const std::string&           observableLabel) const
	{
		if (observableLabel.empty())
			throw std::invalid_argument("The observable label must not be empty");

		const ObservableMap reference
		    = readObservables(referenceFile, observableLabel, "reference");
		const ObservableMap actual = readObservables(actualFile, observableLabel, "actual");

		auto referenceIter = reference.begin();
		auto actualIter    = actual.begin();
		while (referenceIter != reference.end() || actualIter != actual.end()) {
			if (actualIter == actual.end()
			    || (referenceIter != reference.end()
			        && referenceIter->first < actualIter->first))
				throw std::runtime_error(
				    "Observable index missing from actual file: "
				    + formatIndex(referenceIter->first));

			if (referenceIter == reference.end()
			    || actualIter->first < referenceIter->first)
				throw std::runtime_error("Observable index extra in actual file: "
				                         + formatIndex(actualIter->first));

			compareObservableData(
			    referenceIter->first, referenceIter->second, actualIter->second);
			++referenceIter;
			++actualIter;
		}
	}

private:

	using ObservableIndex = std::tuple<std::size_t, RealType, std::string>;

	struct ObservableData {
		std::complex<RealType> value;
		std::complex<RealType> superdensity;
	};

	using ObservableMap = std::map<ObservableIndex, ObservableData>;

	enum class TokenClass
	{
		Word,
		Integer,
		Double,
		Complex
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
			if (!std::isspace(static_cast<unsigned char>(c))) {
				token.push_back(c);
				break;
			}
		}

		if (token.empty())
			return false;

		while (stream.get(c)) {
			if (std::isspace(static_cast<unsigned char>(c)))
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
	 * \brief Test whether text equals an ASCII word without regard to case
	 */
	static bool equalsAsciiNoCase(const std::string& text, std::size_t offset, const char* word)
	{
		std::size_t index = offset;
		for (; *word != '\0'; ++word, ++index) {
			if (index >= text.size()
			    || std::tolower(static_cast<unsigned char>(text[index])) != *word)
				return false;
		}
		return index == text.size();
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Recognize one complex component with deterministic linear scanning
	 */
	static bool isComplexComponentSyntax(const std::string& component)
	{
		if (component.empty())
			return false;

		std::size_t position = 0;
		if (component[position] == '+' || component[position] == '-')
			++position;
		if (position == component.size())
			return false;

		if (equalsAsciiNoCase(component, position, "inf")
		    || equalsAsciiNoCase(component, position, "infinity")
		    || equalsAsciiNoCase(component, position, "nan"))
			return true;

		bool hasDigits = false;
		while (position < component.size()
		       && std::isdigit(static_cast<unsigned char>(component[position]))) {
			hasDigits = true;
			++position;
		}

		if (position < component.size() && component[position] == '.') {
			++position;
			while (position < component.size()
			       && std::isdigit(static_cast<unsigned char>(component[position]))) {
				hasDigits = true;
				++position;
			}
		}

		if (!hasDigits)
			return false;

		if (position < component.size()
		    && (component[position] == 'e' || component[position] == 'E')) {
			++position;
			if (position < component.size()
			    && (component[position] == '+' || component[position] == '-'))
				++position;

			const std::size_t exponentStart = position;
			while (position < component.size()
			       && std::isdigit(static_cast<unsigned char>(component[position])))
				++position;
			if (position == exponentStart)
				return false;
		}

		return position == component.size();
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Recognize strict `(real,imaginary)` complex-token syntax
	 */
	static bool isComplexSyntax(const std::string& token)
	{
		if (token.size() < 5 || token.front() != '(' || token.back() != ')')
			return false;

		const std::size_t comma = token.find(',', 1);
		if (comma == std::string::npos || token.find(',', comma + 1) != std::string::npos)
			return false;

		return isComplexComponentSyntax(token.substr(1, comma - 1))
		    && isComplexComponentSyntax(token.substr(comma + 1, token.size() - comma - 2));
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Recognize strict signed decimal-integer syntax
	 */
	static bool isIntegerSyntax(const std::string& token)
	{
		if (token.empty())
			return false;

		std::size_t position = 0;
		if (token[position] == '+' || token[position] == '-')
			++position;
		if (position == token.size())
			return false;

		for (; position < token.size(); ++position) {
			if (!std::isdigit(static_cast<unsigned char>(token[position])))
				return false;
		}
		return true;
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Recognize strict decimal or exponent real-number syntax
	 */
	static bool isDoubleSyntax(const std::string& token)
	{
		if (token.empty())
			return false;

		std::size_t position = 0;
		if (token[position] == '+' || token[position] == '-')
			++position;

		const std::size_t integerStart = position;
		while (position < token.size()
		       && std::isdigit(static_cast<unsigned char>(token[position])))
			++position;
		const bool hasIntegerDigits = position != integerStart;

		bool hasDecimalPoint = false;
		if (position < token.size() && token[position] == '.') {
			hasDecimalPoint = true;
			++position;
			const std::size_t fractionStart = position;
			while (position < token.size()
			       && std::isdigit(static_cast<unsigned char>(token[position])))
				++position;
			if (!hasIntegerDigits && position == fractionStart)
				return false;
		}

		bool hasExponent = false;
		if (position < token.size() && (token[position] == 'e' || token[position] == 'E')) {
			hasExponent = true;
			++position;
			if (position < token.size()
			    && (token[position] == '+' || token[position] == '-'))
				++position;
			const std::size_t exponentStart = position;
			while (position < token.size()
			       && std::isdigit(static_cast<unsigned char>(token[position])))
				++position;
			if (position == exponentStart)
				return false;
		}

		return position == token.size() && (hasDecimalPoint || hasExponent)
		    && (hasDecimalPoint || hasIntegerDigits);
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Classify a token as a word, integer, real, or complex number
	 *
	 * \param[in] token Token to classify
	 *
	 * \returns The token class
	 */
	static TokenClass classify(const std::string& token)
	{
		if (isIntegerSyntax(token))
			return TokenClass::Integer;
		if (isDoubleSyntax(token))
			return TokenClass::Double;
		if (isComplexSyntax(token))
			return TokenClass::Complex;
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
		case TokenClass::Complex:
			return "complex";
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
	 * \brief Parse a component of a complex-number token
	 *
	 * \param[in] component Component text to parse
	 * \param[in] token     Complete complex token for diagnostics
	 *
	 * \returns The represented finite real value
	 * 	hrows std::runtime_error if the component is invalid, non-finite, or
	 *         out of range
	 */
	static RealType parseComplexComponent(const std::string& component,
	                                      const std::string& token)
	{
		std::size_t parsed = 0;
		try {
			const RealType value = std::stod(component, &parsed);
			if (parsed != component.size())
				throw std::runtime_error("Invalid complex token: " + quote(token));
			if (!std::isfinite(value))
				throw std::runtime_error(
				    "Complex token has a non-finite component: " + quote(token));
			return value;
		} catch (const std::invalid_argument&) {
			throw std::runtime_error("Invalid complex token: " + quote(token));
		} catch (const std::out_of_range&) {
			throw std::runtime_error("Complex component is out of range: "
			                         + quote(token));
		}
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Parse a complex-number token in `(real,imaginary)` form
	 *
	 * \param[in] token Token to parse
	 *
	 * \returns The represented finite complex value
	 * 	hrows std::runtime_error if either component is invalid, non-finite, or
	 *         out of range
	 */
	static std::complex<RealType> parseComplex(const std::string& token)
	{
		const std::size_t comma = token.find(',');
		if (token.size() < 5 || token.front() != '(' || token.back() != ')'
		    || comma == std::string::npos)
			throw std::runtime_error("Invalid complex token: " + quote(token));

		const std::string realPart      = token.substr(1, comma - 1);
		const std::string imaginaryPart = token.substr(comma + 1, token.size() - comma - 2);
		return { parseComplexComponent(realPart, token),
			 parseComplexComponent(imaginaryPart, token) };
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Parse a non-negative site index with full range checking
	 */
	static std::size_t parseSite(const std::string& token)
	{
		if (token.empty() || token.find_first_not_of("0123456789") != std::string::npos)
			throw std::runtime_error("Invalid site token: " + quote(token));

		std::size_t parsed = 0;
		try {
			const unsigned long long value = std::stoull(token, &parsed, 10);
			if (parsed != token.size()
			    || value > std::numeric_limits<std::size_t>::max())
				throw std::runtime_error("Site token is out of range: "
				                         + quote(token));
			return static_cast<std::size_t>(value);
		} catch (const std::invalid_argument&) {
			throw std::runtime_error("Invalid site token: " + quote(token));
		} catch (const std::out_of_range&) {
			throw std::runtime_error("Site token is out of range: " + quote(token));
		}
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Parse an integer-looking or real-number time token
	 */
	static RealType parseObservableTime(const std::string& token)
	{
		const TokenClass tokenClass = classify(token);
		if (tokenClass != TokenClass::Integer && tokenClass != TokenClass::Double)
			throw std::runtime_error("Invalid time token: " + quote(token));
		return parseDouble(token);
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Test whether a token has unsigned decimal-integer syntax
	 */
	static bool isUnsignedDecimalSyntax(const std::string& token)
	{
		if (token.empty())
			return false;
		for (const char character : token) {
			if (!std::isdigit(static_cast<unsigned char>(character)))
				return false;
		}
		return true;
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Identify selected rows even when a missing field shifts the label
	 */
	static bool isObservableCandidate(const std::vector<std::string>& tokens,
	                                  const std::string&              observableLabel)
	{
		bool containsLabel         = false;
		bool labelAtRecordPosition = false;
		bool containsComplex       = false;
		for (std::size_t index = 0; index < tokens.size(); ++index) {
			const std::string& token = tokens[index];
			if (token == observableLabel) {
				containsLabel = true;
				labelAtRecordPosition
				    = labelAtRecordPosition || (index >= 2 && index <= 4);
			}
			containsComplex = containsComplex || classify(token) == TokenClass::Complex;
		}

		return containsLabel
		    && (labelAtRecordPosition
		        || (!tokens.empty() && isUnsignedDecimalSyntax(tokens.front()))
		        || containsComplex);
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Read and collapse selected observable rows from one file
	 */
	static ObservableMap readObservables(const std::filesystem::path& file,
	                                     const std::string&           observableLabel,
	                                     const std::string&           fileRole)
	{
		std::ifstream stream(file);
		if (!stream)
			throw std::runtime_error("Cannot open " + fileRole
			                         + " file: " + file.string());

		ObservableMap data;
		std::string   line;
		std::size_t   lineNumber = 0;
		while (std::getline(stream, line)) {
			++lineNumber;
			std::istringstream       lineStream(line);
			std::vector<std::string> tokens;
			std::string              token;
			while (lineStream >> token)
				tokens.push_back(token);

			if (!isObservableCandidate(tokens, observableLabel))
				continue;

			const std::string location = "Malformed observable record in " + fileRole
			    + " file at line " + std::to_string(lineNumber);
			if (tokens.size() != 5 || tokens[3] != observableLabel)
				throw std::runtime_error(location + ": " + line);

			try {
				if (classify(tokens[1]) != TokenClass::Complex
				    || classify(tokens[4]) != TokenClass::Complex)
					throw std::runtime_error(
					    "value or superdensity is not complex");

				const std::size_t     site = parseSite(tokens[0]);
				const RealType        time = parseObservableTime(tokens[2]);
				const ObservableIndex index { site, time, tokens[3] };
				const ObservableData  observableData { parseComplex(tokens[1]),
                                                                      parseComplex(tokens[4]) };
				data[index] = observableData;
			} catch (const std::exception& error) {
				throw std::runtime_error(location + ": " + error.what());
			}
		}

		if (stream.bad())
			throw std::runtime_error("Error while reading " + fileRole
			                         + " file: " + file.string());

		if (data.empty())
			throw std::runtime_error("No observable records for label "
			                         + observableLabel + " in " + fileRole
			                         + " file: " + file.string());
		return data;
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Format a complete observable index for diagnostics
	 */
	static std::string formatIndex(const ObservableIndex& index)
	{
		std::ostringstream output;
		output << "(site=" << std::get<0>(index) << ", time=" << std::get<1>(index)
		       << ", label=" << std::get<2>(index) << ")";
		return output.str();
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Compare both data fields for one common observable index
	 */
	void compareObservableData(const ObservableIndex& index,
	                           const ObservableData&  reference,
	                           const ObservableData&  actual) const
	{
		const RealType valueError = std::abs(reference.value - actual.value);
		if (valueError > tolerance_)
			failObservable(
			    "observable value", index, reference.value, actual.value, valueError);

		const RealType superdensityError
		    = std::abs(reference.superdensity - actual.superdensity);
		if (superdensityError > tolerance_)
			failObservable("superdensity",
			               index,
			               reference.superdensity,
			               actual.superdensity,
			               superdensityError);
	}

	//---------------------------------------------------------------------------//
	/*!
	 * \brief Throw a detailed observable-data mismatch
	 */
	[[noreturn]] void failObservable(const std::string&            field,
	                                 const ObservableIndex&        index,
	                                 const std::complex<RealType>& reference,
	                                 const std::complex<RealType>& actual,
	                                 RealType                      error) const
	{
		std::ostringstream message;
		message << field << " mismatch at " << formatIndex(index)
		        << "; reference=" << reference << ", actual=" << actual
		        << ", error=" << error << ", tolerance=" << tolerance_;
		throw std::runtime_error(message.str());
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
		case TokenClass::Complex:
		{
			const std::complex<RealType> value1 = parseComplex(token1);
			const std::complex<RealType> value2 = parseComplex(token2);
			if (std::abs(value1 - value2) > tolerance_)
				fail(position,
				     quote(token1),
				     quote(token2),
				     "complex values differ");
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
