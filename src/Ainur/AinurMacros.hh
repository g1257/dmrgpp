#ifndef AINURMACROS_HH
#define AINURMACROS_HH
#include "../PsimagLite.h"
#include <map>
#include <string>

namespace PsimagLite
{

class AinurMacros
{

public:

	struct NativeMacro {
		std::string type;
		std::string name;
		std::string value;
	};

	class AinurFunction
	{

		static constexpr SizeType MAX_LINE = 2048;

	public:

		AinurFunction(const std::string& args)
		    : separator_("space")
		    , columnX_(0)
		    , columnY_(1)
		{
			std::string csl = deleteEnclosing(args, '(', ')');
			std::vector<std::string> tokens;
			split(tokens, csl, ",");
			if (tokens.size() == 0) {
				throw RuntimeError(
				    "Parse error for " + args + ". Expected (filename, separator, "
								"column_x, column_y),"
				    + " where separator and columns are "
				      "optional.\n");
			}

			filename_ = deleteEnclosing(tokens[0], '\"');
			if (tokens.size() > 1) {
				separator_ = tokens[1];
			}

			separatorForSplit_ = setSeparatorForSplit(separator_);

			if (tokens.size() > 2) {
				columnX_ = PsimagLite::atoi(tokens[2]);
			}

			if (tokens.size() > 3) {
				columnY_ = PsimagLite::atoi(tokens[3]);
			}

			if (tokens.size() > 4) {
				throw RuntimeError("Parse error for " + args + ": too many arguments\n");
			}

			readData();
		}

		const double& operator()(const double& x) const
		{
			auto it = std::find_if(
			    data_.begin(), data_.end(), [&x](const std::pair<double, double>& pair) {
				    return (pair.first == x);
			    });
			return it->second;
		}

		static std::string deleteEnclosing(const std::string& content,
		    char b)
		{
			return deleteEnclosing(content, b, b);
		}

		static std::string deleteEnclosing(const std::string& content,
		    char b,
		    char e)
		{
			SizeType length = content.size();
			if (length == 0)
				return content;
			SizeType last = length - 1;
			SizeType total = length;
			SizeType start = 0;

			if (content[0] == b) {
				start = 1;
				--total;
			}

			if (last > 0 && content[last] == e) {
				--total;
			}

			return content.substr(start, total);
		}

	private:

		static std::string
		setSeparatorForSplit(const std::string& separator)
		{
			if (separator == "space") {
				return " ";
			} else if (separator == "comma") {
				return ",";
			}

			throw RuntimeError("Wrong separator" + separator + ": only space and comma allowed.\n");
		}

		static bool emptyLine(const std::string& s)
		{
			for (SizeType i = 0; i < s.size(); ++i) {
				if (s[i] != ' ' || s[i] != '\t')
					return false;
			}

			return true;
		}
		void readData()
		{
			std::ifstream fin(filename_);
			if (!fin || !fin.good() || fin.bad()) {
				throw RuntimeError("Could not open " + filename_ + "\n");
			}

			char s[MAX_LINE];
			while (!fin.eof()) {
				fin.getline(s, MAX_LINE);
				if (emptyLine(std::string(s)))
					continue;
				std::vector<std::string> tokens;
				split(tokens, s);
				if (columnX_ > tokens.size() || columnY_ > tokens.size()) {
					throw RuntimeError("Line too small: " + std::string(s) + "\n");
				}

				double vx = PsimagLite::atof(tokens[columnX_]);
				double vy = PsimagLite::atof(tokens[columnY_]);
				std::pair<double, double> pair(vx, vy);
				data_.emplace_back(pair);
			}
		}

		std::string filename_;
		std::string separator_;
		std::string separatorForSplit_;
		SizeType columnX_;
		SizeType columnY_;
		std::vector<std::pair<double, double>> data_;
	};

	AinurMacros()
	    : AINUR_FROM_FILE("AinurFromFile")
	{
		nativeMacros_.push_back(
		    { "function", AINUR_FROM_FILE, "!" + AINUR_FROM_FILE });
	}

	SizeType total() const { return nativeMacros_.size(); }

	const NativeMacro& nativeMacro(SizeType ind) const
	{
		assert(ind < nativeMacros_.size());
		return nativeMacros_[ind];
	}

	std::string procNativeMacro(const std::string& line)
	{
		if (line.substr(1, AINUR_FROM_FILE.size()) == AINUR_FROM_FILE) {
			SizeType start = AINUR_FROM_FILE.size() + 1;
			SizeType len = line.size();
			assert(len > start);
			SizeType rest = len - start;
			std::string content = line.substr(start, rest);
			return addAinurFromFile(content);
		}

		throw RuntimeError("Unknown native macro for " + line + "\n");
	}

	// Expect ("function.txt") or
	// ("function.txt")(3.0)
	std::string valueFromFunction(const std::string& line) const
	{
		std::pair<std::string, std::string> nameValue = getFunctionNameValue(line);

		if (nameValue.second == "") {
			return nameValue.first;
		}

		if (!isAfloat(nameValue.second)) {
			return line;
		}

		// return value
		std::string functionName = nameValue.first;

		// using and defining in same line not allowed
		if (functionNameToIndex_.count(functionName) == 0) {
			throw RuntimeError("valueFromFunction: " + functionName + " undefined.\n");
		}

		SizeType index = functionNameToIndex_.at(functionName);
		assert(index < ainurFunctions_.size());
		std::string argument = AinurFunction::deleteEnclosing(nameValue.second, '(', ')');
		double x = PsimagLite::atof(argument);
		double y = ainurFunctions_[index](x);
		return std::to_string(y);
	}

private:

	// Expect ("function.txt") or
	// ("function.txt")(3.0)
	static std::pair<std::string, std::string>
	getFunctionNameValue(const std::string& line)
	{
		SizeType length = line.size();
		if (length < 4) {
			return std::pair<std::string, std::string>(line, "");
		}

		SizeType last = length - 1;

		if (line[0] != '(' || line[last] != ')') {
			return std::pair<std::string, std::string>(line, "");
		}

		SizeType i = 1;
		for (; i < line.size(); ++i) {
			if (line[i] == ')' && line[i + 1] == '(')
				break;
		}

		// no value found
		if (i == line.size()) {
			return std::pair<std::string, std::string>(line, "");
		}

		return std::pair<std::string, std::string>(
		    line.substr(0, i + 1), line.substr(i + 2, length - i - 3));
	}

	std::string addAinurFromFile(const std::string& content)
	{
		// Uniqueness
		if (functionNameToIndex_.count(content) == 0) {
			functionNameToIndex_[content] = ainurFunctions_.size();
			AinurFunction ainurFunction(content);
			ainurFunctions_.emplace_back(ainurFunction);
		}

		return content;
	}

	const std::string AINUR_FROM_FILE;
	std::vector<NativeMacro> nativeMacros_;
	std::map<std::string, SizeType> functionNameToIndex_;
	std::vector<AinurFunction> ainurFunctions_;
};

} // namespace PsimagLite
#endif // AINURMACROS_HH
