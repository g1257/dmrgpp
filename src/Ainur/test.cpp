#include <iostream>
#include <fstream>
#include "Vector.h"
#include "TypeToString.h"
#include "PsimagLite.h"
#include "Tokenizer.h"

typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

class Ainur {

public:

	Ainur(PsimagLite::String filename)
	{
		std::ifstream fin(filename.c_str());
		PsimagLite::String str;

		fin.seekg(0, std::ios::end);
		str.reserve(fin.tellg());
		fin.seekg(0, std::ios::beg);

		str.assign((std::istreambuf_iterator<char>(fin)),
		            std::istreambuf_iterator<char>());
		fin.close();

		replaceAndStoreEscaped(escapedChars_, str);

		replaceAndStoreQuotes(vecStr_, str, '"');

		replaceAndStoreQuotes(vecChar_, str, '\'');

		removeComments(str);

		replaceAndStoreBraces(vecBrace_, str);

		std::cout<<str<<"\n";
		printContainers(std::cout);

		VectorStringType statements;
		PsimagLite::tokenizer(str, statements, ";");

		procStatements(statements);
	}

private:

	void printContainers(std::ostream& os) const
	{
		os<<"------------\n";
		os<<"------------\n";
		os<<vecStr_;
		os<<"------------\n";
		os<<vecChar_<<"\n";
		os<<"------------\n";
		os<<vecBrace_<<"\n";
	}

	template<typename T>
	void replaceAndStoreQuotes(T& t, PsimagLite::String& str, char q)
	{
		SizeType l = str.length();
		bool openQuote = false;
		PsimagLite::String newStr("");
		PsimagLite::String buffer("");
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == q) {
				if (openQuote) {
					PsimagLite::String metaString("@ ");
					metaString[1] = getMetaChar(q);
					metaString += ttos(t.size());
					newStr += metaString;
					pushInto(t, buffer);
					buffer = "";
					openQuote = false;
				} else {
					openQuote = true;
				}

				continue;
			}

			if (openQuote) {
				buffer += str[i];
			} else {
				newStr += str[i];
			}
		}

		str = newStr;
	}

	void replaceAndStoreEscaped(PsimagLite::String& t, PsimagLite::String& str)
	{
		SizeType l = str.length();
		PsimagLite::String newStr("");

		char q = '\\';
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == q) {
				newStr += "@e" + ttos(t.length());
				if (l == i + 1)
					err("Syntax Error (escaped): " + getContext(str) + "\n");
				t += str[++i];
				continue;
			}

			newStr += str[i];
		}

		str = newStr;
	}

	void replaceAndStoreBraces(VectorStringType& t, PsimagLite::String& str)
	{
		SizeType l = str.length();
		SizeType openBrace = 0;
		PsimagLite::String newStr("");
		PsimagLite::String buffer("");
		char qOpen = '{';
		char qClose = '}';
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == qOpen) {
				++openBrace;
				buffer += qOpen;
				continue;
			}

			if (str[i] == qClose) {
				if (openBrace == 0)
					err("Syntax Error (closing brace): " + getContext(str) + "\n");
				--openBrace;
				buffer += qClose;
				if (openBrace > 0)
					continue;
				newStr += "@b" + ttos(t.size()) + ";";
				t.push_back(buffer);
				buffer = "";
				continue;
			}

			if (openBrace > 0) {
				buffer += str[i];
			} else {
				newStr += str[i];
			}
		}

		str = newStr;
	}


	PsimagLite::String getContext(const PsimagLite::String& str,
	                              SizeType n = 10) const
	{
		SizeType l = str.length();
		SizeType start = (l < n) ? 0 : l - n;
		return str.substr(start);
	}

	void pushInto(PsimagLite::String& dest, PsimagLite::String src) const
	{
		dest += src;
	}

	void pushInto(VectorStringType& dest,
	              PsimagLite::String src) const
	{
		dest.push_back(src);
	}

	char getMetaChar(char q) const
	{
		if (q == '"') return 's';
		if (q == '\'') return 'q';
		err("getMetaChar\n");
		return 0;
	}

	void removeTrailingWhitespace(PsimagLite::String& s) const
	{
		SizeType start = 0;
		SizeType l = s.length();
		for (SizeType i = 0; i < l; ++i) {
			if (isWhitespace(s[i]) || isEOL(s[i]))
				start = i + 1;
			else
				break;
		}

		if (start == l) {
			s = "";
			return;
		}


		PsimagLite::String newStr = s.substr(start);
		l = newStr.length();
		SizeType end = l;
		for (int i = l - 1; i >= 0; --i) {
			if (isWhitespace(newStr[i]) || isEOL(newStr[i]))
				end = i;
			else
				break;
		}

		s = newStr.substr(0, end);
	}

	void removeComments(PsimagLite::String& str) const
	{
		SizeType l = str.length();
		PsimagLite::String newStr("");
		char qComment = '#';
		bool comment = false;
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == qComment) {
				comment = true;
				continue;
			}

			if (isEOL(str[i]))
				comment = false;

			if (comment)
				continue;

			newStr += str[i];
		}

		str = newStr;
	}

	bool isWhitespace(char c) const
	{
		return (c == ' ' || c == '\t');
	}

	bool isEOL(char c) const
	{
		return (c == '\n' || c == '\r');
	}

	void procStatements(VectorStringType& s)
	{
		SizeType n = s.size();
		for (SizeType i = 0; i < n; ++i) {
			removeTrailingWhitespace(s[i]);
			procStatement(s[i]);
		}
	}

	void procStatement(const PsimagLite::String& s)
	{
		bool bEq = (s.find("=") != PsimagLite::String::npos);
		bool bScOrFunc = (s.find("@b") != PsimagLite::String::npos);

		if (bEq && bScOrFunc) {
			std::cout<<s<<" <---- syntax error\n";
			return;
		}

		if  (bEq) {
			VectorStringType leftAndRight;
			PsimagLite::tokenizer(s, leftAndRight, "=");
			if (leftAndRight.size() != 2)
				err("Syntax error: " + s + "\n");
			VectorStringType lhs;
			PsimagLite::tokenizer(leftAndRight[0], lhs, " ");
			SizeType storageIndex = procLeftEquality(lhs, s);
			std::cerr<<"[" << storageIndex <<"] <--- " << leftAndRight[1] <<"\n";
 		}

		if (bScOrFunc)
			err("Unimplemented: function or scope\n");
	}

	SizeType procLeftEquality(const VectorStringType& lhs,
	                          PsimagLite::String context)
	{
		SizeType l = lhs.size();
		if (l == 0 || l > 2)
			err("Nothing or too much on left? " + context + "\n");
		int x = -1;
		PsimagLite::String name = lhs[0];
		if (l == 1)
			x = storageIndexByName(name);
		if (l == 2) {
			if (lhs[0] != "let" && lhs[0] != "require" && lhs[0] != "const")
				err("Expected let require or const " + context + "\n");
			name = lhs[1];
			x = assignStorageByName(name);
		}

		if (x < 0)
			err("Undeclared variable " + name + "\n");

		return x;
	}

	int assignStorageByName(PsimagLite::String name)
	{
		int x = storageIndexByName(name);
		if (x >= 0)
			err("Already in scope " + name + "\n");
		names_.push_back(name);
		return names_.size() - 1;
	}

	int storageIndexByName(PsimagLite::String name) const
	{
		VectorStringType::const_iterator it = std::find(names_.begin(), names_.end(), name);
		if (it == names_.end())
			return -1;
		return it - names_.begin();
	}

	PsimagLite::String escapedChars_;
	VectorStringType vecStr_;
	PsimagLite::String vecChar_;
	VectorStringType vecBrace_;
	VectorStringType names_;
};

int main(int argc, char** argv)
{
	if (argc == 1) return 1;
	Ainur ainur(argv[1]);
}
