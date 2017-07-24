#ifndef _AINUR_H_
#define _AINUR_H_
#include <iostream>
#include <fstream>
#include "Vector.h"
#include "TypeToString.h"
#include "PsimagLite.h"

namespace PsimagLite {

class Ainur {

public:

	typedef Vector<String>::Type VectorStringType;

	Ainur(String filename, String import)
	{
		std::ifstream fin(filename.c_str());
		String str;

		fin.seekg(0, std::ios::end);
		str.reserve(fin.tellg());
		fin.seekg(0, std::ios::beg);

		str.assign((std::istreambuf_iterator<char>(fin)),
		            std::istreambuf_iterator<char>());
		fin.close();

		str = import + str;
		replaceAtAndCheck(str);

		replaceAndStoreEscaped(escapedChars_, str);

		replaceAndStoreQuotes(vecStr_, str, '"');

		replaceAndStoreQuotes(vecChar_, str, '\'');

		removeComments(str);

		replaceAndStoreBraces(vecBrace_, str);

		std::cout<<str<<"\n";
		printContainers(std::cout);

		VectorStringType statements;
		split(statements, str, ";");

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

	void replaceAtAndCheck(String& str)
	{
		SizeType l = str.length();
		String newStr("");
		char atChar = '@';
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == atChar) {
				newStr += "@a";
				continue;
			}

			if (!allowedChar(str[i]))
				err("Lexical error: Not allowed char " + getContext(str,i) + "\n");

			newStr += str[i];
		}

		str = newStr;
	}

	template<typename T>
	void replaceAndStoreQuotes(T& t, String& str, char q)
	{
		SizeType l = str.length();
		bool openQuote = false;
		String newStr("");
		String buffer("");
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == q) {
				if (openQuote) {
					String metaString("@ ");
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

	void replaceAndStoreEscaped(String& t, String& str)
	{
		SizeType l = str.length();
		String newStr("");

		char q = '\\';
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == q) {
				newStr += "@e" + ttos(t.length());
				if (l == i + 1)
					err("Syntax Error (escaped): " + getContext(str, i) + "\n");
				t += str[++i];
				continue;
			}

			newStr += str[i];
		}

		str = newStr;
	}

	void replaceAndStoreBraces(VectorStringType& t, String& str)
	{
		SizeType l = str.length();
		SizeType openBrace = 0;
		String newStr("");
		String buffer("");
		char qOpen = '{';
		char qClose = '}';
		for (SizeType i = 0; i < l; ++i) {
			if (str[i] == qOpen) {
				++openBrace;
				if (openBrace > 1)
					err("Syntax error (nested braces not allowed) " +
					    getContext(str, i) + "\n");

				buffer += qOpen;
				continue;
			}

			if (str[i] == qClose) {
				if (openBrace == 0)
					err("Syntax Error (closing brace): " + getContext(str, i) + "\n");
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


	String getContext(const String& str,
	                              SizeType start,
	                              SizeType n = 10) const
	{
		SizeType l = str.length();
		SizeType end = start + n;
		if (end >= l) end = l;
		return str.substr(start, end);
	}

	void pushInto(String& dest, String src) const
	{
		dest += src;
	}

	void pushInto(VectorStringType& dest,
	              String src) const
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

	void removeTrailingWhitespace(String& s) const
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


		String newStr = s.substr(start);
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

	void removeComments(String& str) const
	{
		SizeType l = str.length();
		String newStr("");
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

	bool allowedChar(unsigned char c) const
	{
		if (isWhitespace(c) || isEOL(c)) return true;
		if (c < 33 ||c > 126) return false;
		if (c == 95 || c == 96) return false;
		return true;
	}

	void procStatements(VectorStringType& s)
	{
		SizeType n = s.size();
		for (SizeType i = 0; i < n; ++i) {
			removeTrailingWhitespace(s[i]);
			procStatement(s[i]);
		}
	}

	void procStatement(const String& s)
	{
		bool sEq = (s.find("=") != String::npos);
		bool sBrace = (s.find("@b") != String::npos);

		if (sEq && sBrace) {
			std::cout<<s<<" <---- syntax error\n";
			return;
		}

		if  (sEq) {
			VectorStringType leftAndRight;
			split(leftAndRight, s, "=");
			if (leftAndRight.size() != 2)
				err("Syntax error: " + s + "\n");
			VectorStringType lhs;
			split(lhs, leftAndRight[0], " ");
			SizeType storageIndex = procLeftEquality(lhs, s);
			std::cerr<<"[" << storageIndex <<"] <--- " << leftAndRight[1] <<"\n";
 		}

		if (sBrace)
			err("Unimplemented: scope\n");
	}

	SizeType procLeftEquality(const VectorStringType& lhs,
	                          String context)
	{
		SizeType l = lhs.size();
		if (l == 0 || l > 2)
			err("Nothing or too much on left? " + context + "\n");
		int x = -1;
		String name = lhs[0];
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

	int assignStorageByName(String name)
	{
		int x = storageIndexByName(name);
		if (x >= 0)
			err("Already in scope " + name + "\n");
		names_.push_back(name);
		return names_.size() - 1;
	}

	int storageIndexByName(String name) const
	{
		VectorStringType::const_iterator it = std::find(names_.begin(), names_.end(), name);
		if (it == names_.end())
			return -1;
		return it - names_.begin();
	}

	String escapedChars_;
	VectorStringType vecStr_;
	String vecChar_;
	VectorStringType vecBrace_;
	VectorStringType names_;
};

}
#endif // _AINUR_H_
