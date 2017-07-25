#ifndef _AINUR_H_
#define _AINUR_H_
#include <iostream>
#include <fstream>
#include "Vector.h"
#include "TypeToString.h"
#include "PsimagLite.h"
#include "AinurStatements.h"

namespace PsimagLite {

class Ainur {

public:

	typedef Vector<String>::Type VectorStringType;
	typedef AinurStatements AinurStatementsType;
	typedef AinurStatementsType::AinurLexicalType AinurLexicalType;

	Ainur(String str)
	{
		replaceAtAndCheck(str);

		replaceAndStoreEscaped(escapedChars_, str);

		replaceAndStoreQuotes(vecStr_, str, '"');

		replaceAndStoreQuotes(vecChar_, str, '\'');

		removeComments(str);

		replaceAndStoreBraces(vecBrace_, str);

		//std::cout<<str<<"\n";
		//printContainers(std::cout);

		VectorStringType statements;
		split(statements, str, ";");

		procStatements(statements);
	}

	template<typename SomeType>
	void readValue(SomeType& t, String label) const
	{
		return statements_.readValue(t, label);
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

	void removeComments(String& str)
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

			if (AinurLexicalType::isEOL(str[i]))
				comment = false;

			if (comment)
				continue;

			newStr += str[i];
		}

		str = newStr;
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

	bool allowedChar(unsigned char c) const
	{
		if (AinurLexicalType::isWhitespace(c) || AinurLexicalType::isEOL(c))
			return true;
		if (c < 33 ||c > 126) return false;
		if (c == 95 || c == 96) return false;
		return true;
	}

	void procStatements(VectorStringType& s)
	{
		SizeType n = s.size();
		for (SizeType i = 0; i < n; ++i)
			statements_.push(s[i]);
	}

	String escapedChars_;
	VectorStringType vecStr_;
	String vecChar_;
	VectorStringType vecBrace_;
	AinurStatements statements_;
};

}
#endif // _AINUR_H_
