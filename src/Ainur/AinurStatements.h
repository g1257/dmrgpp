#ifndef AINURSTATEMENT_H
#define AINURSTATEMENT_H
#include "AinurStore.h"
#include "AinurReadable.h"

namespace PsimagLite {

class AinurStatements {

public:

	typedef AinurReadable::RealType RealType;
	typedef AinurReadable::VectorStoreType VectorStoreType;
	typedef AinurReadable::StoreType StoreType;
	typedef StoreType::AinurLexicalType AinurLexicalType;
	typedef AinurLexicalType::VectorStringType VectorStringType;

	AinurStatements(const VectorStringType& vecStr,
	                const String& vecChar,
	                const String& escapedChars,
	                const VectorStringType& vecBrace)
	    : vecStr_(vecStr),
	      vecChar_(vecChar),
	      escapedChars_(escapedChars),
	      vecBrace_(vecBrace),
	      readable_(names_, storage_)
	{}

	VectorStringType push(const String& s2, String prefix)
	{
		VectorStringType emptyStringVector;
		String s = s2;
		AinurLexicalType::removeTrailingBlanks(s);
		AinurLexicalType::removeTrailingBlanks(prefix);
		if (s == "") return emptyStringVector;
		SizeType storageIndex = 0;

		VectorStringType leftAndRight;
		split(leftAndRight, s, "=");
		AinurLexicalType::removeTrailingWhitespace(leftAndRight);
		if (leftAndRight.size() != 1 && leftAndRight.size() != 2)
			err("Syntax error: " + s + "\n");

		String identifier = procLeftEquality(storageIndex, leftAndRight[0], prefix, s);
		if (leftAndRight.size() == 1) return emptyStringVector;

		if (storageIndex >= storage_.size())
			err("StorageIndex too big\n");

		unescape(leftAndRight[1]);
		unescape(identifier);
		StoreType& store = storage_[storageIndex];
		if (store.type() == StoreType::SCALAR && store.subType() == StoreType::GROUP)  {
			String right = leftAndRight[1];
			SizeType last = right.length();
			--last;
			bool inBraces = (last < right.length() && right[0] == '{' && right[last] == '}');

			if (!inBraces)
				err("Group must be in braces, " + leftAndRight[0] + "\n");

			leftAndRight[1] = (last < 2) ? "" : right.substr(1,last - 1);
			leftAndRight[0] = identifier + ":";
			return leftAndRight;

		}

		store.setRhs(leftAndRight[1], identifier);
		return emptyStringVector;
	}

	void secondPass()
	{
		SizeType n = storage_.size();
		assert(n == names_.size());
		for (SizeType i = 0; i < n; ++i) {
			const String& name = names_[i];
			StoreType& store = storage_[i];
			SizeType total = store.valueSize();

			for (SizeType j = 0; j < total; ++j) {
				String replacement = solveExpression(store.value(j, name), i);
				store.value(j, name) = replacement;
			}
		}
	}

	AinurReadable& readable() { return readable_; }

	const AinurReadable& readable() const { return readable_; }

private:

	String procLeftEquality(SizeType& y,
	                        String s,
	                        String prefix,
	                        String context)
	{
		VectorStringType dotified;
		VectorStringType lhs;
		split(lhs,s," ");
		SizeType l = lhs.size();
		// require vector.vector.integer FiniteLoops
		// matrix.integer FiniteLoops
		// identifier
		if (l < 1 || l > 3)
			err("Too much or too little on left? -- " + context + " --\n");
		int x = -1;

		if (l == 1) { // identifier
			split(dotified,prefix + lhs[0],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = readable_.storageIndexByName(dotified[0]);
		} else if (l == 2) { // matrix.integer FiniteLoops
			split(dotified,prefix + lhs[1],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = assignStorageByName(dotified[0]);
			storage_.push_back(StoreType(lhs[0],""));
		} else if (l == 3) {
			// require vector.vector.integer FiniteLoops
			split(dotified,prefix + lhs[2],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = assignStorageByName(dotified[0]);
			storage_.push_back(StoreType(lhs[1], lhs[0]));
		}

		if (x < 0)
			err("Undeclared variable " + dotified[0] + "\n");

		y = x;

		return dotified[0];
	}

	StoreType::Attribute getAttribute(String s, String context) const
	{
		if (s == "require") return StoreType::REQUIRED;
		if (s == "const") return StoreType::CONST;

		err("Expected let require or const " + context + "\n");
		return StoreType::NONE;
	}

	int assignStorageByName(String name)
	{
		int x = readable_.storageIndexByName(name);
		if (x >= 0)
			err("Already in scope " + name + "\n");
		names_.push_back(name);
		return names_.size() - 1;
	}

	void unescape(String& s) const
	{
		SizeType l = s.length();
		String newStr("");
		for (SizeType i = 0; i < l; ++i) {
			if (s[i] == '@') {
				newStr += getReplacement(i, s, l);
				continue;
			}

			newStr += s[i];
		}

		s = newStr;
	}

	// i @, i+1 s, i+2 0, i+3 @
	String getReplacement(SizeType& i, String s, SizeType l) const
	{
		assert(i < l);
		String oneChar(" ");
		oneChar[0] = s[i];
		if (i + 3 >= l) return oneChar;
		char c = s[++i];
		String number;
		SizeType j = i + 1;
		for (; j < l; ++j) {
			if (s[j] == '@') break;
			if (s[j] < 48 || s[j] > 57)
				err("Error while replacing string\n");
			number += s[j];
		}

		if (s[j] != '@')
			err("Error while replacing string, no final @ found\n");

		i = j + 1;
		SizeType n = atoi(number.c_str());
		return getReplacement(c, n);
	}

	String getReplacement(char c, SizeType n) const
	{
		if (c == 's') {
			if (n >= vecStr_.size())
				err("Error while replacing string, index too big\n");
			return vecStr_[n];
		}

		if (c == 'b') {
			if (n >= vecBrace_.size())
				err("Error while replacing string, index too big\n");
			return vecBrace_[n];
		}

		if (c == 'q') {
			if (n >= vecChar_.length())
				err("Error while replacing string, index too big\n");
			String oneChar(" ");
			oneChar[0] = vecChar_[n];
			return oneChar;
		}

		if (c == 'e') {
			if (n >= escapedChars_.length())
				err("Error while replacing string, index too big\n");
			String oneChar(" ");
			oneChar[0] = escapedChars_[n];
			return oneChar;
		}

		err("Expected s or b or q or e after replacement\n");
		return "";
	}

	// FIXME: Must be generalized to general expressions
	String solveExpression(String s, SizeType end) const
	{
		assert(end <= names_.size());
		assert(names_.size() == storage_.size());
		for (SizeType i = 0; i < end; ++i) {
			if (s != names_[i]) continue;
			const StoreType& store = storage_[i];
			if (store.valueSize() > 1)
				err("Cannot replace vectors yet\n");

			s = store.value(0, names_[i]);
			break;
		}

		return s;
	}

	const VectorStringType& vecStr_;
	const String& vecChar_;
	const String& escapedChars_;
	const VectorStringType& vecBrace_;
	VectorStringType names_;
	VectorStoreType storage_;
	AinurReadable readable_;
}; // class AinurStatements
} // namespace PsimagLite
#endif // AINURSTATEMENT_H
