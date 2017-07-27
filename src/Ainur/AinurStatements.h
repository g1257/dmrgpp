#ifndef AINURSTATEMENT_H
#define AINURSTATEMENT_H
#include "AinurLexical.h"
#include "AinurStore.h"
#include <complex>

namespace PsimagLite {

class AinurStatements {

public:

	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef AinurLexical AinurLexicalType;
	typedef AinurLexicalType::VectorStringType VectorStringType;
	typedef Vector<Store>::Type VectorStoreType;

	AinurStatements(const VectorStringType& vecStr,
	                const String& vecChar,
	                const String& escapedChars,
	                const VectorStringType& vecBrace)
	    : vecStr_(vecStr),
	      vecChar_(vecChar),
	      escapedChars_(escapedChars),
	      vecBrace_(vecBrace),
	      prefix_("")
	{}

	VectorStringType push(const String& s2, String prefix)
	{
		VectorStringType emptyStringVector;
		String s = s2;
		AinurLexicalType::removeTrailingWhitespace(s);
		AinurLexicalType::removeTrailingWhitespace(prefix);
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
		Store& store = storage_[storageIndex];
		if (store.type() == Store::SCALAR && store.subType() == Store::GROUP)  {
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
			Store& store = storage_[i];
			SizeType total = store.valueSize();

			for (SizeType j = 0; j < total; ++j) {
				String replacement = solveExpression(store.value(j, name), i);
				store.value(j, name) = replacement;
			}
		}
	}

	String& prefix() { return prefix_; }

	const String& prefix() const { return prefix_; }

	void printUnused(std::ostream& os) const
	{
		SizeType n = storage_.size();
		for (SizeType i = 0; i < n; ++i) {
			if (storage_[i].used() > 0)
				continue;
			assert(i < names_.size());
			os<<"WARNING: Unused label "<<names_[i]<<"\n";
		}
	}

	void readValue(long int& t, String s) const
	{
		int t2 = 0;
		readValue(t2, s);
		t = t2;
	}

	void readValue(SizeType& t, String s) const
	{
		int t2 = 0;
		readValue(t2, s);
		t = t2;
	}

	void readValue(int& t, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::INTEGER)
			err("In input, " + s + " must be an integer\n");
		store.increaseUsage();
		t = atoi(store.value(0, names_[x]).c_str());
	}

	void readValue(RealType& t, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::REAL)
			err("In input, " + s + " must be a real\n");
		store.increaseUsage();
		t = atof(store.value(0, names_[x]).c_str());
	}

	void readValue(String& t, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::SCALAR && store.subType() != Store::STRING)
			err("In input, " + s + " must be a string\n");
		store.increaseUsage();
		t = store.value(0, names_[x]);
	}

	// read vectors
	template<typename VectorLikeType>
	typename EnableIf<IsVectorLike<VectorLikeType>::True,void>::Type
	readValue(VectorLikeType& v, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::VECTOR)
			err("In input, " + s + " must be a vector\n");
		store.increaseUsage();
		SizeType n = store.valueSize();
		if (n == 0)
			err("In input, vector " + s + " has 0 entries\n");

		store.increaseUsage();

		String tmp = (n == 2) ? store.value(1, names_[x]) : "";
		AinurLexicalType::removeTrailingWhitespace(tmp);
		if (n == 2 && tmp == "...") {
			assert(static_cast<SizeType>(x) < names_.size());
			if (v.size() < 3)
				err("Ellipsis cannot be used for this vector, " + names_[x] + "\n");
			n = v.size();
			for (SizeType i = 0; i < n; ++i)
				getEntryFromString(v[i], store.value(0, names_[x]));
			return;
		}

		if (v.size() != n) {
			v.clear();
			v.resize(n);
		}

		for (SizeType i = 0; i < n; ++i)
			getEntryFromString(v[i], store.value(i, names_[x]));
	}

	// read matrices
	template<typename FloatingType>
	typename EnableIf<Loki::TypeTraits<FloatingType>::isFloat,void>::Type
	readValue(Matrix<FloatingType>& m, String s) const
	{
		s = prefix_ + s;
		int x = storageIndexByName(s);
		if (x < 0)
			err("Not found " + s + "\n");
		const Store& store = storage_[x];
		if (store.type() != Store::MATRIX)
			err("In input, " + s + " must be a matrix\n");
		store.increaseUsage();
		SizeType n = store.valueSize();
		if (n < 2)
			err("In input, matrix " + s + " internal storage error I\n");

		store.increaseUsage();

		SizeType rows = atoi(store.value(0, names_[x]).c_str());
		SizeType cols = atoi(store.value(1, names_[x]).c_str());

		m.clear();
		if (rows == 0 && cols == 0) return;
		if (rows*cols == 0)
			err("Matrix with one of {rows, cols} 0 must have both 0\n");

		m.resize(rows, cols);

		if (rows*cols +2 != n)
			err("In input, matrix " + s + " internal storage error II\n");

		for (SizeType i = 0; i < rows; ++i)
			for (SizeType j = 0; j < cols; ++j)
				getEntryFromString(m(i,j), store.value(i + j*rows + 2, names_[x]));
	}

private:

	void getEntryFromString(SizeType& entry, String s) const
	{
		entry = atoi(s.c_str());
	}

	void getEntryFromString(RealType& entry, String s) const
	{
		entry = atof(s.c_str());
	}

	// One of double or double + double*i or double - double*i
	// will later replace with general expression evaluator
	void getEntryFromString(ComplexType& entry, String s) const
	{
		err("getEntryFromString not implemented for complex\n");
	}

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
			x = storageIndexByName(dotified[0]);
		} else if (l == 2) { // matrix.integer FiniteLoops
			split(dotified,prefix + lhs[1],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = assignStorageByName(dotified[0]);
			storage_.push_back(Store(lhs[0],""));
		} else if (l == 3) {
			// require vector.vector.integer FiniteLoops
			split(dotified,prefix + lhs[2],".");
			if (dotified.size() != 1)
				err("Dotified failed " + context + "\n");
			x = assignStorageByName(dotified[0]);
			storage_.push_back(Store(lhs[1], lhs[0]));
		}

		if (x < 0)
			err("Undeclared variable " + dotified[0] + "\n");

		y = x;

		return dotified[0];
	}

	Store::Attribute getAttribute(String s, String context) const
	{
		if (s == "require") return Store::REQUIRED;
		if (s == "const") return Store::CONST;

		err("Expected let require or const " + context + "\n");
		return Store::NONE;
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
		VectorStringType::const_iterator it = std::find(names_.begin(),
		                                                names_.end(),
		                                                name);
		if (it == names_.end())
			return -1;
		return it - names_.begin();
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
			const Store& store = storage_[i];
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
	String prefix_;
}; // class AinurStatements
} // namespace PsimagLite
#endif // AINURSTATEMENT_H
